module wave_model_mod

  use mpp_mod, only: mpp_npes, mpp_pe, mpp_get_current_pelist
  use mpp_domains_mod, only: mpp_global_field, mpp_get_compute_domain
  use mpp_domains_mod, only: mpp_define_domains, mpp_define_layout

  use wave_type_mod, only: wave_data_type, atmos_wave_boundary_type, ice_wave_boundary_type

  use time_manager_mod, only: time_type, operator(+), get_date

  USE WMMDATMD, ONLY: MDSI, MDSO, MDSS, MDST, MDSE, &
                      NMPROC, IMPROC, NMPSCR, NRGRD, ETIME


  implicit none

  INTEGER, ALLOCATABLE :: TEND(:,:), TSTRT(:,:)
  INTEGER              :: MPI_COMM = -99

  contains


    !> wave_model_init initializes the wave model interface
    subroutine wave_model_init(Atm2Waves,Ice2Waves,Wav)
      ! Author: Brandon Reichl
      ! Origination: August 2019

      ! USE statements
      use wminitmd, only: wminit, wminitnml
      use w3gdatmd, only: NX, NY
      use w3gdatmd, only: usspf, ussp_wn

      ! Subroutine arguments
      type(atmos_wave_boundary_type), intent(inout)  :: Atm2Waves
      type(ice_wave_boundary_type), intent(inout)  :: Ice2Waves
      type(wave_data_type), intent(inout)            :: Wav

      ! Local parameters
      logical              :: Flag_NML
      INTEGER              :: Ierr_MPI
      integer :: layout(2)

      !----------------------------------------------------------------------
      !This block of code checks if the wave model has been initialized
      if (Wav%Waves_Is_Init) then
        write(*,*)'wave_model_init was called, but it is already registered'
        stop
      end if
      Wav%Waves_Is_Init = .true.
      !----------------------------------------------------------------------

      !----------------------------------------------------------------------
      !This block of code sets up the MPI communicator for WW3
      ! The equivalent calls in ww3_multi are given in comments
      ! -> FMS has already done the MPI_INIT.
      ! CALL MPI_INIT      ( IERR_MPI )
      ! MPI_COMM = MPI_COMM_WORLD
      ! -> We need to access the commID from FMS.  This is the only
      !    way I can figure out using MPP libraries.
      call mpp_get_current_pelist(Wav%pelist,commID=MPI_COMM)
      ! -> FMS does provide an easy way to get the comm_size and rank
      ! CALL MPI_COMM_SIZE ( MPI_COMM, NMPROC, IERR_MPI )
      !CALL MPI_COMM_RANK ( MPI_COMM, IMPROC, IERR_MPI )
      NMPROC = mpp_npes()
      IMPROC = mpp_pe()
      IMPROC = IMPROC + 1
      !----------------------------------------------------------------------

      !----------------------------------------------------------------------
      !This block of code calls the WW3 intializers
      INQUIRE(FILE="ww3_multi.nml", EXIST=Flag_NML)
      IF (Flag_NML) THEN
        CALL WMINITNML ( MDSI, MDSO, MDSS, 10, MDSE, 'ww3_multi.nml', MPI_COMM )
      ELSE
        CALL WMINIT ( MDSI, MDSO, MDSS, 10, MDSE, 'ww3_multi.inp', MPI_COMM )
      END IF
      ALLOCATE ( TEND(2,NRGRD), TSTRT(2,NRGRD) )
      !----------------------------------------------------------------------

      !----------------------------------------------------------------------
      !This block of codes sets up the wave domains
      call mpp_define_layout((/1,NX,1,NY/),mpp_npes(),layout)
      call mpp_define_domains((/1,NX,1,NY/),layout, Wav%domain )
      !----------------------------------------------------------------------

      !
      !This block of code sets up the global coupler domains
      allocate(Atm2Waves%wavgrd_u10_glo(1:NX,1:NY,1))
      Atm2Waves%wavgrd_u10_glo(1:NX,1:NY,1) = 0.0
      allocate(Atm2Waves%wavgrd_v10_glo(1:NX,1:NY,1))
      Atm2Waves%wavgrd_v10_glo(1:NX,1:NY,1) = 0.0
      allocate(Ice2Waves%wavgrd_ucurr_glo(1:NX,1:NY,1))
      Ice2Waves%wavgrd_ucurr_glo(1:NX,1:NY,1) = 0.0
      allocate(Ice2Waves%wavgrd_vcurr_glo(1:NX,1:NY,1))
      Ice2Waves%wavgrd_vcurr_glo(1:NX,1:NY,1) = 0.0

      if (usspf(1).gt.usspf(2)) then
        print*,'usspf(1): ',usspf(1)
        print*,'usspf(2): ',usspf(2)
        print*,'usspf(1) is greater than usspf(2) in mod_def.ww3. ',&
             ' Check ww3_grid.inp file, recreate mod_def.ww3 with ussp <= iussp to continue'
        stop
      endif
      Wav%num_Stk_bands = usspf(2)-usspf(1)+1
      allocate(Wav%stk_wavenumbers(1:Wav%num_Stk_Bands)) ; Wav%stk_wavenumbers(:) = 0.0
      Wav%stk_wavenumbers(:) = USSP_WN(usspf(1):usspf(2))

      allocate(Wav%ustkb_glo(1:NX,1:NY,Wav%num_stk_bands))
      Wav%ustkb_glo(:,:,:) = 0.0
      allocate(Wav%vstkb_glo(1:NX,1:NY,Wav%num_stk_bands))
      Wav%vstkb_glo(:,:,:) = 0.0

      return
    end subroutine wave_model_init

    !> wave_model_timestep integrates the wave fields over one
    !!  coupling timestep
    subroutine update_wave_model(Atm2Waves, Ice2Waves, Wav, Time_start, Time_increment)
      ! Author: Brandon Reichl
      ! Origination: August 2019

      ! USE statements
      use wmwavemd, only: wmwave
      use w3gdatmd, only: NX, NY
      use w3idatmd, only: wx0, wxN, wy0, wyN, TW0, TWN, &
                          cx0, cxN, cy0, cyN, TC0, TCN
      use w3adatmd, only: ussx, ussy, ussp
      use w3gdatmd, only: nseal, mapsf, NK
      use w3odatmd, only: iaproc, naproc
      ! Subroutine arguments
      type(atmos_wave_boundary_type), intent(in) :: Atm2Waves
      type(ice_wave_boundary_type),   intent(in) :: Ice2Waves
      type(wave_data_type),        intent(inout) :: Wav
      type(time_type),                intent(in) :: Time_start,&
                                                    Time_increment

      ! Local parameters
      integer :: I, yr, mo, da, hr, mi, se, is, ie, js, je
      integer :: isea, isea_g, ix, iy, jx, jy, pix, piy
      integer :: b

      integer :: glob_loc_x(NX,NY), glob_loc_y(NX,NY)

      !----------------------------------------------------------------------
      !Convert the ending time of this call into WW3 time format, which
      ! is integer(2) :: (YYYYMMDD, HHMMSS)
      DO I=1, NRGRD
        call get_date(Time_start,&
             yr,mo,da,hr,mi,se)
        TSTRT(1,I) = yr*1e4+mo*1e2+da
        TSTRT(2,I) = hr*1e4+mi*1e2+se
        call get_date(Time_start+Time_increment,&
             yr,mo,da,hr,mi,se)
        TEND(1,I) = yr*1e4+mo*1e2+da
        TEND(2,I) = hr*1e4+mi*1e2+se
      END DO
      !----------------------------------------------------------------------

      !----------------------------------------------------------------------
      !Call WW3 timestepper with an argument for the time to return
      ! back
      TW0(:) = TSTRT(:,1)
      TWN(:) = TEND(:,1)
      TC0(:) = TSTRT(:,1)
      TCN(:) = TEND(:,1)
      call mpp_global_field(Wav%domain,atm2waves%wavgrd_u10_mpp(:,:,:),atm2waves%wavgrd_u10_glo(:,:,:))
      wx0(:,:) = atm2waves%wavgrd_u10_glo(:,:,1)
      wxN(:,:) = wx0(:,:)
      call mpp_global_field(Wav%domain,atm2waves%wavgrd_v10_mpp(:,:,:),atm2waves%wavgrd_v10_glo(:,:,:))
      wy0(:,:) = atm2waves%wavgrd_v10_glo(:,:,1)
      wyN(:,:) = wy0(:,:)
      call mpp_global_field(Wav%domain,ice2waves%wavgrd_ucurr_mpp(:,:,:),ice2waves%wavgrd_ucurr_glo(:,:,:))
      cx0(:,:) = ice2waves%wavgrd_ucurr_glo(:,:,1)
      cxN(:,:) = cx0(:,:)
      call mpp_global_field(Wav%domain,ice2waves%wavgrd_vcurr_mpp(:,:,:),ice2waves%wavgrd_vcurr_glo(:,:,:))
      cy0(:,:) = ice2waves%wavgrd_vcurr_glo(:,:,1)
      cyN(:,:) = cy0(:,:)
      ! write(*,*)'Into wave model U10 max: ',maxval(atm2Waves%U_10_global),maxval(wx0)
      ! write(*,*)'Into wave model V10 max: ',maxval(atm2Waves%V_10_global),maxval(wy0)
      ! write(*,*)'Into wave model UO max: ',maxval(ice2Waves%Ucurr_global),maxval(cx0)
      ! write(*,*)'Into wave model VO max: ',maxval(ice2Waves%Vcurr_global),maxval(cy0)
      CALL WMWAVE ( TEND )

      Wav%glob_loc_X(:,:) = 0
      Wav%glob_loc_Y(:,:) = 0

      call mpp_get_compute_domain( Wav%domain, is, ie, js, je )
      isea = 1
      do ix=is,ie
        do iy=js,je
          if (isea<nseal) then
            ISEA_G   = IAPROC + (ISEA-1)*NAPROC
            jx = MAPSF(ISEA_G,1)
            jy = MAPSF(ISEA_G,2)
            Wav%glob_loc_X(ix,iy) = jx
            Wav%glob_loc_Y(ix,iy) = jy
            do b=1,Wav%num_stk_bands
              Wav%ustkb_mpp(ix,iy,b) = USSP(isea,b)
              Wav%vstkb_mpp(ix,iy,b) = USSP(isea,NK+b)
            enddo
          endif
          isea = isea+1
        end do
      end do

      do b = 1,Wav%num_stk_bands
        call mpp_global_field(Wav%domain,Wav%ustkb_mpp(:,:,b),wav%ustkb_glo(:,:,b))
        call mpp_global_field(Wav%domain,Wav%vstkb_mpp(:,:,b),wav%vstkb_glo(:,:,b))
      enddo
      call mpp_global_field(Wav%domain,Wav%glob_loc_X,glob_loc_X)
      call mpp_global_field(Wav%domain,Wav%glob_loc_Y,glob_loc_Y)

      Wav%ustkb_mpp(:,:,:) = 0.0
      Wav%vstkb_mpp(:,:,:) = 0.0

      isea = 1
      do ix=1,NX
         do iy=1,NY
           Pix = glob_loc_X(ix,iy)
           Piy = glob_loc_Y(ix,iy)
           if (Pix>=is .and. Pix<=ie .and. Piy>=js .and. Piy<=je) then
              do b = 1,Wav%num_stk_bands
                Wav%ustkb_mpp(Pix,Piy,b) = wav%ustkb_glo(ix,iy,b)
                Wav%vstkb_mpp(Pix,Piy,b) = wav%vstkb_glo(ix,iy,b)
              enddo
           endif
         enddo
      enddo

      !----------------------------------------------------------------------

      return
    end subroutine update_wave_model

    !>This subroutine finishes the wave routine and deallocates memory
    subroutine wave_model_end(Waves, Atm2Waves,Ice2Waves)
      ! Use statements
      use wmfinlmd, only: wmfinl

      ! Subroutine arguments
      type(wave_data_type), intent(inout) :: Waves
      type(atmos_wave_boundary_type), intent(inout) :: Atm2Waves
      type(ice_wave_boundary_type), intent(inout) :: Ice2Waves

      ! Local parameters
      INTEGER              :: IERR_MPI
      !----------------------------------------------------------------------
      !Finalize the driver
      CALL WMFINL
      DEALLOCATE ( TEND, TSTRT )
      deallocate(Atm2Waves%wavgrd_u10_glo)
      deallocate(Atm2Waves%wavgrd_v10_glo)
      deallocate(Ice2Waves%wavgrd_ucurr_glo)
      deallocate(Ice2Waves%wavgrd_vcurr_glo)
      CALL MPI_BARRIER ( MPI_COMM, IERR_MPI ) !Do we need this?
      !----------------------------------------------------------------------

      return
    end subroutine wave_model_end


end module wave_model_mod
