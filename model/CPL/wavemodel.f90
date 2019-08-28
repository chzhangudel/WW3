module wave_model_mod

  use mpp_mod, only: mpp_npes, mpp_pe, mpp_get_current_pelist
  use mpp_domains_mod, only: mpp_global_field
  use mpp_domains_mod, only: mpp_define_domains, mpp_define_layout

  use wave_type_mod, only: wave_data_type, atmos_wave_boundary_type

  use time_manager_mod, only: time_type, operator(+), get_date

  USE WMMDATMD, ONLY: MDSI, MDSO, MDSS, MDST, MDSE, &
                      NMPROC, IMPROC, NMPSCR, NRGRD, ETIME


  implicit none

  INTEGER, ALLOCATABLE :: TEND(:,:), TSTRT(:,:)
  INTEGER              :: MPI_COMM = -99

  contains


    !> wave_model_init initializes the wave model interface
    subroutine wave_model_init(Cplr2Waves,Waves)
      ! Author: Brandon Reichl
      ! Origination: August 2019

      ! USE statements
      use wminitmd, only: wminit, wminitnml
      use w3gdatmd, only: NX, NY


      ! Subroutine arguments
      type(atmos_wave_boundary_type), intent(inout)  :: Cplr2Waves
      type(wave_data_type), intent(inout)            :: Waves

      ! Local parameters
      logical              :: Flag_NML
      INTEGER              :: Ierr_MPI
      integer :: layout(2)

      !----------------------------------------------------------------------
      !This block of code checks if the wave model has been initialized
      if (Waves%Waves_Is_Init) then
        write(*,*)'wave_model_init was called, but it is already registered'
        stop
      end if
      Waves%Waves_Is_Init = .true.
      !----------------------------------------------------------------------

      !----------------------------------------------------------------------
      !This block of code sets up the MPI communicator for WW3
      ! The equivalent calls in ww3_multi are given in comments
      ! -> FMS has already done the MPI_INIT.
      ! CALL MPI_INIT      ( IERR_MPI )
      ! MPI_COMM = MPI_COMM_WORLD
      ! -> We need to access the commID from FMS.  This is the only
      !    way I can figure out using MPP libraries.
      call mpp_get_current_pelist(Waves%pelist,commID=MPI_COMM)
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
      call mpp_define_domains((/1,NX,1,NY/),layout, Waves%domain )
      !----------------------------------------------------------------------

      !
      !This block of code sets up the Cplr2Wave domains
      ! this is done in atm_wave_exchange_init
      allocate(Cplr2Waves%U_10_global(1:NX,1:NY,1)) ; Cplr2Waves%U_10_global(1:NX,1:NY,1) = 0.0
      allocate(Cplr2Waves%V_10_global(1:NX,1:NY,1)) ; Cplr2Waves%V_10_global(1:NX,1:NY,1) = 0.0

      return
    end subroutine wave_model_init

    !> wave_model_timestep integrates the wave fields over one
    !!  coupling timestep
    subroutine update_wave_model(Cplr2Waves,Waves,Time_start,Time_increment)
      ! Author: Brandon Reichl
      ! Origination: August 2019

      ! USE statements
      use wmwavemd, only: wmwave
      use w3gdatmd, only: NX, NY
      use w3idatmd, only: wx0, wxN, wy0, wyN, TW0, TWN
      ! Subroutine arguments
      type(atmos_wave_boundary_type), intent(in) :: Cplr2Waves
      type(wave_data_type), intent(inout) :: Waves
      type(time_type), intent(in)         :: Time_start, Time_increment

      ! Local parameters
      integer :: I, yr, mo, da, hr, mi, se

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
      TW0 = TSTRT(:,1)
      TWN = TEND(:,1)
      call mpp_global_field(Waves%domain,cplr2waves%u_10_mpp(:,:,:),cplr2waves%u_10_global(:,:,:))
      wx0 = cplr2waves%u_10_global(:,:,1)
      wxN = wx0
      call mpp_global_field(Waves%domain,cplr2waves%v_10_mpp(:,:,:),cplr2waves%v_10_global(:,:,:))
      wy0 = cplr2waves%v_10_global(:,:,1)
      wyN = wy0
      write(*,*)'Into wave model U10 min: ',minval(Cplr2Waves%U_10_global),minval(wx0)
      write(*,*)'Into wave model U10 max: ',maxval(Cplr2Waves%U_10_global),maxval(wx0)
      write(*,*)'Into wave model V10 min: ',minval(Cplr2Waves%V_10_global),minval(wy0)
      write(*,*)'Into wave model V10 max: ',maxval(Cplr2Waves%V_10_global),maxval(wy0)
      CALL WMWAVE ( TEND )
      !----------------------------------------------------------------------

      return
    end subroutine update_wave_model

    !>This subroutine finishes the wave routine and deallocates memory
    subroutine wave_model_end(Waves)
      ! Use statements
      use wmfinlmd, only: wmfinl

      ! Subroutine arguments
      type(wave_data_type), intent(inout) :: Waves

      ! Local parameters
      INTEGER              :: IERR_MPI
      !----------------------------------------------------------------------
      !Finalize the driver
      CALL WMFINL
      DEALLOCATE ( TEND, TSTRT )
      CALL MPI_BARRIER ( MPI_COMM, IERR_MPI ) !Do we need this?
      !----------------------------------------------------------------------

      return
    end subroutine wave_model_end


end module wave_model_mod
