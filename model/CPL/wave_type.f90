module wave_type_mod

use mpp_domains_mod,  only: domain2D
use MOM_time_manager,  only : time_type

public :: wave_data_type, atmos_wave_boundary_type, ice_wave_boundary_type

private

type wave_data_type
  type(domain2D)          :: Domain
  type(time_type)         :: Time
  logical                 :: pe
  logical                 :: Waves_Is_Init = .false.
  integer                 :: xtype
  integer                 :: num_Stk_bands
  integer, pointer, dimension(:)   :: pelist   =>NULL() !< Used for flux-exchange.
  logical, pointer, dimension(:,:) :: ocean_pt =>NULL() !< An array that indicates ocean points as true.

  real, pointer, dimension(:) :: stk_wavenumbers => NULL() !< Used for stokes drift

  ! These fields are used to provide information about the waves to the atmosphere/ocean/ice
  real, pointer, dimension(:,:) :: &
       HS => NULL() !< The significant wave height [m]

  real, pointer, dimension(:,:,:) :: &
       ustkb_mpp => NULL(), &
       vstkb_mpp => NULL(), &
       ustkb_glo => NULL(), &
       vstkb_glo => NULL()

  integer, dimension(:,:), pointer :: & ! (lon, lat,tile)
       glob_loc_X => NULL(), & !
       glob_loc_Y => NULL()

  ! These fields provide information from the atmosphere/ocean/ice to the waves
  real, pointer, dimension(:,:) :: &
       U10 => NULL(), &
       V10 => NULL(), &
       ICE_CONCENTRATION => NULL()

end type wave_data_type

type atmos_wave_boundary_type
   real, dimension(:,:,:), pointer :: & ! (lon, lat,tile)
        wavgrd_u10_mpp => NULL(), & !
        wavgrd_v10_mpp => NULL()
   real, dimension(:,:,:), pointer :: & ! (lon, lat,tile)
        wavgrd_u10_glo => NULL(), & !
        wavgrd_v10_glo => NULL()

   integer :: xtype             !REGRID, REDIST or DIRECT

end type atmos_wave_boundary_type

type ice_wave_boundary_type

   real, dimension(:,:,:), pointer :: & ! (lon, lat,tile)
        wavgrd_Ucurr_mpp => NULL(), & !
        wavgrd_Vcurr_mpp => NULL(), & !
        wavgrd_ucurr_glo => NULL(), & !
        wavgrd_vcurr_glo => NULL()

   real, dimension(:,:,:,:), pointer :: & ! (lon, lat,tile,Nstk)
        icegrd_ustkb_mpp => NULL(), & !
        icegrd_vstkb_mpp => NULL()

   integer :: xtype             !REGRID, REDIST or DIRECT

end type ice_wave_boundary_type

end module wave_type_mod
