module wave_type_mod

use mpp_domains_mod,  only: domain2D
use MOM_time_manager,  only : time_type

public :: wave_data_type

type wave_data_type
  type(domain2D)          :: Domain
  type(time_type)         :: Time
  logical                 :: pe
  logical                 :: Waves_Is_Init = .false.
  integer                 :: xtype
  integer, pointer, dimension(:)   :: pelist   =>NULL() !< Used for flux-exchange.
  logical, pointer, dimension(:,:) :: ocean_pt =>NULL() !< An array that indicates ocean points as true.

  ! These fields are used to provide information about the waves to the atmosphere/ocean/ice
  real, pointer, dimension(:,:) :: &
       HS => NULL() !< The significant wave height [m]


  ! These fields provide information from the atmosphere/ocean/ice to the waves
  real, pointer, dimension(:,:) :: &
       U10 => NULL(), &
       V10 => NULL(), &
       ICE_CONCENTRATION => NULL()

end type wave_data_type

type atmos_wave_boundary_type
   real, dimension(:,:,:), pointer :: & ! (lon, lat,tile)
        U_10_mpp => NULL(), & !
        V_10_mpp => NULL()
   real, dimension(:,:,:), pointer :: & ! (lon, lat,tile)
        U_10_global => NULL(), & !
        V_10_global => NULL()

   integer :: xtype             !REGRID, REDIST or DIRECT

end type atmos_wave_boundary_type

end module wave_type_mod
