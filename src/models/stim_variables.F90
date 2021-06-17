!-----------------------------------------------------------------------
!BOP
!
! !MODULE: stim_variables - global variables for the ice models
!
! !INTERFACE:
   module stim_variables
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !PUBLIC DATA MEMBERS:
   integer, parameter, public :: rk = kind(1.d0)
   real(rk), public :: ice_cord = 1._rk
!KB   integer, parameter, public :: rk = kind(1.0)
   real(rk), public :: Cw = 4.18e+6_rk    ! volumetric heat capacity of water (J K-1 m-3)
   real(rk), public :: L_ice = 333500._rk ! latent heat of freezing (J kg-1)
   real(rk), public :: K_ice = 2.1_rk     ! ice heat conduction coefficient (W m-1 K-1)
   real(rk), public :: lambda_ice = 5._rk ! 
   real(rk), public :: rho_ice = 910._rk  ! kg m-3

   real(rk), public :: init_ice_energy = 0._rk    ! 
   real(rk), target, public :: surface_ice_energy = 0._rk    ! 
   real(rk), target, public :: bottom_ice_energy = 0._rk    ! 

   real(rk), target, public :: Tf = 0._rk ! freezing temperature
   real(rk), target, public :: Tice_surface = 0._rk
   real(rk), allocatable, target, public :: Tice(:)
   real(rk), target, public :: Hfrazil = 0._rk ! Total frazil ice thickness
   real(rk), target, public :: Hsnow = 0._rk   ! Total snow thickness
   real(rk), target, public :: Hice = 0._rk    ! Total ice thickness
   real(rk), target, public :: dHis = 0._rk    ! surface ice growth
   real(rk), target, public :: dHib = 0._rk    ! bottom ice growth
   real(rk), target, public :: albedo_ice = 0._rk
   real(rk), target, public :: attenuation_ice = 0._rk
   real(rk), target, public :: ocean_ice_flux = 0._rk
   real(rk), target, public :: transmissivity = 1._rk

   ! Lebedev
   real(rk), target, public :: fdd

   !Flato 
! Flato PUBLIC DATA MEMBERS: 
! parameters
   public :: hlaymin,rhoice,Tfreezi,rCpmix,Hfi

! !DEFINED PARAMETERS:
!   hsmin               - minimum snow thickness required to have a separate snow
!                       layer. If 'ice_hs' is less than this snow is ignored in
!                       heat conduction scheme (m)
   real(rk), parameter :: hsmin = 0.001D+00  
!   theta               - a parameter between 0.5 and 1. which determines how  
!                       'implicit' the scheme is. theta=0.5 is the conventional
!                       Crank-Nicholson algorithm and theta=1. is the fully
!                       implicit scheme (the method of choice in this model because
!		                  C-R scheme allows spurious temperature oscillations when
!		                  layer thickness becomes small - an instability in the
!		                  numerical scheme which is avoided by using theta=1).
   real(rk), parameter :: theta = 1.D+00 
!   sigma            - Stefan-Boltzmann constant (W m-2 K-4)
   real(rk), parameter :: sigma = 5.67D-08  
!   epsilon          - emissivity of ice (dimensionless)
   real(rk), parameter :: epsilon = 0.99D+00
!   PenFrac          - fraction of incoming short wave radiation that penetrates  the surface
   real(rk), parameter :: PenFrac = 0.17D+00 
!   hlaymin          - thickness below which a linear temperature profile is assumed (m)
   real(rk), parameter :: hlaymin = 0.2D+00
!   rhoscold         - specified 'cold' snow density (kg m-3)
   real(rk), parameter :: rhoscold = 330.D+00
!   rhoswarm         - specified 'warm' snow density (kg m-3)
   real(rk), parameter :: rhoswarm = 450.D+00
!   rhowaterfresh    - fresh water density (kg m-3)
   real(rk), parameter :: rhowaterfresh = 1000.D+00
!   rhoice           - ice density (kg m-3)
   real(rk), parameter :: rhoice = 913.D+00
!************
!   kelvin           - zero deg Celsius (K)
   real(rk), parameter :: kelvin = 273.15D+00
!************   !kelvin exists in air_sea variables.f90 --> as 273.16 
!   Tmelts            - melting temperature of snow (fresh water) (K)
   real(rk), parameter :: Tmelts = 273.05D+00
!   Tmelti      - melting temperature of sea ice (K)
   real(rk), parameter  :: Tmelti=273.05D+00
!   Condfi            - conductivity of pure ice (W m-1 K-1)
   real(rk), parameter :: Condfi = 2.034D+00
!   rhoCpfi           - heat capacity of pure ice (J m-3 K-1)
   real(rk), parameter :: rhoCpfi = 1.883D+06
!   rCpmix            - volumetric heat capacity of sea water (J m-3 K-1)
   real(rk), parameter :: rCpmix = 1025.D+00*3.99D+03
!   Hfi               - latent heat of fusion of sea ice (J kg-1)
   real(rk), parameter :: Hfi = 2.93D+05
!   Hfw               - latent heat of fusion of fresh water (J kg-1)
   real(rk), parameter :: Hfw = 3.335D+05
!   swkappa           - bulk short-wave extinction coefficient (m-1)
   real(rk), parameter :: swkappa = 1.5
!   Tfreezi           - freezing temperature of sea water (K)
   real(rk), parameter :: Tfreezi = 271.2D+00

! !Maximum snow and ice layers
   integer, public, parameter :: nlmax=99

!-----------------------------------------------------------

!  !initializing yaml variables to reasonable defaults 
!   nilay              - number of snow & ice layers
   integer, public :: nilay = 0             
!   sfall_method       - define how snow fall is determined 
!                        1:constant snow fall
!                        2:calculate snowfall from precipitation 
   integer, public :: sfall_method = 0       
!   const_sfall        - constant snow fall rate (m d^-1)
   real(rk), public :: const_sfall = 0._rk 
!   dfact              - drift factor allowing a factor to increase snow fall from 
!                        precipitation via drifing snow (only for sfall_method=2) 
   real(rk), public :: dfact = 0._rk
!   depmix             - prescribed mixed layer depth for open_water
!                         calculation. If set to zero top layer thickness is used.   
   real(rk), public :: depmix = 0._rk
!   sice_method        - define how sea-ice salinity is to be calculated
!                        1: constant ice salinity
!                        2: Simple ice salinity profile (Vancoppenolle et al. 2009?)
   integer, public :: sice_method = 0
!   snow_dist          - logical switch between uniform and Weibull-distributed snow
!                        if true a snow distribution function will be applied
!                        if false a uniform snow thickness will be applied  
   logical, public :: snow_dist = .false.
!   const_Sice         - prescribed sea ice salinity (ppt)
   real(rk), public :: const_Sice = 0._rk
!   distr_type         - integer to chose the type of distribution 
   integer, public :: distr_type = -1
!   meltpond           - logical to switch on meltpond
   logical, public :: meltpond = .false.
!   Ameltmax           - maximum meltpond area
   real(rk), public :: Ameltmax = 0._rk
!   drainrate          - fixed meltpond drainrate - (m/d) converted to (kg/m^2/s)
   real(rk), public :: drainrate = 0._rk
!   hh0                -  initial thickness for S calculation        
   real(rk), public :: hh0 = 0._rk
!  ice_hi_i            - initial ice thickness (m) 
   real(rk), public :: ice_hi_i = 0.2D+00
!  ice_hs_i            - initial snow thickness (m) 
   real(rk), public :: ice_hs_i = 0.01D+00
!  albice_method       - 1:constant albedo (albice_f and albsnow_f)
!                        2:albice dependent on ice thickness - attempt to
!                          capture the evolution of melt ponds (Flato&Brown 1996)
!                          albsnow is set to albsnow_f (constant)
!                        3:albice and albsnow based on eq12&13 of Flato&Brown1996
!                        4:albice and albsnow dependent on temperature (same as transs/transi formulation)
   integer, public :: albice_method = 0
!  albice_f            - freezing ice albedo 
   real(rk), public :: albice_f = 0._rk
!  albmelt             - melt pond albedo   
   real(rk), public :: albmelt = 0._rk
!  albsnow_f           - freezing snow albedo
   real(rk), public :: albsnow_f = 0._rk
!  albice_m            - melting ice albedo
   real(rk), public :: albice_m = 0._rk
!  albsnow_m           - melting snow albedo
   real(rk), public :: albsnow_m = 0._rk
!  transsf             - freezing snow transmission coefficient
   real(rk), public :: transsf = 0._rk
!  transsm             - melting snow transmission coefficient
   real(rk), public :: transsm = 0._rk
!  transif             - freezing ice transmission coefficient
   real(rk), public :: transif = 0._rk
!  transim             - melting ice transmission coefficient
   real(rk), public :: transim = 0._rk
!  transm              - melt pond transmision coefficient
   real(rk), public :: transm = 0._rk
!  swkappasm           - melting snow extinction coefficient
   real(rk), public :: swkappasm = 0._rk
!  swkappasf           - freezing snow extinction coefficient
   real(rk), public :: swkappasf = 0._rk
!  swkappaim           - melting ice extinction coefficient
   real(rk), public :: swkappaim = 0._rk
!  swkappaif           - freezing ice extinction coefficient
   real(rk), public :: swkappaif = 0._rk

!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_stim_variables - initialise 2D related stuff.
!
! !INTERFACE:
   subroutine init_stim_variables(ice_model)
   IMPLICIT NONE
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
   integer, intent(in)     :: ice_model
!
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_stim_variables() # ',Ncall
#endif

   !LEVEL2 'init_stim_variables'
   !LEVEL3 'using ice model ',ice_model

#ifdef DEBUG
   write(debug,*) 'Leaving init_stim_variables()'
   write(debug,*)
#endif
   return
   end subroutine init_stim_variables
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: clean_stim_variables - cleanup after
!
! !INTERFACE:
   subroutine clean_stim_variables()
   IMPLICIT NONE
!
! !DESCRIPTION:
!  This routine is currently empty.
!
! !LOCAL VARIABLES:
!
!EOP
!-----------------------------------------------------------------------
!BOC

   return
   end subroutine clean_stim_variables
!EOC

!-----------------------------------------------------------------------

   end module stim_variables

!-----------------------------------------------------------------------
! Copyright (C) 2019 - Karsten Bolding (BB)                            !
!-----------------------------------------------------------------------
