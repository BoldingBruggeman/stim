!-----------------------------------------------------------------------
!BOP
!
! !MODULE: ice_variables - global variables for the ice models
!
! !INTERFACE:
   module ice_variables
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
   real(rk), public :: melt_ice_energy = 0._rk    ! 

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
   real(rk), target, public :: sensible_ice_water = 0._rk

   ! Lebedev
   real(rk), target, public :: fdd
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
! !IROUTINE: init_ice_variables - initialise 2D related stuff.
!
! !INTERFACE:
   subroutine init_ice_variables(ice_model)
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
   write(debug,*) 'init_ice_variables() # ',Ncall
#endif

   !LEVEL2 'init_ice_variables'
   !LEVEL3 'using ice model ',ice_model

#ifdef DEBUG
   write(debug,*) 'Leaving init_ice_variables()'
   write(debug,*)
#endif
   return
   end subroutine init_ice_variables
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: clean_ice_variables - cleanup after
!
! !INTERFACE:
   subroutine clean_ice_variables()
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
   end subroutine clean_ice_variables
!EOC

!-----------------------------------------------------------------------

   end module ice_variables

!-----------------------------------------------------------------------
! Copyright (C) 2019 - Karsten Bolding (BB)                            !
!-----------------------------------------------------------------------
