#include"cppdefs.h"

   MODULE gotm_stim_driver
     !! This module provides the GOTM interface to STIM
     !!
     !! author: Karsten Bolding

   use stim_variables, only: rk
   use stim_models
   IMPLICIT NONE

   private

   public init_ice, post_init_ice, do_ice, clean_ice
   integer, public :: ice_cover=0

   interface init_ice
      module procedure init_stim_yaml
   end interface

   interface post_init_ice
      module procedure post_init_stim
   end interface

   interface do_ice
      module procedure do_stim
   end interface

   interface clean_ice
      module procedure clean_stim
   end interface

   integer, public :: ice_model
     !! select ice model to apply

   ENUM, BIND(C)
      ENUMERATOR :: NONE=0
      ENUMERATOR :: SIMPLE=1
      ENUMERATOR :: BASAL_MELT=2
      ENUMERATOR :: LEBEDEV=3
      ENUMERATOR :: MYLAKE=4
      ENUMERATOR :: WINTON=5
   END ENUM

!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------

   SUBROUTINE init_stim_yaml()
     !! Initialize configuration via entries in a YAML-file
     !! Using the GOTM settings module

   use settings
   IMPLICIT NONE

   class (type_gotm_settings), pointer :: branch
     !! GOTM settings variable
!-----------------------------------------------------------------------
   LEVEL1 'init_stim_yaml'
   ice_model = 0
   branch => settings_store%get_typed_child('surface/ice')
   call branch%get(ice_model, 'model', 'model', default=0, &
                   options=&
                   (/option(0, 'none'), &
                     option(1, 'Simple'), &
                     option(2, 'Basal_Melt'), &
                     option(3, 'Lebedev (1938)'), &
                     option(4, 'MyLake'), &
                     option(5, 'Winton')/))
   call branch%get(Hice, 'H', 'initial ice thickness', 'm',default=0._rk)
   call branch%get(ocean_ice_flux, 'ocean_ice_flux', &
                   'ocean->ice heat flux','W/m^2',default=0._rk, display=display_hidden)
   LEVEL2 'done.'
allocate(Tice(2))
   END SUBROUTINE init_stim_yaml

!-----------------------------------------------------------------------

   SUBROUTINE post_init_stim(Ta,S)

   real(rk), intent(in) :: Ta
      !! Air temperature [C]
   real(rk), intent(in) :: S
      !! Salinity [g/kg]
!-----------------------------------------------------------------------
   LEVEL1 'post_init_stim()'

   if(Hice .gt. 0._rk .and. ice_model /= 0) then
      ice_cover=2
   end if

   Tf = -0.0575_rk*S

   select case (ice_model)
      case(NONE)
         LEVEL1 'no ice'
#if 0
      case(SIMPLE)
#endif
#ifdef STIM_BASAL_MELT
      case(BASAL_MELT)
!KB         call init_stim_basal_melt()
#endif
#ifdef STIM_LEBEDEV
      case(LEBEDEV)
         call init_stim_lebedev(ice_cover)
#endif
#ifdef STIM_MYLAKE
      case(MYLAKE)
         call init_stim_mylake()
#endif
#ifdef STIM_WINTON
      case(WINTON)
#if 1
!KB         allocate(Tice(2))
         call init_stim_winton(Ta)
#else
         LEVEL0 "Winton model is compiled - but execution is disabled"
         LEVEL0 "change line 138 in gotm_stim.F90 - then recompile - "
         LEVEL0 "then do some work to make the Winton ice model work ...."
         LEVEL0 ".... in STIM"
         stop 'post_init_stim(): init_stim_winton()'
#endif
#endif
      case default
         stop 'invalid ice model'
   end select

   call init_stim_variables(ice_model)

   LEVEL2 'done.'
   END SUBROUTINE post_init_stim

!-----------------------------------------------------------------------

   SUBROUTINE do_stim(dz,dt,ustar,Tw,S,Ta,precip,Qsw,Qfluxes)

   !! Arguments
   real(rk), intent(inout)    :: dz
      !! layer thickness [m]
   real(rk), intent(inout)    :: dt
      !! time step [s]
   real(rk), intent(inout)    :: ustar
      !! surface friction velocity [m/s]
   real(rk), intent(inout) :: Tw
      !! water temperature [C]
   real(rk), intent(inout)    :: Ta
      !! air temperature [C]
   real(rk), intent(inout)    :: S
      !! salinity [g/kg]
   real(rk), intent(inout)    :: precip
      !! precipitation [mm?]
   real(rk), intent(inout)    :: Qsw
      !! short wave radiation [W/m^2]
!
   interface
      SUBROUTINE Qfluxes(T,qh,qe,qb)
         use, intrinsic :: iso_fortran_env
         integer, parameter :: rk=real64

         real(rk), intent(in) :: T
           !! temperature [C]
         real(rk), intent(out) :: qh
           !! latent heat [W/m^2]
         real(rk), intent(out) :: qe
           !! sensible heat [W/m^2]
         real(rk), intent(out) :: qb 
           !! net longwave radiation [W/m^2]
      END SUBROUTINE
   END interface

   real(rk) :: Tf
!-----------------------------------------------------------------------
   select case (ice_model)
      case(SIMPLE)
         Tf = -0.0575_rk*S
         if (Tw .lt. Tf) then
            Tw = Tf
            Hice = 0.1_rk
         else
            Hice = 0._rk
         end if
#ifdef STIM_BASAL_MELT
      case(BASAL_MELT)
         call do_stim_basal_melt(dz,ustar,Tw,S)
#endif
#ifdef STIM_LEBEDEV
      case(LEBEDEV)
         call do_stim_lebedev(ice_cover,dt,Tw,S,Ta,precip)
#endif
#ifdef STIM_MYLAKE
      case(MYLAKE)
         call do_stim_mylake(ice_cover,dz,dt,Tw,S,Ta,precip,Qsw,Qfluxes)
#endif
#ifdef STIM_WINTON
      case(WINTON)
         if (S .lt. 0.01) then
            LEVEL0 'The Winton ice model is developed for oceanic conditions.'
            LEVEL0 'Very low salinity is not supported - and the principle'
            LEVEL0 'advantage of the model (brine contribution to latent'
            LEVEL0 'heat calculation) is not met.'
            LEVEL0 'Please select another ice model.'
            stop 'do_stim()'
         else
            call do_stim_winton(ice_cover,dz,dt,Tw,S,Ta,precip,Qsw,Qfluxes)
         end if
#endif
      case default
         stop 'invalid ice model'
   end select
   END SUBROUTINE do_stim

!-----------------------------------------------------------------------

   SUBROUTINE clean_stim()
     !!  De-allocates all memory allocated via init\_ice()
!-----------------------------------------------------------------------
   LEVEL1 'clean_ice'

   LEVEL2 'de-allocation ice memory ...'
   LEVEL2 'done.'
   END SUBROUTINE clean_stim
!EOC

!-----------------------------------------------------------------------

   END MODULE gotm_stim_driver

!-----------------------------------------------------------------------
! Copyright by the STIM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
