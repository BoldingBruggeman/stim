#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: Main ice module
!
! !INTERFACE:
   module ice
!
! !DESCRIPTION:
!  This module provides all variables necessary for the ice
!  calculation and also makes the proper initialisations.
!
! !USES:
   use stim_models
   IMPLICIT NONE
!
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public init_ice, post_init_ice, do_ice, clean_ice
   integer, public :: ice_cover=0 ! 0=no ice, 1=frazil ice, 2=solid ice
!
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
!
! !PUBLIC DATA MEMBERS:

!  the 'ice' namelist
   integer, public                    :: ice_model
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
! !IROUTINE: Initialisation of the ice variables
!
! !INTERFACE:
   subroutine init_stim_yaml()
!
! !DESCRIPTION:
!
! !USES:
   use settings
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  See log for the ice module
!
! !LOCAL VARIABLES:
   class (type_gotm_settings), pointer :: branch
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL1 'init_stim_yaml'
   ice_model    = 0
   branch => settings_store%get_typed_child('surface/ice')
   call branch%get(ice_model, 'model', 'model', default=0, &
                   options=&
                   (/option(0, 'none', 'none'), &
                     option(1, 'Lebedev (1938)', 'Lebedev'), &
                     option(2, 'MyLake', 'MyLake'), &
                     option(3, 'Winton', 'Winton'), &
                     option(4,'Flato', 'Flato')/))

   call branch%get(Hice, 'H', 'initial ice thickness', 'm',default=0._rk)
   call branch%get(ocean_ice_flux, 'ocean_ice_flux', &
                   'ocean->ice heat flux','W/m^2',default=0._rk, display=display_hidden)
   !flato
   call branch%get(nilay, 'nilay', 'number of ice layers', '', default=0)
   !call branch%get(sfall_method, 'sfall_method', 'define how snow fall is determined ','', default=0, &
                 ! options=(/option(1, 'constant snow fall', 'constant'), option(2, 'calculate snowfall from precipitation', &
                 ! 'precipitation')/))
   call branch%get(sfall_method, 'sfall_method', 'define how snow fall is determined ','', default=0)
   call branch%get(const_sfall,'const_sfall ', 'constant snow fall rate', 'm d^-1', default=0._rk)
   call branch%get(dfact,'dfact', 'drift factor', '', default=0._rk)
   call branch%get(depmix ,'depmix', 'prescribed mixed layer depth', '', default=0._rk)
  ! call branch%get(sice_method, 'sice_method','define how sea-ice salinity is to be calculated','', default=0, &
                 ! options=(/option(1, 'cconstant ice salinity', 'constant'), option(2, 'Simple ice salinity profile', 'simple')/))
   call branch%get(sice_method, 'sice_method','define how sea-ice salinity is to be calculated','', default=0)
   call branch%get(const_Sice, 'const_Sice', 'prescribed sea ice salinity', 'ppt', default=0._rk)
   call branch%get(snow_dist,'snow_dist ','logical switch between uniform and Weibull-distributed snow', default=.false.)
   !call branch%get(distr_type,'distr_type', 'integer to chose the type of distribution', '', default=-1, &
                !  options=(/option(-1, 'no specific distribution chosen', 'no'), option(0, 'Weibull distribution', 'Weibull'), &
                !  option(1, 'Reighleigh distribution', 'Reighleigh')/))
   call branch%get(distr_type,'distr_type', 'integer to chose the type of distribution', '', default=-1)
   call branch%get(meltpond,'meltpond', 'If true meltponds are included If false only bare ice is included', default=.false.)
   call branch%get(Ameltmax,'Ameltmax ', 'Maximum meltpond area fraction allowed', '', default=0._rk)
   call branch%get(drainrate,'drainrate', 'Melt pond drainage rate in ','m/d',  default=0._rk)
   call branch%get(hh0, 'hh0', 'initial thickness for S calculation','', default=0._rk)
   call branch%get(ice_hi_i ,'ice_hi_i ', 'initial ice thickness', '', default=0._rk)
   call branch%get(ice_hs_i ,'ice_hs_i ',' initial snow thickness', '',  default=0._rk)
   !call branch%get(albice_method,'albice_method ', 'define how ice albedo is determined ', '', default=0, &
           ! options=(/option(1, 'constant albedo (albice_f and albsnow_f)', 'constant'), &
           ! option(2, 'albice dependent on ice thickness', 'ice thickness constant'), &
           ! option(3, 'albice and albsnow based on eq12&13 of Flato&Brown1996', 'Flato_Brown1996'), &
           ! option(4, 'albice and albsnow dependent on temperature (same as transs/transi formulation)', &
           ! 'transs/transi formulation')/))
   call branch%get(albice_method,'albice_method ', 'define how ice albedo is determined ', '', default=0)
   call branch%get(albice_f, 'albice_f', 'freezing ice albedo','', default=0._rk)
   call branch%get(albmelt,'albmelt',  'melt pond albedo','', default=0._rk)
   call branch%get(albsnow_f,'albsnow_f', 'freezing snow albedo','', default=0._rk)
   call branch%get(albice_m ,'albice_m ', 'melting ice albedo','', default=0._rk)
   call branch%get(albsnow_m, 'albsnow_m',  'melting snow albedo','', default=0._rk)
   call branch%get(transsf,'transsf',  'freezing snow transmission coefficient','', default=0._rk)
   call branch%get(transsm ,'transsm ',  'melting snow transmission coefficient','', default=0._rk)
   call branch%get(transif, 'transif ',  'freezing ice transmission coefficient','', default=0._rk)
   call branch%get(transim, 'transim ', 'melting ice transmission coefficient','', default=0._rk)
   call branch%get(transm,'transm', 'melt pond transmision coefficient','', default=0._rk)
   call branch%get(swkappasm,'swkappasm',  'melting snow extinction coefficient','', default=0._rk)
   call branch%get(swkappasf,'swkappasf ', 'freezing snow extinction coefficient','',  default=0._rk)
   call branch%get(swkappaim,'swkappaim ', 'melting ice extinction coefficient','',  default=0._rk)
   call branch%get(swkappaif, 'swkappaif', 'freezing ice extinction coefficient','', default=0._rk)

   LEVEL2 'done.'
allocate(Tice(2))
   return
   end subroutine init_stim_yaml
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the air--sea interaction module \label{sec:init-air-sea}
!
! !INTERFACE:
   subroutine post_init_stim(Ta,S)
!
! !DESCRIPTION:
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: Ta,S
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL1 'post_init_stim'
   print *, 'ice_model: ', ice_model

   if(Hice .gt. _ZERO_ .and. ice_model /= 0) then
      ice_cover=2
   end if

   Tf = -0.0575_rk*S

   select case (ice_model)
      case(0)
         LEVEL1 'no ice'
#ifdef STIM_LEBEDEV
      case(1)
         call init_stim_lebedev(ice_cover)
#endif
#ifdef STIM_MYLAKE
      case(2)
         call init_stim_mylake()
#endif
#ifdef STIM_WINTON
      case(3)
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
#ifdef STIM_FLATO
      case(4)
#if 1
!KB         allocate(Tice(2))
         call init_stim_flato(Ta)
#else
         LEVEL0 "Flato model is compiled - but execution is disabled"
         LEVEL0 "change line 138 in gotm_stim.F90 - then recompile - "
         LEVEL0 "then do some work to make the Flato ice model work ...."
         LEVEL0 ".... in STIM"
         stop 'post_init_stim(): init_stim_flato()'
#endif
#endif
      case default
         stop 'invalid ice model'
   end select

   call init_stim_variables(ice_model)

   LEVEL2 'done.'
   return
   end subroutine post_init_stim
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: do the ice calculations
!
! !INTERFACE:
   subroutine do_stim(dz,dt,Tw,S,Ta,precip,Qsw,Qfluxes)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)    :: dz,dt,Ta,S,precip,Qsw
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE, intent(inout) :: Tw
!
   interface
      subroutine Qfluxes(T,qh,qe,qb)
         REALTYPE, intent(in)                 :: T
         REALTYPE, intent(out)                :: qh,qe,qb 
      end subroutine
   end interface
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  See log for the ice module
!
! !LOCAL VARIABLES:
   REALTYPE                  :: Tf
!EOP
!-----------------------------------------------------------------------
!BOC

   select case (ice_model)
      case(0)
         Tf = -0.0575*S
         if (Tw .lt. Tf) then
            Tw = Tf
            Hice = 0.1
         else
            Hice = _ZERO_
         end if
#ifdef STIM_LEBEDEV
      case(1)
         call do_stim_lebedev(ice_cover,dt,Tw,S,Ta,precip)
#endif
#ifdef STIM_MYLAKE
      case(2)
         call do_stim_mylake(ice_cover,dz,dt,Tw,S,Ta,precip,Qsw,Qfluxes)
#endif
#ifdef STIM_WINTON
      case(3)
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
#ifdef STIM_FLATO
      case(4)
         if (S .lt. 0.01) then
            LEVEL0 'The FLato ice model is developed for oceanic conditions.'
            LEVEL0 'Very low salinity is not supported - and the principle'
            LEVEL0 'advantage of the model (brine contribution to latent'
            LEVEL0 'heat calculation) is not met.'
            LEVEL0 'Please select another ice model.'
            stop 'do_stim()'
         else
            call do_stim_flato(ice_cover,dz,dt,Tw,S,Ta,precip,Qsw,Qfluxes)
         end if
#endif
      case default
         stop 'invalid ice model'
   end select

   return
   end subroutine do_stim
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Cleaning up the mean flow variables
!
! !INTERFACE:
   subroutine clean_stim()
!
! !DESCRIPTION:
!  De-allocates all memory allocated via init\_ice()
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  See log for the ice module
!
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL1 'clean_ice'

   LEVEL2 'de-allocation ice memory ...'
   LEVEL2 'done.'

   return
   end subroutine clean_stim
!EOC

!-----------------------------------------------------------------------

   end module ice

!-----------------------------------------------------------------------
! Copyright by the STIM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
