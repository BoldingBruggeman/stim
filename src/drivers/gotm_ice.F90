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
!KB   use ice_variables
   use ice_models
   IMPLICIT NONE
!
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public init_ice, post_init_ice, do_ice, clean_ice
   integer, public :: ice_cover=0 ! 0=no ice, 1=frazil ice, 2=solid ice
!
   interface init_ice
      module procedure init_ice_nml
!KB      module procedure init_ice_yaml
   end interface
!
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
   subroutine init_ice_nml(namlst,fn)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                      :: namlst
   character(len=*), intent(in)             :: fn
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  See log for the ice module
!
! !LOCAL VARIABLES:
   integer                   :: rc
   namelist /ice/  ice_model,Hice,sensible_ice_water
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL1 'init_ice'
   ice_model    = 0

!  Read namelist from file.
   open(namlst,file=fn,status='old',action='read',err=80)
   LEVEL2 'reading ice namelists..'
   read(namlst,nml=ice,err=81)
   close (namlst)
   LEVEL2 'done.'

   return
80 FATAL 'I could not open: ',trim(fn)
   stop 'init_ice'
81 FATAL 'I could not read "ice" namelist'
   stop 'init_ice'

   end subroutine init_ice_nml
!EOC

#if 0
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialisation of the ice variables
!
! !INTERFACE:
   subroutine init_ice_yaml(S)
!
! !DESCRIPTION:
!
! !USES:
   use settings
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                     :: S
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  See log for the ice module
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL1 'init_ice_yaml'
   ice_model    = 0
   branch => settings_store%get_typed_child('ice')
   call branch%get(ice_model, 'ice_model', '' ) !&
                !options=(/type_option(1, 'Kondo (1975)'), type_option(2, 'Fairall et al. (1996)')/), default=1)
   call branch%get(Hice, 'Hice', 'total ice thickness', 'm',default=0._rk)
   call branch%get(sensible_ice_water, 'sensible_ice_water','sensible heat flux ice/water','W',default=0._rk)
   LEVEL2 'done.'
   return
   end subroutine init_ice_yaml
!EOC
#endif

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the air--sea interaction module \label{sec:init-air-sea}
!
! !INTERFACE:
   subroutine post_init_ice(Ta,S)
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
   LEVEL1 'post_init_ice'

   if(Hice .gt. _ZERO_) then
      ice_cover=2
   end if

   Tf = -0.0575_rk*S

   select case (ice_model)
      case(0)
         LEVEL1 'no ice'
#ifdef STIM_LEBEDEV
      case(1)
         call init_ice_lebedev(ice_cover)
#endif
#ifdef STIM_MYLAKE
      case(2)
         call init_ice_mylake()
#endif
#ifdef STIM_WINTON
      case(3)
         allocate(Tice(2))
         call init_ice_winton(Ta)
#endif
      case default
         stop 'invalid ice model'
   end select

   call init_ice_variables(ice_model)

   LEVEL2 'done.'
   return
   end subroutine post_init_ice
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: do the ice calculations
!
! !INTERFACE:
   subroutine do_ice(dz,dt,Tw,S,Ta,precip,Qsw,Qfluxes)
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
            ice_cover = 2
         else
            Hice = _ZERO_
            ice_cover = 0
         end if
#ifdef STIM_LEBEDEV
      case(1)
         call do_ice_lebedev(ice_cover,dt,Tw,S,Ta,precip)
#endif
#ifdef STIM_MYLAKE
      case(2)
         call do_ice_mylake(ice_cover,dz,dt,Tw,S,Ta,precip,Qsw,Qfluxes)
#endif
#ifdef STIM_WINTON
      case(3)
         call do_ice_winton(ice_cover,dz,dt,Tw,S,Ta,precip,Qsw,Qfluxes)
#endif
      case default
         stop 'invalid ice model'
   end select

   return
   end subroutine do_ice
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Cleaning up the mean flow variables
!
! !INTERFACE:
   subroutine clean_ice()
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
   end subroutine clean_ice
!EOC

!-----------------------------------------------------------------------

   end module ice

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
