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
   use stim_variables
   use stim_models
   use domain, only: imin,imax,jmin,jmax,az
   IMPLICIT NONE
!
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public init_ice, post_init_ice, do_ice, clean_ice
!
   interface init_ice
      module procedure init_stim_nml
!KB      module procedure init_stim_yaml
   end interface
!
!
! !PUBLIC DATA MEMBERS:
!
!  the 'ice' namelist
   integer, public                    :: stim_model
   ! 0=no ice, 1=frazil ice, 2=solid ice
   integer, public, allocatable       :: ice_cover(:,:)
   REALTYPE, public, allocatable      :: Hice_2d(:,:)
   REALTYPE, public, allocatable      :: Tf_2d(:,:)
   ! Lebedev
   REALTYPE, public, allocatable      :: fdd_2d(:,:)
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
   subroutine init_stim_nml(namlst,fn)
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

   end subroutine init_stim_nml
!EOC

#if 0
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialisation of the ice variables
!
! !INTERFACE:
   subroutine init_stim_yaml(S)
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
   LEVEL1 'init_stim_yaml'
   ice_model    = 0
   branch => settings_store%get_typed_child('ice')
   call branch%get(ice_model, 'ice_model', '' ) !&
                !options=(/option(1, 'Kondo (1975)'), option(2, 'Fairall et al. (1996)')/), default=1)
   call branch%get(Hice, 'Hice', 'total ice thickness', 'm',default=0._rk)
   call branch%get(sensible_ice_water, 'sensible_ice_water','sensible heat flux ice/water','W',default=0._rk)
   LEVEL2 'done.'
   return
   end subroutine init_stim_yaml
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
   REALTYPE, intent(in)                :: Ta(:,:),S(:,:)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
! !LOCAL VARIABLES:
   integer                             :: i,j,rc
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL1 'post_init_ice'

   if(Hice .gt. _ZERO_) then
      ice_cover=2
   end if

   allocate(Tf_2d(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'post_init_ice:: Error allocating memory (Tf_2d)'
   Tf_2d = -0.0575_rk*S

   allocate(Hice_2d(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'post_init_ice:: Error allocating memory (Hice_2d)'
   Hice_2d = _ZERO_

   select case (ice_model)
      case(0)
         LEVEL1 'no ice'
      case(1)
return
         allocate(fdd_2d(E2DFIELD),stat=rc)
         if (rc /= 0) stop 'post_init_ice:: Error allocating memory (fdd_2d)'
         fdd_2d = _ZERO_
         do j=jmin,jmax
            do i=imin,imax
              call init_stim_lebedev(ice_cover(i,j))
              fdd_2d(i,j) = fdd
              Hice_2d(i,j) = Hice
            end do
         end do
#if 0
      case(2)
         call init_stim_mylake()
      case(3)
         allocate(Tice(2))
         call init_stim_winton(Ta)
#endif
      case default
         stop 'invalid ice model'
   end select

   call init_stim_variables(ice_model)

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
#define _TESTING_
#ifdef _TESTING_
   subroutine do_ice(dz,dt,Tw,S,Ta)
#else
   subroutine do_ice(dz,dt,Tw,S,Ta,precip,Qsw,Qfluxes,julianday,secondsofday,longitude, &
                     latitude,I_0,airt,airp,hum,u10,v10,cloud,rho,rho_0,longwave_radiation, &
                     hum_method,fluxes_method,albedo,heat)
#endif
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
#ifdef _TESTING_
   REALTYPE, intent(in)    :: dz(:,:),dt,S(:,:),Ta(:,:)
#else
   REALTYPE, intent(in)    :: dz(:,:),dt,Ta(:,:),S(:,:),precip(:,:),Qsw(:,:)
   REALTYPE, intent(in)    :: longitude,latitude,I_0,airt,airp,hum,u10,v10,cloud,rho,rho_0,albedo,heat
   integer, intent(in)     :: julianday,secondsofday
   integer, intent(in)     :: longwave_radiation,hum_method,fluxes_method
  
#endif
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE, intent(inout) :: Tw(:,:)
!
#ifndef _TESTING_
   interface
      subroutine Qfluxes(T,qh,qe,qb)
         REALTYPE, intent(in)                 :: T
         REALTYPE, intent(out)                :: qh,qe,qb 
      end subroutine
   end interface
#endif
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  See log for the ice module
!
! !LOCAL VARIABLES:
   integer                             :: i,j,rc
   REALTYPE                            :: airt
!EOP
!-----------------------------------------------------------------------
!BOC
   select case (ice_model)
      case(0)
!KB         Tf_2d = -0.0575_rk*S
         do j=jmin,jmax
            do i=imin,imax
STDERR i,j
               if (Tw(i,j) .lt. Tf_2d(i,j)) then
                  Tw(i,j) = Tf_2d(i,j)
                  Hice_2d(i,j) = 0.1
                  ice_cover(i,j) = 2
               else
                  Hice_2d(i,j) = _ZERO_
                  ice_cover(i,j) = 0
               end if
            end do
         end do
      case(1)
STDERR Tw(imax/2,jmax/2)
return
         do j=jmin,jmax
            do i=imin,imax
               fdd = fdd_2d(i,j)
               if (i .lt. imax/2) then
                  airt = -2.
               end if
#if 0
               call do_stim_lebedev(ice_cover(i,j),dt,Tw(i,j),S(i,j),Ta(i,j),_ZERO_)
#else
               call do_stim_lebedev(ice_cover(i,j),dt,Tw(i,j),S(i,j),airt,_ZERO_)
#endif
              fdd_2d(i,j) = fdd
              Hice_2d(i,j) = Hice
            end do
         end do
#if 0
      case(2)
         call do_stim_mylake(ice_cover,dz,dt,Tw,S,Ta,precip,Qsw,Qfluxes)
      case(3)
         call do_stim_winton(ice_cover,dz,dt,Tw,S,Ta,precip,Qsw,Qfluxes)
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
! Copyright by the STIM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
