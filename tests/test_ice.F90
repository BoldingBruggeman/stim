!
   program test_ice
!
! !DESCRIPTION:
!  program to test the ice models
!
!  To build:
!  make test_ice
!  To execute:
!  ./test_ice
!  To plot:
!  python $GOTM_BASE/scr/plot_airsea.py
!
! !USES:
   use ice_variables
   use ice_models
   implicit none
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding

! !LOCAL VARIABLES
   integer  :: i
!EOP
!-----------------------------------------------------------------------
!BOC

write(*,*) 'Ice models included:'
#ifdef ICE_LEBEDEV
   write(*,*) 'Lebedev: on'
#else
   write(*,*) 'Lebedev: off'
#endif
#ifdef ICE_MYLAKE
   write(*,*) 'MyLake:  on'
#else
   write(*,*) 'MyLake:  off'
#endif
#ifdef ICE_WINTON
   write(*,*) 'Winton:  on'
#else
   write(*,*) 'Winton:  off'
#endif

   end program test_ice
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
