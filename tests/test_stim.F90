!
   program test_stim
!
! !DESCRIPTION:
!  program to test the ice models
!
! !USES:
   use stim_variables
   use stim_models
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
#ifdef STIM_LEBEDEV
   write(*,*) 'Lebedev: on'
#else
   write(*,*) 'Lebedev: off'
#endif
#ifdef STIM_MYLAKE
   write(*,*) 'MyLake:  on'
#else
   write(*,*) 'MyLake:  off'
#endif
#ifdef STIM_WINTON
   write(*,*) 'Winton:  on'
#else
   write(*,*) 'Winton:  off'
#endif
#ifdef STIM_FLATO
   write(*,*) 'Flato:  on'
#else
   write(*,*) 'Flato:  off'
#endif

   end program test_stim
!EOC

!-----------------------------------------------------------------------
! Copyright by the STIM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
