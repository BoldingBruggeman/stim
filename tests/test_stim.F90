!> program to test the ice models
!>
!> author: Karsten Bolding
   program test_stim

   use stim_variables
   use stim_models

   implicit none

   integer  :: i
!-----------------------------------------------------------------------

write(*,*) 'Ice models included in compilation:'
#ifdef STIM_BASAL_MELT
   write(*,*) 'Basal_Melt: on'
#else
   write(*,*) 'Basal_Melt: off'
#endif
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

   end program test_stim

!-----------------------------------------------------------------------
! Copyright by the STIM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
