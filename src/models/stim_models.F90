   module stim_models
     !! Collection of the various ice models in STIM
     !!
     !! Author: Karsten Bolding

   use stim_variables

#if 0
#ifdef STIM_OBSICE
   use stim_obsice, only: init_stim_obsice, do_stim_obsice
#endif
#endif
#ifdef STIM_BASAL_MELT
   use stim_basal_melt, only: init_stim_basal_melt, do_stim_basal_melt, clean_stim_basal_melt
#endif
#ifdef STIM_LEBEDEV
   use stim_lebedev, only: init_stim_lebedev, do_stim_lebedev, clean_stim_lebedev
#endif
#ifdef STIM_MYLAKE
   use stim_mylake, only: init_stim_mylake, do_stim_mylake, clean_stim_mylake
#endif
#ifdef STIM_WINTON
   use stim_winton, only: init_stim_winton, do_stim_winton
#endif
   IMPLICIT NONE

!-----------------------------------------------------------------------

   end module stim_models

!-----------------------------------------------------------------------
! Copyright by the STIM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
