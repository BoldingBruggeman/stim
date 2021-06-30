!-----------------------------------------------------------------------
!BOP
!
! !MODULE: ice model module
!
! !INTERFACE:
   module stim_models
!
! !DESCRIPTION:
!
! !USES:
   use stim_variables
#ifdef STIM_LEBEDEV
   use stim_lebedev, only: init_stim_lebedev, do_stim_lebedev, clean_stim_lebedev
#endif
#ifdef STIM_MYLAKE
   use stim_mylake, only: init_stim_mylake, do_stim_mylake, clean_stim_mylake
#endif
#ifdef STIM_WINTON
   use stim_winton, only: init_stim_winton, do_stim_winton
#endif
#ifdef STIM_FLATO
   use stim_flato, only: init_stim_flato, do_stim_flato, do_ice_uvic ! jp added : do_ice_uvic 
#endif
   IMPLICIT NONE
!
   public

!-----------------------------------------------------------------------

   end module stim_models

!-----------------------------------------------------------------------
! Copyright by the STIM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
