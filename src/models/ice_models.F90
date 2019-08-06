!-----------------------------------------------------------------------
!BOP
!
! !MODULE: ice model module
!
! !INTERFACE:
   module ice_models
!
! !DESCRIPTION:
!
! !USES:
   use ice_variables
#ifdef STIM_LEBEDEV
   use ice_lebedev, only: init_ice_lebedev, do_ice_lebedev, clean_ice_lebedev
#endif
#ifdef STIM_MYLAKE
   use ice_mylake, only: init_ice_mylake, do_ice_mylake, clean_ice_mylake
#endif
#ifdef STIM_WINTON
   use ice_winton, only: init_ice_winton, do_ice_winton
#endif
   IMPLICIT NONE
!
   public

!-----------------------------------------------------------------------

   end module ice_models

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
