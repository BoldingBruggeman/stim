!-----------------------------------------------------------------------
!BOP
!
! !MODULE: stim_winton --- Winton thermodynamic ice model
! \label{sec:stim_winton}
!
! !INTERFACE:
   module stim_winton
!
! !DESCRIPTION:
!  The model consists of a zero heat capacity snow layer overlying two equally
!  thick sea ice layers. The upper ice layer has a variable heat capacity to
!  represent brine pockets. The lower ice layer has a fixed heat capacity.
!  The prognostic variables are hs (snow layer thickness), hi (ice layer
!  thickness), T1 and T2, the upper and lower ice layer temperatures located
!  at the midpoints of the layers. The ice model performs two functions, the
!  first is to calculate the ice temperature and the second is to calculate
!  changes in the thickness of ice and snow.
!
!------------------------------------------------------------------------------!
!                                                                              !
!                       THREE-LAYER VERTICAL THERMODYNAMICS                    !
!                                                                              !
! Reference:  M. Winton, 2000: "A reformulated three-layer sea ice model",     !
!            Journal of Atmospheric and Oceanic Technology, 17, 525-531.       !
!                                                                              !
!                                                                              !
!        -> +---------+ <- Ts - diagnostic surface temperature ( <= 0C )       !
!       /   |         |                                                        !
!     hs    |  snow   | <- 0-heat capacity snow layer                          !
!       \   |         |                                                        !
!        => +---------+                                                        !
!       /   |         |                                                        !
!      /    |         | <- T1 - upper 1/2 ice temperature; this layer has      !
!     /     |         |         a variable (T/S dependent) heat capacity       !
!   hi      |...ice...|                                                        !
!     \     |         |                                                        !
!      \    |         | <- T2 - lower 1/2 ice temp. (fixed heat capacity)      !
!       \   |         |                                                        !
!        -> +---------+ <- Tf - base of ice fixed at seawater freezing temp.   !
!                                                                              !
!                                                     Mike Winton (mw@gfdl.gov)!
!------------------------------------------------------------------------------!
!  Note: in this implementation the equations are multiplied by hi to improve
!  thin ice accuracy
!
!  The code is based on the open source sea ice model included in the Modular
!  Ocean Model.
!
! !USES:
   use stim_variables
   use ice_thm_mod
   IMPLICIT NONE
!
!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public                              :: init_stim_winton
   public                              :: do_stim_winton
!KB   public                              :: ice_optics
!
! !PRIVATE DATA MEMBERS:
   real(rk), pointer :: Ts,T1,T2
   real(rk), pointer :: hi,hs,dh1,dh2
!
! !REVISION HISTORY:
!  Original author: Michael Winton
!  Author(s): Adolf Stips, Jesper Larsen and Karsten Bolding
!
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Calculate ice thermodynamics \label{sec:do_ice_winton}
!
! !INTERFACE:
!KB   subroutine init_stim_winton(ice_cover,dz,dt,Tw,S,Ta,precip)
   subroutine init_stim_winton(Ta)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
   real(rk), intent(in)    :: Ta
#if 0
! !INPUT PARAMETERS:
   real(rk), intent(in)    :: dz,dt,Ta,S,precip
!
! !INPUT/OUTPUT PARAMETERS:
   integer, intent(inout)  :: ice_cover
   real(rk), intent(inout) :: Tw
#endif
!
! !LOCAL VARIABLES:
#if 0
   real(rk)        :: Tf    ! seawater freezing temperature (deg-C)
   real(rk)        :: fb    ! heat flux from ocean to ice bottom (W/m^2)
   real(rk)        :: I     ! solar absorbed by upper ice (W/m^2)
   real(rk)        :: evap  ! evaporation of ice (m/s)
   real(rk)        :: snow
   real(rk)        :: trn, pen
   real(rk)        :: dtemp
   real(rk)        :: rho_0 = 1025.
   logical         :: has_ice
!
   real(rk), pointer :: Ts,T1,T2,hi,hs
   real(rk)        :: A, B, dts=0.01
   real(rk)        :: qh,qe,qb
!
   real(rk)        ::  K12,K32
   real(rk)        ::  A1,B1,C1
   real(rk)        ::  Ms,Mb
   real(rk)        ::  h1,h2
   real(rk)        ::  dhs
   real(rk)        ::  T1new,T2new
#endif
!
! !LOCAL PARAMETERS:
   real(rk)        :: kelvin=273.16 ! absolute zero
!EOP
!-----------------------------------------------------------------------
!BOC

   Ts => Tice_surface
   T1 => Tice(1)
   T1 = Tf
   T2 => Tice(2)
   T2 = Tf
!   Tf => Tf
   hs => Hsnow
   hi => Hice
   dh1 => dHis
   dh2 => dHib

   return
!EOC
end subroutine init_stim_winton

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Calculate ice thermodynamics \label{sec:do_ice_winton}
!
! !INTERFACE:
   subroutine do_stim_winton(ice_cover,dz,dt,Tw,S,Ta,precip,Qsw,Qfluxes)
!KB   subroutine do_ice_winton(dt,h,S,sst,T,hs,hi,Ti,surface_melt,bottom_melt)
!
! !DESCRIPTION:
!  This subroutine updates the sea ice prognostic variables. The updated
!  variables are upper ice layer temperature (T1), lower ice layer temperature
!  (T2), snow thickness (hs), and ice thickness (hi).
!
!  The ice model performs this in two steps. First the temperatures are updated
!  and secondly the changes in ice and snow thickness are calculated.
!
!  Any surplus energy that is not used for melting is returned in tmelt and
!  bmelt.
!
!  Evaporation and bottom ablation formation are not included in
!  this version of the model. Furthermore we do not keep an explicit water
!  and salt budget for the sea ice and how that affects the water and salt
!  budgets in the ocean.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   real(rk), intent(in)    :: dz,dt,Ta,S,precip,Qsw
!
! !INPUT/OUTPUT PARAMETERS:
   integer, intent(inout)  :: ice_cover
   real(rk), intent(inout) :: Tw
!
   interface
      subroutine Qfluxes(T,qh,qe,qb)
         integer, parameter                   :: rk = kind(1.d0)
         real(rk), intent(in)                 :: T
         real(rk), intent(out)                :: qh,qe,qb
      end subroutine
   end interface
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
! !OUTPUT PARAMETERS:
!
! !LOCAL VARIABLES:
   real(rk)        :: fb    ! heat flux from ocean to ice bottom (W/m^2)
   real(rk)        :: I     ! solar absorbed by upper ice (W/m^2)
   real(rk)        :: evap  ! evaporation of ice (m/s)
   real(rk)        :: snow
   real(rk)        :: trn, pen
   real(rk)        :: dtemp
   real(rk)        :: rho_0 = 1025.
   logical         :: has_ice
!
   real(rk)        :: A, B, dts=0.01
   real(rk)        :: qh,qe,qb
!
   real(rk)        ::  K12,K32
   real(rk)        ::  A1,B1,C1
   real(rk)        ::  Ms,Mb
   real(rk)        ::  h1,h2
   real(rk), pointer ::  dh1,dh2
   real(rk)        ::  dhs
   real(rk)        ::  T1new,T2new
!
! !LOCAL PARAMETERS:
   real(rk)        :: kelvin=273.16 ! absolute zero
!EOP
!-----------------------------------------------------------------------
!BOC
   !LEVEL0 'do_stim_winton'

   fb = 0._rk; I = 100._rk

!  Calculate seawater freezing temperature
   Tf = -0.0575_rk*S
   h1 = hi/2._rk
   h2 = h1

#if 0
   call ice_optics(albedo_ice, pen, trn, hs, hi, ts, Tf)
#endif

! check this out
   call Qfluxes(Ts,Qh,qh,qb)
   A = qh+qe-qb ! (7-)
   call Qfluxes(Ts+dts,qh,qh,qb)
   B = qh+qe-qb
   B = (B-A)/dts ! (8)
   A = A+I - Ts*B ! (-7)

write(*,*) 'ALB ',albedo_ice,pen,trn

!KB   call ice3lay_temp()
write(*,*) 'AAA ',hi,t1,t2
!https://github.com/mom-ocean/MOM5/blob/08266af73b04d2334be4a52d0c45c174f447cee4/src/ice_sis/ice_model.F90
!                                      hf hfd 
   call ice3lay_temp(0._rk,hi,t1,t2,ts,A,B,I,Tf,fb,dt,Ms,Mb)
#if 1
write(*,*) 'AB ',A,B
write(*,*) 'Hi ',hi
write(*,*) 'T  ',t1,t2,ts
write(*,*) 'M  ',Ms,Mb
#endif
!KBstop 'egon'

#if 0
   call ice3lay_resize(hs, hi, t1, t2, snow, frazil, evap, tmelt, bmelt, &
                       tfw, heat_to_ocn, h2o_to_ocn, h2o_from_ocn,       &
                       snow_to_ice, bablt)
!   call e_to_melt(hs, h1, t1, h2, t2)
#endif

   return

end subroutine do_stim_winton

!-----------------------------------------------------------------------

   end module stim_winton

!-----------------------------------------------------------------------
! Copyright by the GETM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
