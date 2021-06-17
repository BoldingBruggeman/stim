!-----------------------------------------------------------------------
!BOP
!
! !MODULE: stim_flato --- flato thermodynamic ice model
! \label{sec:stim_flato}
!
! !INTERFACE:
   module stim_flato
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
! Reference:  M. Winton , 2000: "A reformulated three-layer sea ice model",     !
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
   public                              :: init_stim_flato
   public                              :: do_stim_flato
!

! !PRIVATE DATA MEMBERS:
   real(rk), pointer :: Ts,T1,T2
   real(rk), pointer :: hi,hs,dh1,dh2
   real(rk), pointer :: trn
   real(rk)          :: pen
   real(rk), pointer :: tmelt,bmelt
   real(rk), pointer :: fb    ! heat flux from ocean to ice bottom (W/m^2)
!

! LOCAL VARIABLES: 
!   rhosnow     - snow density (kg m-3)
   real(rk), public :: rhosnow

!fixed size variables
!   Iceflux     - a 2-element array of time-step averaged boundary fluxes
!                  returned from heat conduction scheme, defined as positive
!                 downward in accordance with vertical coordinate sign
!                 convention (W m-2)
   real(rk), dimension(2), public :: Iceflux
!    bctype 	- 2 element array defining upper and lower boundary         
!                  condition type. bctype(1) refers to upper boundary,
!                  bctype(2) refers to lower boundary. Negative value
!                  indicates temperature boundary condition, positive
!                 value indicates flux boundary condition.
   real(rk), dimension(2) :: bctype
!    bcs 	    - 2 element array containing boundary condition values.  
!                  bcs(1) refers to upper boundary, bcs(2) refers to lower
!                  boundary. If temperature boundary condition, units: (K),
!                  if flux boundary condition, units: (W m-2)
   real(rk), dimension(2) :: bcs
!    dti        - timestep in the ice model (usually set to ocean time step)
   real(rk) :: dti
!
!*****fluxes
!   qb           - long wave back radiation (in-out) (W m-2)	
   real(rk) :: qb
!   qh           - latent heat flux into ice (W m-2)		
   real(rk) :: qh
!   qe           - sensible heat flux into ice (W m-2)		
   real(rk) :: qe
!   tx,ty        - surface stress components in x and y direction (Pa)!check
   real(rk) :: tx,ty
!   PenSW        - short wave radiation that penetrates the surface (W m-2)
   real(rk), public :: PenSW
!   fluxt        - net flux at surface of ice/snow slab calculated from
!                  surface energy budget and used as upper boundary condition
!                  in heat conduction solution. NOTE: this does not include
!                  'PenSW' which is taken as a distributed source within the
!                  slab (W m-2)
   real(rk) :: fluxt
!   simass       - ice mass per unit area (kg m-2)
   real(rk), public :: simass
!   snmass       - snow mass per unit area (kg m-2)
   real(rk), public :: snmass
!   simasso      - ice mass per unit area at previous timestep(kg m-2)
   real(rk) :: simasso
!   snmasso      - ice mass per unit area at previous timestep (kg m-2) !snow mass? 
   real(rk) :: snmasso
!   Ts           - upper surface temperature (K) ! jpnote 
   !real(rk) :: Ts
!   Tsav         - average snow layer temperature (K)
   real(rk) :: Tsav  

! Ice salinity
!   ice_salt     - logical variable to turn on/off the salt profile scheme
    logical,  public :: ice_salt=.false.
! Atmospheric forcing
!   sfall        - snow fall rate (m s-1)
   real(rk), public :: sfall
!   airtk        - surface air temperature (K)
   real(rk) :: airtk
! coefficients for 
   real(rk), dimension(5,3) :: C
! coefficients for RHS vector
   real(rk), dimension(5,3) :: R
!  dto           - time step from GOTM (dt)
   real(rk) :: dto
! nslay number of snow layers
   integer, public :: nslay
!Snow_dist: Snowdistribution variables:
!  Asnow         - Area which is covered with snow 
   real(rk), public :: Asnow
!  Aice          - Area which is covered with ice
   real(rk), public :: Aice
!  Amelt         - Area which is covered with melt pond
   real(rk), public :: Amelt
!  hsmax         - Maximal height of snow for the calculations of Weibull-distributed snow
   real(rk) :: hsmax
!  albice        - Albedo of ice
   real(rk) :: albice
!  albsnow       - Albedo of snow
   real(rk):: albsnow
!   meltmass      - melt pond mass per unit area (kg m-2)
   real(rk), public  :: meltmass
!   meltmasso     - melt pond mass per unit area at previous timestep(kg m-2)
   real(rk) :: meltmasso
!  Pi !NSnote read from gotm?
   real(rk), parameter :: pi=4.D+00*atan(1.D+00)


! !REVISION HISTORY:
!  Original author: Michael Winton 
!  Author(s): Adolf Stips, Jesper Larsen and Karsten Bolding
!
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Calculate ice thermodynamics \label{sec:do_ice_flato}
!
! !INTERFACE:
!KB   subroutine init_stim_flato(ice_cover,dz,dt,Tw,S,Ta,precip)
   subroutine init_stim_flato(Ta)
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
!
! !LOCAL PARAMETERS:
!EOP
!-----------------------------------------------------------------------
!BOC
   trn => transmissivity
   Ts => Tice_surface
   Ts = Ta
   T1 => Tice(1)
   T1 = Ta
   T2 => Tice(2)
   T2 = Tf
   hs => Hsnow
   hi => Hice
   dh1 => dHis
   dh2 => dHib
   tmelt => surface_ice_energy
   bmelt => bottom_ice_energy
   fb => ocean_ice_flux


   print *, hlaymin
   print *, rhoice
   print *, Tfreezi
   print *, rCpmix
   print *, Hfi
   print *, hsmin
   print *, theta
   print *, sigma
   print *, epsilon
   print *, PenFrac
   print *, hlaymin
   print *, rhoscold
   print *, rhoswarm
   print *, rhowaterfresh
   print *, rhoice
   print *, kelvin
   print *, Tmelts
   print *, Tmelti
   print *, Condfi
   print *, rhoCpfi 
   print *, rCpmix
   print *, Hfi
   print *, Hfw
   print *, swkappa
   print *, Tfreezi
   print *, nlmax
   print *, 'print nilay from yaml', nilay
   print *, 'print sfall_method from yaml', sfall_method
   print *, 'print const_sfall from yaml', const_sfall
   print *, 'print dfact from yaml', dfact
   print *, 'print depmix from yaml', depmix
   print *, 'print sice_method from yaml', sice_method
   print *, 'print dist_type from yaml', distr_type
   print *, 'print snow_dist from yaml', snow_dist
   print *, 'print const_Sice from yaml', const_Sice
   print *, 'print meltpond from yaml', meltpond
   print *, 'print Ameltmax from yaml', Ameltmax
   print *, 'print drainrate from yaml', drainrate
   print *, 'print hh0 from yaml', hh0
   print *, 'print ice_hi_i from yaml', ice_hi_i
   print *, 'print ice_hs_i from yaml', ice_hs_i
   print *, 'print albice_method from yaml', albice_method
   print *, 'print albice_f from yaml', albice_f
   print *, 'print albmelt from yaml', albmelt
   print *, 'print albsnow_f from yaml', albsnow_f
   print *, 'print albice_m from yaml', albice_m
   print *, 'print albsnow_m from yaml', albsnow_m
   print *, 'print transsf from yaml', transsf
   print *, 'print transsm  from yaml', transsm 
   print *, 'print swkappasm from yaml', swkappasm
   print *, 'print swkappasf from yaml', swkappasf
   print *, 'print swkappaim from yaml',  swkappaim
   print *, 'print  swkappaif  from yaml',  swkappaif 



   return
!EOC
end subroutine init_stim_flato

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Calculate ice thermodynamics \label{sec:do_ice_flato}
!
! !INTERFACE:
   subroutine do_stim_flato(ice_cover,dz,dt,Tw,S,Ta,precip,Qsw,Qfluxes)
!KB   subroutine do_ice_flato(dt,h,S,sst,T,hs,hi,Ti,surface_melt,bottom_melt)
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
   real(rk)        :: I     ! solar absorbed by upper ice (W/m^2)
   real(rk)        :: evap  ! evaporation of ice (m/s)
   real(rk)        :: snow
!
   real(rk)        :: A, B, dts=0.01
   real(rk)        :: qe,qh,qb
!
   real(rk)        :: h1,h2
   real(rk)        :: ts_new
   real(rk)        :: frazil
   real(rk)        :: heat_to_ocn, h2o_to_ocn, h2o_from_ocn, snow_to_ice
!
! !LOCAL PARAMETERS:
!EOP
!-----------------------------------------------------------------------
!BOC
   !LEVEL0 'do_stim_flato'
   tmelt = 0._rk
   bmelt = 0._rk

   ! Calculate seawater freezing temperature
   Tf = -0.0575_rk*S

   if (ice_cover .gt. 0) then
      call ice_optics(albedo_ice, pen, trn, hs, hi, ts, Tf)
      I = Qsw*(1._rk-trn)
      h1 = hi/2._rk
      h2 = h1

      ! check this out
      call Qfluxes(Ts,qe,qh,qb)
      A = -(qe+qh+qb) ! (7-)
      call Qfluxes(Ts+dts,qe,qh,qb)
      B = -(qe+qh+qb)
      B = (B-A)/dts ! (8)
      A = A-I - Ts*B ! (-7)

      !https://github.com/mom-ocean/MOM5/blob/08266af73b04d2334be4a52d0c45c174f447cee4/src/ice_sis/ice_model.F90
      call ice3lay_temp(hs,hi,t1,t2,ts_new,A,B,pen*I,Tf,fb,dt,tmelt,bmelt)
      ts = ts_new
 !     frazil = 0._rk
      Hfrazil = 0._rk
   else
      frazil = -(Tw-Tf)*dz*Cw
      if (frazil .gt. 0._rk) Hfrazil = frazil/(rho_ice*L_ice)
   end if

   if (ice_cover .gt. 0 .or. frazil .gt. 0._rk) then
      call ice3lay_resize(hs, hi, t1, t2, snow, frazil, evap, tmelt, bmelt, &
                          tf, heat_to_ocn, h2o_to_ocn, h2o_from_ocn,       &
                          snow_to_ice)
   !                       snow_to_ice, bablt)
!write(*,*) 'CC ',heat_to_ocn, h2o_to_ocn, h2o_from_ocn
   end if

hs = 0._rk
   if (hi .gt. 0._rk) then
      ice_cover = 2
   else
      ice_cover = 0
   end if

   return

end subroutine do_stim_flato

!-----------------------------------------------------------------------

   end module stim_flato

!-----------------------------------------------------------------------
! Copyright by the GETM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
