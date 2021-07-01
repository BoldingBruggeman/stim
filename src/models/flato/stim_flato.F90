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
   use ice_thm_mod_flato
   IMPLICIT NONE
!
!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
!WINTON
   public                              :: init_stim_flato
   public                              :: do_stim_flato  !winton specific jp 

!FLATO 
   public                              :: do_ice_uvic !jp added 
!
! !PRIVATE DATA MEMBERS:
!WINTON 
   real(rk), pointer :: Ts,T1,T2
   !real(rk), pointer :: hi,hs,dh1,dh2 !jpnote: commented orig
   real(rk), pointer :: dh1,dh2
   real(rk), pointer :: trn
   real(rk)          :: pen
   real(rk), pointer :: tmelt,bmelt
   real(rk), pointer :: fb    ! heat flux from ocean to ice bottom (W/m^2)
!Flato 
   real(rk), pointer :: hi,hs  ! jpnote: keep for flato

   !for pointer 
   !real(rk), pointer :: ice_hi, ice_hs 
   !real(rk), pointer :: swr_0,precip_i,sfall_i 
   !real(rk), pointer :: TopMelt,Botmelt,TerMelt,TopGrowth,BotGrowth 
   !real(rk), pointer :: Fh,Ff,Fs 
   !real(rk), pointer :: Sice_bulk,Hmix,Aice_i,Asnow_i,Amelt_i,ice_hm 
   !real(rk), pointer :: Cond,rhoCp,Sint,dzi,zi,Told,Pari 
   !real(rk), pointer :: uvic_Tice
   !parb,parui
   !Tice => ice_uvic_Tice
  

!FLATO
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
   real(rk) :: hsmax = hsmin
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

   real(rk) :: dummy = 0 ! test variables for testing functions and subroutines -jp

!-----------------end of flato vars-------------------------------------------------- 


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
   real(rk), intent(in)    :: Ta    !winton specific 
#if 0
! !INPUT PARAMETERS:
   real(rk), intent(in)    :: dz,dt,Ta,S,precip
!
! !INPUT/OUTPUT PARAMETERS:
   integer, intent(inout)  :: ice_cover
   real(rk), intent(inout) :: Tw
#endif


!Flato 
! !LOCAL VARIABLES:
!
   !integer             :: k,rc     !from init_ice_uvic - jp 
! !LOCAL PARAMETERS:
!EOP
!-----------------------------------------------------------------------
!BOC
   if (runwintonflato .eq. 1) then  !winton specific 


   hs => Hsnow !keep for flato
   hi => Hice  !keep for flato    !should I change to ice_hi ? 
   trn => transmissivity
   Ts => Tice_surface
   Ts = Ta
   T1 => Tice(1)
   T1 = Ta
   T2 => Tice(2)
   T2 = Tf
   dh1 => dHis
   dh2 => dHib
   tmelt => surface_ice_energy
   bmelt => bottom_ice_energy
   fb => ocean_ice_flux

   

   else

   !hs => Hsnow !keep for flato  ???
   !hi => Hice  !keep for flato 
  
   ! for flato ??? jpnote 
   !--------------
   !ts => Tice(1)                     !ice_uvic_ts=ice_uvic_Tice(1)
   !tb => Tice(nilay)                 !ice_uvic_tb=ice_uvic_Tice(nilay)
   !ice_uvic_parb => Pari(nilay)      !ice_uvic_parb=ice_uvic_Pari(nilay)
   !ice_uvic_parui => Pari(nilay+1)   !ice_uvic_parui=ice_uvic_Pari(nilay+1)
   !-------------
!
!  FLATO 
!
!-------------------------------------------------------------------------------------
!      copy paste from init_ice_uvic FROM ice_uvic.F90
!-------------------------------------------------------------------------------------
 
! !DESCRIPTION: !copy paste from ice_uvic.f90 without the namelist initialization 
!
! !USES:
  ! IMPLICIT NONE
!
! !INPUT PARAMETERS:

! !REVISION HISTORY:
!  Original author(s): Karsten Bolding, Nadja Steiner based on original model from Greg Flato
!
!
! !LOCAL VARIABLES:
  ! integer             :: k,rc  !jpnote: added aboce

   ! check if the snow distribution is set
   if(snow_dist .and. distr_type .eq. -1) then
     print*,'No specific snow distribution was chosen'
     stop
   end if

!  The snow fall rate 
   select case (sfall_method)
      case (1)
        ! LEVEL2 'Using snow fall rate= ', const_sfall
         sfall=const_sfall/86400.D+00 !convert to [m/s]
      case (2) 
         sfall = 0.00D+00 ! initialize snowfall from precipitation to 0.0
      case default
   end select
 
      snmass=ice_hs_i*rhoscold
      simass=ice_hi_i*rhoice
      meltmass=0.0
      drainrate=drainrate*rhowaterfresh/86400.D+00 
                                   !conversion from [m/d] to [kg/m^2/s]

! The sea ice salinity calculation
   select case (sice_method)
      case (1)
!         LEVEL2 'Using constant bulk salinity',const_Sice
      case (2) 
!         LEVEL2 'Using simple ice salinity calculation (Vancoppenolle 2009)'
         ice_salt=.true.
      case default
   end select
! check if the snow distribution is set snow_dist
   if(snow_dist .and. distr_type .eq. -1) then
     !LEVEL2 'No specific snow distribution was chosen'
     stop
   end if
!return

!-------------------------------------------------------------------------
!
! --> commented out parts: not relevant to new gotm 
!
   !  initialize namelist variables to reasonable defaults.
  ! ice_method=0

   !  The different ice models
   !select case (ice_method)
   !case (0)
     ! LEVEL2 'No ice calculations included'
   !case (1)
      !LEVEL2 'Clip heat-fluxes if SST < freezing point (function of S)'
    !  ice_layer=0
   !case (2)
      !LEVEL2 'Thermodynamic ice model adapted from Winton'
    !  hsnow=0;hice=0;ice_T1=0;ice_T2=0
     ! ice_ts=0;ice_tmelt=0;ice_bmelt=0
  ! case (3)
      !LEVEL2 'Thermodynamic ice model adapted from Flato&Brown, 1992, UVic'
     ! call init_ice_uvic(namlst)
#if 0
   allocate(ice_uvic_Tice(nilay+1),stat=rc)
   if (rc /= 0) STOP 'init_ice: Error allocating (ice_uvic_Tice)'
   ice_uvic_Tice=0
      do k=1,nilay+1
         ice_uvic_Tice(k)=245.+(Tfreezi-245.)*float(k-1)/float(nilay)
      enddo
   allocate(ice_uvic_Cond(nilay),stat=rc)
   if (rc /= 0) STOP 'init_ice: Error allocating (ice_uvic_Cond)'
   ice_uvic_Cond =0
   allocate(ice_uvic_rhoCp(nilay),stat=rc)
   if (rc /= 0) STOP 'init_ice: Error allocating (ice_uvic_rhoCp)'
   ice_uvic_rhoCp =0
   allocate(ice_uvic_Sint(nilay+1),stat=rc)
   if (rc /= 0) STOP 'init_ice: Error allocating (ice_uvic_Sint)'
   ice_uvic_Sint =0
   allocate(ice_uvic_dzi(nilay),stat=rc)
   if (rc /= 0) STOP 'init_ice: Error allocating (ice_uvic_dzi)'
   ice_uvic_dzi =0
   allocate(ice_uvic_zi(nilay+1),stat=rc)
   if (rc /= 0) STOP 'init_ice: Error allocating (ice_uvic_zi)'
   ice_uvic_zi =0
   allocate(ice_uvic_Told(nilay+1),stat=rc)
   if (rc /= 0) STOP 'init_ice: Error allocating (ice_uvic_Told)'
   ice_uvic_Told =0
   allocate(ice_uvic_Pari(nilay+1),stat=rc)
   if (rc /= 0) STOP 'init_ice: Error allocating (ice_uvic_Pari)'
   ice_uvic_Pari =0
   allocate(ice_uvic_dum(nilay+1),stat=rc)
   if (rc /= 0) STOP 'init_ice: Error allocating (ice_uvic_dum)'
   allocate(ice_uvic_dzice(nilay),stat=rc)
   if (rc /= 0) STOP 'init_ice: Error allocating (ice_uvic_dzi)'
      do k=1,nilay
         ice_uvic_dzice(k)=float(k)
      enddo
   allocate(ice_uvic_zice(nilay+1),stat=rc)
   if (rc /= 0) STOP 'init_ice: Error allocating (ice_uvic_zice)'
      do k=1,nilay+1
         ice_uvic_zice(k)=float(k)
      enddo

   ice_uvic_dum =0
   hsnow=0;hice=0;ice_uvic_hm=0
   ice_uvic_ts=273.16D+00;ice_uvic_tb=273.16D+00;ice_uvic_Fh=0
   ice_uvic_swr_0=0;ice_uvic_precip_i=0;ice_uvic_sfall_i=0
   ice_uvic_parb=0;ice_uvic_parui=0;
   ice_uvic_Ff=0;ice_uvic_Fs=0
   ice_uvic_Sicebulk=6.0D+00
   ice_uvic_topmelt=0;ice_uvic_botmelt=0
   ice_uvic_termelt=0;ice_uvic_topgrowth=0
   ice_uvic_botgrowth=0;ice_uvic_Hmix=0 
   ice_uvic_Aice=0;ice_uvic_Asnow=0 
   ice_uvic_Amelt=0 

#endif 

      !not needed as long as ice_uvic is passed into subroutines called from do_ice_uvic 
      !pointers from flato  --> so I dnot  
      
      !ice_hi => Hsnow
      !ice_hs => Hice
      !swr_0 => ice_uvic_swr_0
      !precip_i => ice_uvic_precip_i
      !sfall_i => ice_uvic_sfall_i
      !parb,parui
      !TopMelt => ice_uvic_topmelt
      !Botmelt => ice_uvic_botmelt 
      !TerMelt => ice_uvic_termelt 
      !TopGrowth => ice_uvic_topgrowth
      !BotGrowth => ice_uvic_botgrowth
      !Fh => ice_uvic_Fh
      !Ff => ice_uvic_Ff
      !Fs => ice_uvic_Fs
      !Sice_bulk => ice_uvic_Sicebulk
      !Hmix => ice_uvic_Hmix
      !Aice_i => ice_uvic_Aice
      !Asnow_i => ice_uvic_Asnow
      !Amelt_i => ice_uvic_Amelt
      !ice_hm =>  ice_uvic_hm
      
      !uvic_Tice => ice_uvic_Tice(nilay+1) !test
      !Tice => ice_uvic_Tice(nilay+1)
      !Cond => ice_uvic_Cond(nilay)
      !rhoCp => ice_uvic_rhoCp(nilay)
      !Sint => ice_uvic_Sint(nilay+1) 
      !dzi => ice_uvic_dzi(nilay) 
     ! zi => ice_uvic_zi(nlmax)
      !Told => ice_uvic_Told(nilay+1)
      !Pari => ice_uvic_Pari(nilay+1)

   !!case default
!end select


!-----------------------------------------------------------------------
   endif
!-------------------------------------------------------------------------------

#if 0 
   !printing Vars for Testing purposes
   print *, '----------------------------------------------'
   print *, 'public data members from ice_uvic'
   print *, '----------------------------------------------'
   print *, 'hlaymin', hlaymin
   print *, 'rhoice',rhoice
   print *, 'Tfreezi',Tfreezi
   print *, 'rCpmix',rCpmix
   print *, 'Hfi',Hfi
   print *, 'hsmin',hsmin
   print *, 'theta',theta
   print *, 'sigma',sigma
   print *, 'epsilon',epsilon
   print *, 'PenFrac',PenFrac
   print *, 'hlaymin',hlaymin
   print *, 'rhoscold',rhoscold
   print *, 'rhoswarm',rhoswarm
   print *, 'rhowaterfresh',rhowaterfresh
   print *, 'rhoice',rhoice
   print *, 'kelvin',kelvin
   print *, 'Tmelts',Tmelts
   print *, 'Tmelti',Tmelti
   print *, 'Condfi',Condfi
   print *, 'rhoCpfi',rhoCpfi 
   print *, 'rCpmix',rCpmix
   print *, 'Hfi',Hfi
   print *, 'Hfw',Hfw
   print *, 'swkappa',swkappa
   print *, 'freezi',Tfreezi
   print *, 'nlmax',nlmax
   print *, '----------------------------------------------'
   print *, 'local Variables'
   print *, '----------------------------------------------'
   print *, 'rhosnow',rhosnow 
   print *, 'Iceflux', Iceflux
   print *, 'bctype', bctype
   print *, 'bcs', bcs
   print *, 'dti', dti
   print *, 'qb', qb
   print *, 'qh', qh
   print *, 'qe', qe
   print *, 'tx', tx
   print *, 'ty', ty
   print *, 'PenSW', PenSW
   print *, 'fluxt', fluxt
   print *, 'simass', simass
   print *, 'snmass',  snmass
   print *, 'simasso', simasso
   print *, 'snmasso', snmasso
   print *, 'Ts',  Ts
   print *, 'Tsav', Tsav
   print *, 'ice_salt', ice_salt
   print *, 'sfall', sfall
   print *, 'dfact', dfact
   print *, 'depmix', depmix
   print *, 'airtk', airtk
   print *, 'C', C
   print *, 'R', R
   print *, 'dto', dto
   print *, 'nslay', nslay
   print *, 'Asnow', Asnow
   print *, 'Aice', Aice
   print *, 'hsmax',  hsmax
   print *, 'meltmass', meltmass
   print *, 'meltmasso', meltmasso
   print *, 'pi', pi
   print *, '----------------------------------------------'
   print *, 'Yaml Variables'
   print *, '----------------------------------------------'
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
   print *, '----------------------------------------------'
   print *, 'ice.F90 vars'
   print *, '----------------------------------------------'
   print *, ''
   print *, 'ice_uvic_dum', ice_uvic_dum
   print *,'hsnow', hsnow
   print *,'hice', hice
   print *,'ice_uvic_hm', ice_uvic_hm
   print *,'ice_uvic_ts', ice_uvic_ts
   print *,'ice_uvic_tb', ice_uvic_tb
   print *,'ice_uvic_Fh', ice_uvic_Fh
   print *, 'ice_uvic_Ff ',ice_uvic_Ff 
   print *, 'ice_uvic_Fs',ice_uvic_Fs
   print *,'ice_uvic_swr_0', ice_uvic_swr_0
   print *,'ice_uvic_precip_i', ice_uvic_precip_i
   print *,'ice_uvic_sfall_i', ice_uvic_sfall_i
   print *, 'ice_uvic_parb', ice_uvic_parb
   print *,'ice_uvic_parui', ice_uvic_parui
   print *,'ice_uvic_Ff',ice_uvic_Ff
   print *,'ice_uvic_Fs',ice_uvic_Fs
   print *,'ice_uvic_Sicebulk',ice_uvic_Sicebulk
   print *,'ice_uvic_topmelt',ice_uvic_topmelt
   print *,'ice_uvic_botmelt',ice_uvic_botmelt
   print *,'ice_uvic_termelt',ice_uvic_termelt
   print *,'ice_uvic_topgrowth',ice_uvic_topgrowth
   print *,'ice_uvic_botgrowth',ice_uvic_botgrowth
   print *,'ice_uvic_Hmix',ice_uvic_Hmix
   print *,' ice_uvic_Aice', ice_uvic_Aice
   print *, 'ice_uvic_Asnow', ice_uvic_Asnow
   print *, 'ice_uvic_Amelt',ice_uvic_Amelt
   print *, 'ice_uvic_Tice ',ice_uvic_Tice 
   print *, 'Tice', Tice
   print *, 'Tice(nilay)', Tice(nilay)
   print *, 'ice_uvic_Cond',ice_uvic_Cond
   print *, 'ice_uvic_rhoCp ',ice_uvic_rhoCp 
   print *, 'ice_uvic_Sint',ice_uvic_Sint
   print *, 'ice_uvic_dzi ',ice_uvic_dzi 
   print *, 'ice_uvic_zi ',ice_uvic_zi 
   print *, 'ice_uvic_Told ',ice_uvic_Told 
   print *, 'ice_uvic_Pari',ice_uvic_Pari
   print *, 'ice_uvic_dzice',ice_uvic_dzice
   print *, 'ice_uvic_zice',ice_uvic_zice
   print *, ''

#endif 

   return
!EOC
end subroutine init_stim_flato

!-----------------------------------------------------------------------

!BOP
!-----------------------------------------------------------------------
!   !WINTON SPECIFIC SUBROUTINE 
!-----------------------------------------------------------------------
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




!----------------------------------------------------------------------
!----------------------------------------------------------------------
!  !BEGINNING OF FLATO SPECIFIC 
!----------------------------------------------------------------------
!----------------------------------------------------------------------


!BOP 
!
! !ROUTINE: Calculate ice thermodynamics \label{sec:do_ice_uvic}
!
! !INTERFACE: 
!subroutine do_ice_uvic(dt,dz,julianday,secondsofday,longitude,latitude, &
                      ! I_0,airt,airp,hum,u10,v10,,precip,cloud,sst,sss,rhowater,rho_0,back_radiation_method, &
                     !  hum_method,fluxes_method,alb,heat)
subroutine do_ice_uvic(dto,h,julianday,secondsofday,lon,lat, &
                        I_0,airt,airp,rh,u10,v10,precip,cloud, &
                        TSS,SSS,rhowater,rho_0, &
                        back_radiation_method,hum_method,fluxes_method, &
                        ice_hi,ice_hs,ice_hm,Tice,Cond,rhoCp,Sint,dzi,zi, &
                        Pari,Told,alb,heat,Fh,Ff,Fs,Sice_bulk,TopMelt,BotMelt,&
                        TerMelt,TopGrowth,BotGrowth,Hmix,Aice_i,Asnow_i,Amelt_i,swr_0,precip_i,sfall_i)
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   real(rk), intent(in)     :: dto ! ocean timestep (sec) ! dto ??? 
   !real(rk), intent(in)     :: dz !--> h from meanflow ???
   real(rk), intent(in)    :: h  ! sea surface layer thickness --> is this the same as dz ???
   
   integer, intent(in)     :: julianday ! this julian day --> from time
   integer, intent(in)     :: secondsofday ! seconds for this day -->  from time
   real(rk), intent(in)    :: lon    ! longitude for this point --> longitude from gotm
   real(rk), intent(in)    :: lat ! latitude for this point --> latitude from gotm 
   real(rk), intent(inout)  :: I_0   ! shortwave radiation at sea surface  
   real(rk), intent(in)     :: airt  ! 2m temperature
   real(rk), intent(in)     :: airp  ! sea surface pressure
   !real(rk), intent(in)     :: hum   ! relative humidity from airsea
   real(rk), intent(in)    :: rh    ! relative humidity --> is this the same as hum ??? 
   real(rk), intent(inout)  :: u10   ! 10 m wind u-component
   real(rk), intent(inout)  :: v10   ! 10 m wind v-component
   real(rk), intent(inout)  :: precip! freshwater precipitation (m/s)
   real(rk), intent(in)     :: cloud ! cloud cover
   real(rk), intent(inout)  :: TSS     ! sea surface temperature
   real(rk), intent(in)     :: SSS     ! sea surface salinity
   real(rk), intent(in)      :: rhowater   ! sea surface layer density --> called rho(nlev) in new code
   real(rk), intent(in)      :: rho_0 ! reference density --> from meanflow
   integer, intent(in)       :: back_radiation_method ! method for LW   !read in from namelist in airsea --> defined as a local variable in airsea
   integer, intent(in)       :: hum_method ! method for humidity
   integer, intent(in)       :: fluxes_method ! method for fluxes
    
   ! !INPUT/OUTPUT PARAMETERS:
   real(rk), intent(inout)   :: ice_hi    ! ice thickness (m)
   real(rk), intent(inout)   :: ice_hs    ! snow thickness (m)
   real(rk), intent(inout)   :: ice_hm    ! meltpond thickness (m)
   real(rk), intent(inout)   :: Tice(nilay+1)  ! ice layer temperature Tice(nilay +1)(deg-C)
   real(rk), intent(inout)   :: Cond(nilay)  ! thermal conductivities defined at the 
                               ! centre of each layer Cond(nilay)(W m-1 K-1)
   real(rk), intent(inout)   :: rhoCp(nilay) ! volumetric heat capacities defined at 
                               ! the centre of each layer rhoCp(nilay)(J m-3 K-1)
   real(rk), intent(inout)   :: Sint(nilay+1) ! internal heat source due to penetrating 
                                ! short wave radiation Sint(nilay)(W m-3)
   real(rk), intent(inout)   :: dzi(nilay) !layer thicknesses dzi(nilay)(m)
   real(rk), intent(inout)   :: zi(nlmax) !layer interface depths zi(nilay+1)(m)
   real(rk), intent(inout)   :: Told(nilay+1) !ice temperature two time steps 
                                ! previous to calculation of outgoing 
                                ! long-wave flux in SEBUDGET Told (nilay+1) (K)
   real(rk), intent(inout)   :: alb   ! surface albedo - water, ice/snow-total
   real(rk), intent(inout)   :: heat  ! surface heat flux

!NSnote, check - maybe adjust units...
   real(rk), intent(inout)  ::  Sice_bulk   ! bulk ice salinity (ppt)
   real(rk), intent(inout)  ::  Fh   ! interface heat flux (W/m2)
   real(rk), intent(out)  ::  Ff   ! interface freshwater flux (m s-1)
   real(rk), intent(out)  ::  Fs   ! interface salt flux - (ppt m s-1)
   real(rk), intent(out)   :: Pari(nilay+1) 
                       !photosynthatically available radiation in ice (W m-2)

! !OUTPUT PARAMETERS:
   real(rk), intent(out)     :: TopMelt ! top melting - ice mass melted at the surface (snow+ice)  at time step(m)
   real(rk), intent(out)     :: Botmelt ! bottom melting - ice mass melted at the ice bottom at time step(m) 
   real(rk), intent(out)     :: TerMelt ! internal melting - ice mass melted in the ice interior at time step (m) 
   real(rk), intent(out)     :: TopGrowth ! top growth ice mass growth at slab surface due to snow submersion (m)
   real(rk), intent(out)     :: BotGrowth ! bottom growth - ice mass growth at the ice bottom at time step (m)     
   real(rk), intent(out)     :: Hmix !  transferred energy - check  (m)
!   Hmix        - mixed layer heat storage (J m-2)	=======> accounts only for 
! the SWR which crosses the ice slab and reach the water. keep it for now
   real(rk), intent(out)     :: Aice_i,Asnow_i,Amelt_i ! ice area fraction which is : open ice, snow and 
                                                ! meltpond, respectively
   real(rk), intent(out)     :: swr_0,precip_i,sfall_i !H! incidental SWR,snowfall



   ! LOCAL VARIABLES 
   !real(rk)        :: h1,h2  !point to hice and hsnow whose values are taken from stim_variables.F90 

   real(rk)          :: dmsi ! dmsi - new ice formation at open water [kg m-2]
   integer           :: yy,mm,dd,j,k
   real(rk)          :: p_tmp ! p_tmp ! for rain to snow conversion
!KB
   real(rk)        :: evap  ! evaporation of ice (m/s)
   real(rk)        :: ohflux !heat flux from ocean into ice underside (W m-2)
!NSnote, not sure what evap is for


   !for fh .. as pointers 
   !declare as local 
   !real(rk) :: ice_hi, ice_hs 
   !real(rk) :: Fh,Ff,Fs 
   !real(rk) :: swr_0,precip_i,sfall_i 
   !real(rk) :: TopMelt,Botmelt,TerMelt,TopGrowth,BotGrowth 
   !real(rk) :: Sice_bulk,Hmix,Aice_i,Asnow_i,Amelt_i,ice_hm 
   !real(rk), dimension(:), allocatable  :: Cond,rhoCp,Sint,dzi,zi,Told,Pari,uvic_Tice


#if 0
   print *,'Hice',Hice
   print *,'Hsnow',Hsnow
   print *,'ice__uvic_hm',ice_uvic_hm
   print *,'ice_uvic_Tice(nilay+1)',ice_uvic_Tice(nilay+1)
   print *,'ice_uvic_Cond(nilay)',ice_uvic_Cond(nilay)
   print *,'ice_uvic_rhoCp(nilay)',ice_uvic_rhoCp(nilay)
   print *,'ice_uvic_Sint(nilay+1) ',ice_uvic_Sint(nilay+1) 
   print *,'ice_uvic_dzi(nilay)',ice_uvic_dzi(nilay)
   print *,'ice_uvic_zi(nlmax)',ice_uvic_zi(nlmax)
   print *, 'ice_uvic_Pari(nilay+1)',ice_uvic_Pari(nilay+1)
   print *,'ice_uvic_Told(nilay+1)',ice_uvic_Told(nilay+1)
   print *,'ice_uvic_Fh',ice_uvic_Fh
   print *,'ice_uvic_Ff',ice_uvic_Ff
   print *,'ice_uvic_Fs',ice_uvic_Fs
   print *,'ice_uvic_Sicebulk',ice_uvic_Sicebulk
   print *,'ice_uvic_TopMelt',ice_uvic_TopMelt
   print *,'ice_uvic_BotMelt', ice_uvic_BotMelt
   print *, 'ice_uvic_TerMelt',ice_uvic_TerMelt
   print *,'ice_uvic_TopGrowth',ice_uvic_TopGrowth
   print *,'ice_uvic_BotGrowth',ice_uvic_BotGrowth
   print *,'ice_uvic_Hmix',ice_uvic_Hmix
   print *,'ice_uvic_Aice',ice_uvic_Aice
   print *,'ice_uvic_Asnow',ice_uvic_Asnow
   print *,'ice_uvic_Amelt',ice_uvic_Amelt
   print *,'ice_uvic_swr_0',ice_uvic_swr_0
   print *,'ice_uvic_precip_i',ice_uvic_precip_i
   print *,'ice_uvic_sfall_i',ice_uvic_sfall_i
#endif 


   !call open_water(nilay,I_0,Sice_bulk,Hmix,Tice,depmix,sst,Fh,heat,precip,precip_i)
   call open_water(nilay,I_0,Sice_bulk,Hmix,Tice,depmix,TSS,Fh,heat,precip,precip_i)
                     

   !call nr_iterate(hum_method,back_radiation_method,fluxes_method,nilay,&
                  !airt,hum,cloud,I_0,Told,Tice,Pari,&
                 ! Sice_bulk,ice_hi,ice_hs,dzi,Cond,rhoCp,zi,Sint,&
                  !latitude,u10,v10,precip,airp,evap,alb)
   call nr_iterate(hum_method,back_radiation_method,fluxes_method,nilay,&
                  airt,rh,cloud,I_0,ice_uvic_Told,ice_uvic_Tice,ice_uvic_Pari,&
                  Sice_bulk,ice_hi,ice_hs,dzi,Cond,rhoCp,zi,Sint,&
                  lat,u10,v10,precip,airp,evap,alb)

   
   !call cndiffus()

end subroutine do_ice_uvic 
! EOC
!-----------------------------------------------------------------------



subroutine therm1d(nilay,Sice_bulk,ice_hi,ice_hs,dzi,Cond,rhoCp,zi,Sint,Pari,Tice,I_0)   
   !nilay,Sice_bulk,ice_hi,ice_hs,dzi,Cond,rhoCp,zi,Sint,Pari,Tice,I_0)
! !USES:
   
   IMPLICIT NONE 

! !INPUT PARAMETERS:
    
   ! !INPUT PARAMETERS:
   integer, intent(in)     :: nilay
   real(rk), intent(in)    :: I_0   !incoming swr for snow_dist 

! !OUTPUT PARAMETERS:
   real(rk),intent(out)   :: ice_hi,ice_hs
   real(rk),intent(out)   :: dzi(nilay),zi(nilay+1)
   real(rk),intent(out)   :: Sint(nilay+1),Cond(nilay),rhoCp(nilay) 
   real(rk),intent(out)   :: Pari(nilay+1) 

! !INOUT PARAMETERS:
   real(rk),intent(inout)  :: Tice(nilay+1)
   real(rk),intent(inout)  :: Sice_bulk
!
!BOC
! !LOCAL VARIABLES:
!   bctype 	- 2 element array defining upper and lower boundary
!                 condition type. bctype(1) refers to upper boundary,
!                 bctype(2) refers to lower boundary. Negative value
!                 indicates temperature boundary condition, positive
!                 value indicates flux boundary condition.
   real(rk)                 ::      bctype(2)
!   bcs 	- 2 element array containing boundary condition values.
!                 bcs(1) refers to upper boundary, bcs(2) refers to lower
!                 boundary. If temperature boundary condition, units: (K),
!                 if flux boundary condition, units: (W m-2)
!                 NOTE: flux is positive downward due to sign convention.
   real(rk)                  ::      bcs(2)
!   hlay	- ice or ice/snow thickness partitioned by ice layers (m)
   real(rk)                   ::      hlay
!   Ticeav - average ice temperature (interim)
   real(rk)                   ::      Ticeav
!   nslay	- number of snow layers
!   integer                   ::      nslay
!   Sice_bulk		- ice salinity used in parameterization of heat capacity
!    		  and conductivity (ppt)
   real(rk)                   ::      snlaythick,Tbot,thicklay,Ttop
!   fluxtsplit  - averaged surface flux over all 3 surface conditions
   real(rk)                   :: fluxtsplit
   integer                   ::      k
!   hzero   - minimum ice thickness to avoid dividing by zero elsewhere in code
   real(rk), parameter           :: hzero=1.D-03


   call cndiffus(bctype,bcs,nilay,dzi,rhoCp,Cond,Sint,Tice)

end subroutine therm1d

!-----------------------------------------------------------------------

! !IROUTINE: 
!
! !INTERFACE:
subroutine cndiffus(bctype,bcs,nilay,dzi,rhoCp,Cond,Sint,Tice)
   !bctype,bcs,nilay,dzi,rhoCp,Cond,Sint,Tice)
! !USES:
   IMPLICIT NONE

! !INPUT PARAMETERS:
   real(rk), intent(in)        :: bctype(2),bcs(2)
   integer, intent(in)         :: nilay
   real(rk), intent(in)        :: dzi(nilay),rhoCp(nilay),Cond(nilay)
   real(rk), intent(in)        :: Sint(nilay+1)
! !INOUT PARAMETERS:
   real(rk), intent(inout)        :: Tice(nilay+1)
! !LOCAL VARIABLES:
! coefficients for tridiagonal matrix
   real(rk)                  ::      C(nlmax+1,3)
! coefficients for RHS vector
   real(rk)                  ::      R(nlmax+1)    
! ratio of conductivity to heat capacity
   real(rk)                  ::      rCpav(nlmax+1)      
! array to store previous temperature 
   real(rk)                  ::      To(nlmax+1)     
   integer                   ::      j,l


   call trisol(C,R,nilay+1,Tice)

end subroutine cndiffus


!-----------------------------------------------------------------------
subroutine trisol(C,R,n,Tice)

! !USES:
   IMPLICIT NONE

   ! !INPUT PARAMETERS:
   real(rk), intent(in)        :: C(nlmax+1,3)
   real(rk), intent(in)        :: R(nlmax+1)
!       n -  size of system of equations
   integer, intent(in)         :: n   !--> corresponds to nilay+1 
   real(rk), intent(inout)     :: Tice(n)

! !LOCAL VARIABLES:
! coefficients for 
   integer,parameter         :: nmax=100
   real(rk),parameter        :: eps=1.D-12
   real(rk)                  :: D1      
   real(rk)                  :: temp(nmax)
   integer                   :: j  

!integer :: n
   !n = nilay+1  ??? jpnote 
   !print *,'trisol n', n 

end subroutine trisol
!-----------------------------------------------------------------------

subroutine growthtb(rhowater,nilay,rhoCp,dzi,Tice,TopGrowth,TerMelt, &
                     TopMelt,BotGrowth,BotMelt,Hmix,ohflux)

 ! !USES:
   IMPLICIT NONE
!
!INPUT PARAMETERS:
      real(rk), intent(in)   :: rhowater
      integer, intent(in)    :: nilay
      real(rk), intent(in)   :: rhoCp(nilay),dzi(nilay)
      real(rk), intent(in)   :: ohflux
!OUTPUT PARAMETERS:
   
!INOUT PARAMETERS:
      real(rk), intent(inout)   :: BotGrowth,BotMelt
      real(rk), intent(inout)   :: TopGrowth,TerMelt,TopMelt,Hmix
      real(rk), intent(inout)   :: Tice(nilay+1)                    


! !LOCAL VARIABLES:
! coefficients for 
   real(rk)                  ::  dhib,dhst,dhs,dhi,Ts
! coefficients for 
   real(rk)                  ::  Tl    
!   snsub       - submerged snow mass per unit area (kg m-2)
   real(rk)                  ::  snsub    
!   porewat     - mass per unit area or pore water frozen into submerged
!                 snow (kg m-2)
   real(rk)                 ::  porewat   
   integer                   ::  k  


   call surfmelt(Tice(1),TopMelt)

end subroutine growthtb

!-----------------------------------------------------------------------

subroutine sebudget(hum_method,back_radiation_method,fluxes_method,&
                     TTss,airt,rh,cloud,ice_hi,ice_hs,&
                     lat,u10,v10,precip,airp,evap)
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
      integer, intent(in)       :: back_radiation_method ! method for LW
      integer, intent(in)       :: hum_method ! method for humidity
      integer, intent(in)       :: fluxes_method ! method for fluxes
      real(rk), intent(in)      :: TTss,airt,rh,cloud
      real(rk), intent(in)      :: lat   ! latitude for this point
      real(rk), intent(in)      :: u10   ! 10 m wind u-component
      real(rk), intent(in)      :: v10   ! 10 m wind v-component
      real(rk), intent(in)      :: airp  ! sea surface pressure
! !INOUT PARAMETERS:
      real(rk), intent(inout)   :: precip! freshwater precipitation (m/s)
      real(rk), intent(inout)   :: evap! freshwater evaporation (m/s)
      real(rk), intent(inout)   :: ice_hi,ice_hs
! !REVISION HISTORY:

      ! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
!      LEVEL1'sebudget'


!...Calculate ice and snow thickness from mass and density
      ice_hs=snmass/rhosnow
      ice_hi=simass/rhoice
!
!...Calculate surface fluxes from meteorological data using bulk formulae
!     Recalculate surface fluxes for sea ice conditions
! NSnote  TTss is in Kelvin!!!
      call humidity(hum_method,rh,airp,TTss-kelvin,airt)
      call longwave_radiation(back_radiation_method, &
                          lat,TTss,airt+kelvin,cloud,qb)
      call airsea_fluxes(fluxes_method,.false.,.false., &
                         TTss-kelvin,airt,u10,v10,precip,evap,tx,ty,qe,qh)

! evap is set to zero for now....as for ice_winton, not sure what it is...
      evap = 0

!...Calculate energy flux into surface
!
!.... long wave flux
      fluxt=qb
!
!....sensible heat flux
      fluxt=fluxt+qe
!
!....latent heat flux
      fluxt=fluxt+qh
!
!
!    LEVEL1'end sebudget'

   return

end subroutine sebudget

!-----------------------------------------------------------------------

subroutine surfmelt(Ts,TopMelt)

   real(rk), intent(in)  :: Ts  !Tice(1)
   real(rk)              :: hmeltinput !???

! !OUTPUT PARAMETERS:

! !INOUT PARAMETERS:
   real(rk), intent(inout) :: TopMelt

! !LOCAL VARIABLES:
! coefficients for 
   real(rk)                 ::   Fdiv,dmsw,dmi,Fdivice,Fdivsnow  

end subroutine surfmelt

!-----------------------------------------------------------------------



subroutine nr_iterate(hum_method,back_radiation_method,fluxes_method,&
                        nilay,airt,rh,cloud,I_0,Told,Tice,Pari,&
                        Sice_bulk,ice_hi,ice_hs,dzi,Cond,rhoCp,zi,Sint,&
                        lat,u10,v10,precip,airp,evap,alb)

                        !                        
! !DESCRIPTION:
!
!*************************************************************************
!  Subroutine to perform Newton-Raphson iteration to solve for a
!  surface temperature consistent with the surface energy budget
!  including conductive heat flux. If no solution is found in
!  initial Newton-Raphson iteration, a bisection algorithm is
!  used.  NOTE: if uppermost layer is greater than 0.1m thick,
!  surface temperature is defined as the temperature at the surface
!  Tice(1); however, if uppermost layer is less than 0.1m thick, surface
!  temperature is defined as the average temperature of the top layer.
!  This prevents rapid fluctuations in surface temperature and so
!  allows a longer time step to be taken.
!*************************************************************************

   !hum_method,back_radiation_method,fluxes_method,&
   !nilay,airt,rh,cloud,I_0,Told,Tice,Pari,&
   !Sice_bulk,ice_hi,ice_hs,dzi,Cond,rhoCp,zi,Sint,&
   !lat,u10,v10,precip,airp,evap,alb)

   !Told -> ice_uvic_Told     Told(nilay+1)
   !Tice -> ice_uvic_Tice     Tice(nilay+1)
   !Pari -> ice_uvic_Pari     Pari(nilay+1)
   !Sice_bulk -> ice_uvic_Sicebulk
   !ice_hi -> hice
   !ice_hs -> hsnow
   !dzi -> ice_uvic_dzi      dzi(nilay)
   !Cond -> ice_uvic_cond     Cond(nilay)
   !rhoCp -> ice_uvic_rhoCp    rhoCp(nilay) 
   !!zi -> ice_uvic_zi        zi(nilay+1)
   ! Sint -> ice_uvic_sint    Sint(nilay+1)
   !lat -> latitude
   !alb -> labedo 
   !rh --> hum 

! !USES:
   IMPLICIT NONE

! !INPUT PARAMETERS:  
   integer, intent(in)       :: back_radiation_method ! method for LW
   integer, intent(in)       :: hum_method ! method for humidity
   integer, intent(in)       :: fluxes_method ! method for fluxes
   integer, intent(in)         :: nilay
   real(rk), intent(in)      :: airt
   real(rk), intent(in)      :: rh
   real(rk), intent(in)      :: cloud 
   real(rk), intent(in)      :: lat   ! latitude for this point
   real(rk), intent(in)      :: u10   ! 10 m wind u-component
   real(rk), intent(in)      :: v10   ! 10 m wind v-component
   real(rk), intent(in)      :: airp  ! sea surface pressure
! !OUTPUT PARAMETERS:
   real(rk), intent(out)       :: Sice_bulk
   real(rk), intent(out)       :: Told(nilay+1)
   real(rk), intent(out)       :: dzi(nilay),zi(nilay+1)
   real(rk), intent(out)       :: Sint(nilay+1),Cond(nilay),rhoCp(nilay) 
   real(rk), intent(out)       :: alb 
   real(rk), intent(out)       :: Pari(nilay+1)
! !INOUT PARAMETERS:
   real(rk), intent(inout)     :: Tice(nilay+1)
   real(rk), intent(inout)     :: ice_hi,ice_hs
   real(rk), intent(inout)     :: I_0
   real(rk), intent(inout)     :: precip! freshwater precipitation (m/s)
   real(rk), intent(inout)     :: evap! freshwater evaporation (m/s)


! !LOCAL VARIABLES:
! coefficients for 
   real(rk)                 ::    fTs,dTemp,Tsp,fTsdT1,fTsdT2
   real(rk)                  ::    fprime,Tsnew,Error,Ts1,fTs1
   real(rk)                  ::    dTs,Tsu,fTsu,Ts2,Tsm,fTsm,fTsl, Tsl,Ts  
   integer                   ::    nrit,ksearch,kb,l
   real(rk),parameter        ::    toler=1.D-02 

#if 0

!      LEVEL1'nr_iterate'
!
!
!...Save initial temperature profile
!
   do l=1,nilay+1
      Told(l)=Tice(l)
   enddo
   
!...First guess at surface temperature
!
   Ts=Tice(1)

!-----------------------------------------------------------------------------
!--- Start initial Newton-Raphson iteration - do a maximum of 5 iterations ---
!-----------------------------------------------------------------------------
   
   do nrit=1,5
      !      
      !......calculate surface energy budget terms
               call sebudget(hum_method,back_radiation_method,fluxes_method,&
                             Ts,airt,rh,cloud,ice_hi,ice_hs,&
                             lat,u10,v10,precip,airp,evap)
               call  albedo_ice(fluxt,I_0,PenSW,alb,ice_hs,ice_hi,Ts)
      !
      !
      !......restore intial temperature profile as therm1d takes a forward time step
               do l=1,nilay+1
                  Tice(l)=Told(l)
               enddo
      !
      !......solve unsteady heat conduction equation to get new temperature profile 
      !......and bottom flux after one time step
               call therm1d(nilay,Sice_bulk,ice_hi,ice_hs,dzi,Cond,&
                     rhoCp,zi,Sint,Pari,Tice,I_0)
      !
      !......now find mismatch in assumed and computed surface temperature
               fTs=Ts-Tice(1)
      !
      !......stop iterating if mismatch is less than tolerance
               if(abs(fTs).lt.toler) then
      !           print*,'Newton-Raphson scheme met tolerance after', &
      !                          nrit-1,'iterations'
                 go to 987
               endif
      !
      !......re-do above calculations with slightly increased and decreased
      !......initial temperature to calculate first derivative by centred
      !......finite difference
               dTemp=0.1
               Tsp=Ts+dTemp
               call sebudget(hum_method,back_radiation_method,fluxes_method,&
                             Tsp,airt,rh,cloud,ice_hi,ice_hs,&
                             lat,u10,v10,precip,airp,evap)
               call  albedo_ice(fluxt,I_0,PenSW,alb,ice_hs,ice_hi,Tsp)
      !
               do l=1,nilay+1
                  Tice(l)=Told(l)
               enddo
               call therm1d(nilay,Sice_bulk,ice_hi,ice_hs,dzi, &
                  Cond,rhoCp,zi,Sint,Pari,Tice,I_0)
               fTsdT1=Tsp-Tice(1)
               Tsp=Ts-dTemp
               call sebudget(hum_method,back_radiation_method,fluxes_method,&
                             Tsp,airt,rh,cloud,ice_hi,ice_hs,&
                             lat,u10,v10,precip,airp,evap)
               call  albedo_ice(fluxt,I_0,PenSW,alb,ice_hs,ice_hi,Tsp)
      !
               do l=1,nilay+1
                  Tice(l)=Told(l)
               enddo
               call therm1d(nilay,Sice_bulk,ice_hi,ice_hs,dzi, &
                Cond,rhoCp,zi,Sint,Pari,Tice,I_0)
               fTsdT2=Tsp-Tice(1)
      !
      !......if derivative vanishes, take solution obtained at this point
      !......and don't iterate any further
               if(abs(fTsdT1-fTsdT2).lt.1.D-02) then
                 print*,'fTsdT1~=fTsdT2 after ',nrit-1,' iterations'
                 go to 987
               endif
      !
      !......otherwise, compute derivative and hence next guess from
      !......Newton-Raphson formula
               fprime=(fTsdT1-fTsdT2)/(2.D+00*dTemp)
               Tsnew=Ts-fTs/fprime
      !
      !......check error and stop iterating if tolerance is met
               Error=Tsnew-Ts
               if(abs(Tsnew-Ts).lt.toler) then
      !           print*,'Newton-Raphson scheme met tolerance after', &
      !                          nrit,'iterations'
                 go to 987
               endif
      !
      !......if another iteration is required, update surface temperature
      !......and try again
               Ts=Tsnew
            enddo
      !
      !...If iteration has not reached a solution at this point, (which is
      !...rare) we have to resort to a cruder, brute force approach. So, ...
      !
            print*,'***Newton-Raphson scheme failed: Error =',Error
            print*,'***Starting method of Bisection'
  

            987 continue
            
#endif        
   return

end subroutine nr_iterate 

!-----------------------------------------------------------------------

! !INTERFACE:
subroutine open_water(nilay,I_0,Sice_bulk,Hmix,Tice,hSS,TSS,Fh,heat,precip,precip_i)

! !USES:                 
   IMPLICIT NONE
! !INPUT PARAMETERS:
   integer,intent(in)   :: nilay
   real(rk),intent(in) :: hSS !--> [[ depmix ]] just named something different when passed to subroutine 
!INOUT PARAMETERS
   real(rk), intent(inout) :: Sice_bulk
   real(rk), intent(inout) :: Hmix
   real(rk), intent(inout) :: Tice(nilay+1)
   real(rk), intent(inout) :: Fh ! interface heat flux (W/m2)
   real(rk), intent(inout) :: TSS ! ocean surface temperature (C)   ! sea surface temperature
   real(rk), intent(inout) :: I_0
   real(rk), intent(inout) :: heat
   real(rk), intent(in)    :: precip !H!
   real(rk), intent(out)   :: precip_i !H!
! !LOCAL VARIABLES:
! coefficients for 
   real(rk)                  ::  Tmix,dmsi
   real(rk)                  ::  Sni   ! Sni - salinity of new ice
   integer                   ::  k,cont,i

   !print *, 'open_water vars', hSS, I_0, heat, precip, Tmix, dmsi, Sni, k, cont, i


   !LEVEL1'open_water'
   !
   !
   !   Calculate the surface energy budget components
   !   We use the thickness of the first layer instead
   !   of 'depmix' in calculations
   ! check, can we get mixed layer and Tmix from GOTM? 
   ! watch consistency for removal of heat
   ! NSnote why set to Tice(1), ok if there is ice, but what if not? 
   ! Answer:should be Tmix if open water before - 
   ! used so Tmix (or TSS) does not need to be transferred 
   
         Sni=6.0D00
   
   !      Tmix=Tfreezi+Hmix/(rCpmix*depmix)
         Tmix=(TSS+kelvin)+ice_uvic_Hmix/(rCpmix*hSS)
   !  Set all 'ice' temperatures to be the mixed layer temperature
   ! (only the surface value has any real meaning in this case, 
   ! since there is actually no ice yet)
   !  
         do k=1,nilay+1
            Tice(k)=Tmix
         enddo
   
   !   Check if mixed layer temperature is below freezing. If so,
   !   create some ice which is assumed to be isothermal at the 
   !   freezing temperature initially 
   
         if(Tmix .lt. Tfreezi) then
            
            dmsi=(Tfreezi-Tmix)*hSS*rCpmix/Hfi
            simass=simass+dmsi
            Hmix=0.D+00
            Tmix=Tfreezi 
   !reset TSS here!!!  
            TSS=Tmix-kelvin
   
   !   New ice salinity
            
            if(ice_salt) then
               Sni=4.606D+00+0.91603D+00/hh0
            endif
            Sice_bulk=Sni
   
   !   New ice temperature
   
            do k=1,nilay+1
               Tice(k)=Tfreezi
            enddo
   
   ! need to set heat and I_0 to zero since this negtive heat has already 
   ! been used for ice growth, otherwise it will be applied twice - second time in temperature.F90 in gotm  - NScheck I_0 for light!
            heat=0.D00
            I_0=0.D00 
   
   ! Note, what is with Fh??? even if only for accounting
   !         Fh=
          else
            !just open water, no ice existing, no melt or growth
            Tice=Tfreezi
            dmsi=0 
            Fh=0 
            snmass=0
            simass=0
            Iceflux(2)=0
            Sice_bulk=Sni
            I_0=I_0
            heat=heat
            Hmix=0
         endif
         precip_i=precip !H!
      
      return

end subroutine open_water

!-----------------------------------------------------------------------
!call saltice_prof_simple(dti,nilay,SSS,simasso,snmasso,meltmasso, &
!rhoice,rhosnow,rhowater,rhowaterfresh,zi,Tice,Sice_bulk,Ff,Fs, &
 !TopGrowth,BotGrowth,TopMelt,BotMelt,TerMelt)
subroutine saltice_prof_simple(dti,nilay,SSS,simasso,snmasso,meltmasso, &
                                 rhoice,rhosnow,rhow,rhof,zi,Tice,Sice_bulk,Ff,Fs, &
                                 higs,higb,himb,hims,himi)
  
! !USES:

   IMPLICIT NONE

   !INPUT PARAMETERS:
      
      integer, intent(in)  :: nilay
      real(rk), intent(in) :: rhow,rhof,dti
      real(rk), intent(in) :: rhoice,rhosnow
      real(rk), intent(in) :: SSS
      real(rk), intent(in) :: Tice(nilay+1)
      real(rk), intent(in) :: zi(nilay+1)
      real(rk), intent(in) :: simasso,snmasso,meltmasso
      real(rk), intent(inout) :: Sice_bulk
      real(rk), intent(out) :: Fs,Ff  
! Partitionated growth/melt ice (thicknesses time-step diff in m):
! higs - thickness of ice grown at surface (snow-ice formation, from TopGrowth)
! higb - thickness of ice grown at bottom (abblation, from Botgrowth)
! hims - thickness of ice melt surface (snow+ice, from TopMelt)
! himb - thickness of ice melt bottom (from BotMelt)
! himi - thickness of ice melt internal (internal ice melting, from TerMelt)
      real(rk), intent(in) :: higs,higb,himb,hims,himi
         
! !LOCAL VARIABLES:
      logical               ::   first=.true.
      integer               ::   k
      real(rk)             ::   V,f=0,ratei,rates,alpha
! Salinity local variables [ppt]:
! So -  old timestep bulk salinity
! Ssi - salinity in the ice formed from snow depression
! Sg  - equilibrium salinity for winter gravity drainage
! Sfl - equilibrium salinity for summer flushing
! S1  - first salinity transition
! S2  - second salinity transition
! Szero - linear salinity profile
!
      real(rk)                    ::   So,Ssi
      real(rk),dimension(nilay+1) ::   Szero, Sice
      real(rk), parameter          ::   Sg=5.D+00,Sfl=2.D+00
      real(rk),parameter          ::   S1=3.5D+00,S2=4.5D+00
! Tg, Tfl - V09 time scales for gravity drainage (20 days) and
! summer flushing (10 days) [sec]
!
      real(rk), parameter    :: Tg=1728000,Tfl=864000

! Terms of the salinity conservation equation [ppt]
!
      real(rk)             ::   x1,x2,x3,x4
   
! Terms for the equivalent salt flux equation
      real(rk)              ::   xx1,xx2,xx3,xx4
   
! Terms for salt and freshwater surface flux
! Ff - freshwater flux [kg m-2 s-1]
! Fs - Total surface salt flux (Fs=Fb+Feq) [ppt m s-1 kg m-3] 
! Fb - brine drainage []
! Feq - equivalent salt flux due to growth/melt sea ice processes [ppt m s-1 kg m-3]
      real(rk)              ::   Fb,Feq


end subroutine saltice_prof_simple

!-----------------------------------------------------------------------

! call  albedo_ice_uvic(fluxt,I_0,PenSW,alb,ice_hs,ice_hi,Ts)
subroutine albedo_ice_uvic(fluxt,I_0,PenSW,alb,ice_hs,ice_hi,TTss)
! !USES:

      IMPLICIT NONE
!
! !INPUT PARAMETERS:
      real(rk), intent(in)      :: I_0,ice_hs,ice_hi,TTss
! !OUTPUT PARAMETERS:
      real(rk), intent(out)      :: alb,PenSW
! !INPUT/OUTPUT PARAMETERS:
      real(rk), intent(inout)      :: fluxt
!
! !LOCAL VARIABLES:
! snow_dist
      real(rk)                  :: I_0_tmp
!

end subroutine albedo_ice_uvic


!-----------------------------------------------------------------------


!-----------------------------------------------------------------------

!BOP
!
! !IROUTINE: Cleaning up the mean flow variables
!
! !INTERFACE:
   subroutine clean_ice_uvic()  !still need to initialize/call this subroutine from outside this file
!
! !DESCRIPTION:
!  De-allocates all memory allocated via init\_icemodel()
!
! !USES:
!  use ncdfout, only: close_ncdf

      IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): 
!
!  See log for the icemodel module
!
!EOP
!-----------------------------------------------------------------------
!BOC
   !LEVEL1 'clean_ice_uvic'  !commented for now 

   !LEVEL3 'closing ice.nc file...'  !commented for now 
!   call close_ncdf()
   
      return

   end subroutine clean_ice_uvic
!EOC
!----------------------------------------------------------------------
! end internal subroutines
!----------------------------------------------------------------------     
!
! 
!
!
!
!-----------------------------------------------------------------------
! begin internal functions
!-----------------------------------------------------------------------
!BOF
! !IROUTINE: inverse complementary error function
!
! !INTERFACE:
real function erfc(x)

!DESCRIPTION: Note: for FORTRAN2008 and higher ERFC is an 
!                intrinsic Fortran function!

!USES:
     implicit none

!INPUT PARAMETERS:
    real(rk), intent(in) :: x 

!LOCAL VARIABLES:
    real(rk)            :: t,z
! !REVISION HISTORY:
!  Original author(s): 
!
!  See log for the icemodel module
!     
!-----------------------------------------------------------------------
!BOC
      z = abs(x)
      t = 1.0 / ( 1.0 + 0.5 * z )

      erfc =  t * exp( -z * z - 1.26551223 + t *               &
             ( 1.00002368 + t * ( 0.37409196 + t *             &
             ( 0.09678418 + t * (-0.18628806 + t *             &
             ( 0.27886807 + t * (-1.13520398 + t *             &
             ( 1.48851587 + t * (-0.82215223 + t * 0.17087277 )))))))))

      if ( x.lt.0.0 ) erfc = 2.0 - erfc

      print *, 'erfc', erfc
      return
end function erfc 

!EOF
!----------------------------------------------------------------------

!BOF
! !IROUTINE: inverse complementary error function
!
! !INTERFACE:
real function inverfc(y)

!DESCRIPTION:

!USES:
     implicit none

!INPUT PARAMETERS:
    real(rk), intent(in) :: y 

!LOCAL VARIABLES:
    real(rk), parameter :: qa  =  9.16461398268964d-01, & 
                         & qb  =  2.31729200323405d-01, &
                         & qc  =  4.88826640273108d-01, &
                         & qd  =  1.24610454613712d-01, &
                         & q0  =  4.99999303439796d-01, &
                         & q1  =  1.16065025341614d-01, &
                         & q2  =  1.50689047360223d-01, &
                         & q3  =  2.69999308670029d-01, &
                         & q4  = -7.28846765585675d-02, &

                         & pa  =  3.97886080735226000d+00, &
                         & pb  =  1.20782237635245222d-01, &
                         & p0  =  2.44044510593190935d-01, &
                         & p1  =  4.34397492331430115d-01, &
                         & p2  =  6.86265948274097816d-01, &
                         & p3  =  9.56464974744799006d-01, &
                         & p4  =  1.16374581931560831d+00, &
                         & p5  =  1.21448730779995237d+00, &
                         & p6  =  1.05375024970847138d+00, &
                         & p7  =  7.13657635868730364d-01, &
                         & p8  =  3.16847638520135944d-01, &
                         & p9  =  1.47297938331485121d-02, &
                         & p10 = -1.05872177941595488d-01, &
                         & p11 = -7.43424357241784861d-02, &

                         & p12 =  2.20995927012179067d-03, &
                         & p13 =  3.46494207789099922d-02, &
                         & p14 =  1.42961988697898018d-02, &
                         & p15 = -1.18598117047771104d-02, & 
                         & p16 = -1.12749169332504870d-02, &
                         & p17 =  3.39721910367775861d-03, &
                         & p18 =  6.85649426074558612d-03, & 
                         & p19 = -7.71708358954120939d-04, & 
                         & p20 = -3.51287146129100025d-03, &
                         & p21 =  1.05739299623423047d-04, &
                         & p22 =  1.12648096188977922d-03
                    

     real(rk)           :: z, w, u, s, t 
! !REVISION HISTORY:
!  Original author(s): 
!
!  See log for the icemodel module
!     
!-----------------------------------------------------------------------
!BOC
         z = y
         if (y .gt. 1d0) then
         z = 2d0 - y
         end if
         w = qa - log(z)
         u = sqrt(w)
         s = (qc + log(u)) / w
         t = 1d0 / (u + qb)
         inverfc = u * (1d0 - s * (0.5d0 + s * qd)) - &
            ((((q4 * t + q3) * t + q2) * t + q1) * t + q0) * t
         t = pa / (pa + inverfc)
         u = t - 0.5d0
         s = (((((((((p22*u+p21)*u+p20)*u+p19)*u+p18)*u+p17)*u+p16)*u+ &
            p15)*u+p14)*u+p13)*u+p12
         s = ((((((((((((s*u+p11)*u+p10)*u+p9)*u+p8)*u+p7)*u+ &
            p6)*u+p5)*u+p4)*u+p3)*u+p2)*u+p1)*u+p0)*t- &
            z*exp(inverfc*inverfc-pb)
         inverfc = inverfc + s * (1d0 + inverfc * s)
         if (y .gt. 1d0) then
         inverfc = -inverfc
         end if

             return
end function inverfc 

!EOF
!-----------------------------------------------------------------------

!BOF
! !IROUTINE: 
!
! !INTERFACE:
real function W(y)
!
!DESCRIPTION: W???
!
!USES: 
      implicit none
!INPUT PARAMETERS:
      real(rk), intent(in) :: y
!LOCAL VARIABLES:
      real(rk), parameter  :: M1=0.3361D+00, M2=-0.0042D+00 
      real(rk), parameter  :: M3=-0.0201D+00
      real(rk)             :: sigma
! !REVISION HISTORY:
!  Original author(s): 
!
!  See log for the icemodel module
!
!-----------------------------------------------------------------------
!BOC
      sigma=-1.D+00-log(-y)

      W=-1.D+00-sigma-2.D+00/M1*(1.D+00-1.D+00/(1.D+00+ &
      ((M1*sqrt(sigma/2.D+00))/(1.D+00+M2*sigma*exp(M3*sqrt(sigma))))))

end function W
!EOF
!----------------------------------------------------------------------

real function swkappai(Tin)
!DESCRIPTION:

!USES:
      implicit none
!INPUT PARAMETERS:
      real(rk), intent(in) :: Tin
! !REVISION HISTORY:
!  Original author(s): 
!
!  See log for the icemodel module
!
!-----------------------------------------------------------------------
!BOC
      swkappai=((swkappaim+swkappaif)+(swkappaim-swkappaif)* &
      tanh(Tin-273.15D+00))/2.D+00
!test
!             swkappai=1.5           

end function swkappai
!EOF
!----------------------------------------------------------------------

!BOF
! !IROUTINE: snow transmissivity
!
! !INTERFACE:

real function transs(Tin)
!DESCRIPTION:

!USES:
      implicit none
!INPUT PARAMETERS:
      real(rk), intent(in) :: Tin
! !REVISION HISTORY:
!  Original author(s): 
!
!  See log for the icemodel module
!             
!-----------------------------------------------------------------------
!BOC
      transs=((transsm+transsf)+(transsm-transsf)* &
      tanh(Tin-273.15D+00))/2.D+00
   !test
!            transs=0.05

end function transs
!EOF
!----------------------------------------------------------------------
!BOF
! !IROUTINE: ice transmissivity
!
! !INTERFACE:

real function transi(Tin)
!DESCRIPTION:

!USES:

!INPUT PARAMETERS:
      real(rk), intent(in) :: Tin
             
!LOCAL VARIABLES:
! !REVISION HISTORY:
!  Original author(s): 
!
!  See log for the icemodel module
!
!-----------------------------------------------------------------------
!BOC

!H!begin
!             transi=0.5
      transi=((transim+transif)+(transim-transif)*tanh(Tin-273.15D+00))/2.D+00
!H!end

end function transi
!EOF


   end module stim_flato

!-----------------------------------------------------------------------
! Copyright by the GETM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
