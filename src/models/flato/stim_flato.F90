!-----------------------------------------------------------------------
!BOP
!
! !MODULE: stim_flato --- flato thermodynamic ice model
! \label{sec:stim_flato}
!
! !INTERFACE
module stim_flato
!
! !DESCRIPTION:
!
!
!
!
!
! !USES:
   use stim_variables
  
   IMPLICIT NONE
!
!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public                              :: init_stim_flato
   public                              :: do_ice_uvic 
!
! !PRIVATE DATA MEMBERS:
   !real(rk), pointer :: ice_hi,ice_hs  !jpnote 
!
!
! LOCAL VARIABLES: 
!
!   rhosnow     - snow density (kg m-3)
   real(rk), public :: rhosnow

!fixed size variables
!
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
!*****fluxes jpnote: added public to vars registered in registered_all_variables.F90
  ! qb           - long wave back radiation (in-out) (W m-2)	
   real(rk), public :: qb
   !qh           - latent heat flux into ice (W m-2)		
   real(rk) :: qh
  ! qe           - sensible heat flux into ice (W m-2)		
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
!   Ts           - upper surface temperature (K) ! jpnote - already declared for use in winton
   real(rk) :: Ts_uvic
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
subroutine init_stim_flato() 
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
! 
! !LOCAL VARIABLES:
!
   integer             :: k,rc     !from init_ice_uvic - jp 
! !LOCAL PARAMETERS:

   !ice_hs => Hsnow
   !ice_hi => Hice !jpnote 
!EOP
!-----------------------------------------------------------------------
!BOC
!
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
   print *, 'Ts_uvic',  Ts_uvic
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
!
! !ROUTINE: Calculate ice thermodynamics \label{sec:do_ice_uvic}
!
! !INTERFACE: 
subroutine do_ice_uvic(dto,h,julianday,secondsofday,lon,lat, &
                        I_0,airt,airp,rh,u10,v10,precip,cloud, &
                        TSS,SSS,rhowater,rho_0, &
                        longwave_radiation_method,hum_method,fluxes_method, &
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
   integer, intent(in)       :: longwave_radiation_method ! method for LW   !read in from namelist in airsea --> defined as a local variable in airsea
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
   !print *,'in do_ice_uvic ice_uvic_zi',ice_uvic_zi
   !for testing jpnote 
   !print *,'Cond(nilay)',Cond(nilay)
   !print *,'rhoCp(nilay)',rhoCp(nilay)
   !print *,'Sint(nilay+1) ',Sint(nilay+1) 
  ! print *,'dzi(nilay)',dzi(nilay)
#if 0
   print *,' dto ', dto! ocean timestep (sec) ! dto ??? 
   !real(rk), intent(in)     :: dz !--> h from meanflow ???
   print *,' h ', h ! sea surface layer thickness --> is this the same as dz ???
   
   print *,'julianday',julianday ! this julian day --> from time
   print *,'secondsofday',secondsofday ! seconds for this day -->  from time
   print *,'lon',  lon   ! longitude for this point --> longitude from gotm
   print *,'lat ',lat! latitude for this point --> latitude from gotm 
   print *,'I_0',  I_0  ! shortwave radiation at sea surface  
   print *,' airt', airt ! 2m temperature
   print *,'airp', airp  ! sea surface pressure
   !real(rk), intent(in)     :: hum   ! relative humidity from airsea
   print *,' rh', rh    ! relative humidity --> is this the same as hum ??? 
   print *,' u10 ', u10  ! 10 m wind u-component
   print *,'v10',  v10  ! 10 m wind v-component
   print *,' precip',precip! freshwater precipitation (m/s)
   print *,' cloud',cloud ! cloud cover
   print *,' TSS ',  TSS   ! sea surface temperature
   print *,' SSS ',  SSS  ! sea surface salinity
   print *,'rhowater', rhowater  ! sea surface layer density --> called rho(nlev) in new code
   print *,' rho_0 ', rho_0 ! reference density --> from meanflow
   print *,' longwave_radiation_method ',longwave_radiation_method! method for LW   !read in from namelist in airsea --> defined as a local variable in airsea
   print *,' hum_method', hum_method ! method for humidity
   print *,' fluxes_method', fluxes_method ! method for fluxes


   print *,'ice_hi',ice_hi
   print *,'ice_hs',ice_hs
   print *,'hm',ice_uvic_hm
   print *,'Tice(nilay+1)',Tice(nilay+1)
   print *,'Cond(nilay)',Cond(nilay)
   print *,'rhoCp(nilay)',rhoCp(nilay)
   print *,'Sint(nilay+1) ',Sint(nilay+1) 
   print *,'dzi(nilay)',dzi(nilay)
   print *,'zi(nlmax)',zi(nlmax)
   print *, 'Pari(nilay+1)',Pari(nilay+1)
   print *,'Told(nilay+1)',Told(nilay+1)
   print *,'Fh',Fh
   print *,'Ff',Ff
   print *,'Fs',Fs
   print *,'Sicebulk',Sice_bulk
   print *,'TopMelt',TopMelt
   print *,'BotMelt', BotMelt
   print *, 'TerMelt',TerMelt
   print *,'TopGrowth',TopGrowth
   print *,'BotGrowth',BotGrowth
   print *,'Hmix',Hmix
   print *,'Aice',Aice_i
   print *,'Asnow',Asnow_i
   print *,'Amelt',Amelt_i
   print *,'swr_0',swr_0
   print *,'precip_i',precip_i
   print *,'sfall_i',sfall_i
   print *, '-----------------------------------------------------------------------'
   print *, '-----------------------------------------------------------------------'
!-----------------------------------------------------------------------
#endif
!

! ------------------------------------ jpnote for testing 
!call growthtb(rhowater,nilay,rhoCp,dzi,Tice,TopGrowth,TerMelt,TopMelt,&
!BotGrowth,BotMelt,Hmix,ohflux)
! ------------------------------------
!#if 0
!BOC
!   LEVEL0 'do_ice_uvic'
!  Calculate seawater freezing temperature
!   Tfreezi = (-0.0575D00*SSS)+kelvin
!KB   STDERR T,S
!
! set ice timestep to ocean timestep

dti=dto
! initialize new ice formation in open water, melt/growth for this timestep
      dmsi=0 
      TopMelt = 0 
      Botmelt = 0 
      TerMelt = 0 
      TopGrowth = 0 
      BotGrowth = 0 
      Ff= 0
      Fs= 0
! degC to K unit conversion
   airtk=airt+kelvin

 if(sfall_method .eq. 2) then
!get snowfall from precipitation
         if((airtk-Tmelts).lt.(-5.0)) then
                 p_tmp = 1.0
         elseif((airtk-Tmelts).ge.(-5.0).and.(airtk-Tmelts).le.(5.0)) then
                     p_tmp = 1.0-(airtk-Tmelts+5.0)*0.1
         else
                     p_tmp = 0.0
         endif
         ! - add factor for drifting snow...
         p_tmp=p_tmp*dfact
      sfall=precip*p_tmp*rhowaterfresh/rhoscold ! m/s  
 endif
! set old snow and ice mass and meltpond mass

!nsnote, are all the mass variables known from the last timestep? if not, set here...???
      snmasso=snmass
      simasso=simass
      meltmasso=meltmass

!  add snowfall to existing snow if ice is present (no accumulation over open water!) 

!  and if air temp. is below freezing.
      if (simass.gt.rhoice*hsmin .and. airtk.lt. kelvin) then               
         snmass=snmass+rhoscold*sfall*dti 
!  Remove snow falling on ice from precip
       precip_i=precip !H! Save precip before set to zero.
       precip=precip-sfall*rhoscold/rhowaterfresh
       precip =max(0.0D0,precip)
      !NSnote, make sure no double counting of precip
      ! if no ice then precip  (or if precip not given, water equivalent for 
      ! snow fall) should go in the ocean, if ice then sfall, but no 
      ! additional precip in the ocean !!!
      endif

!  calculate snow density - depends on mean snow layer temperature

   Tsav=(Tice(1)+Tice(2))/2.D+00
   if(Tsav.lt.Tmelti) then
      rhosnow=rhoscold
   else
      rhosnow=rhoswarm
   endif


   swr_0 = I_0/(1-alb) !H! save incidental swr.
!  if ice is present, calculate growth or melt, if not calculate potential new growth
!
   if(simass.eq.0.D+00) then
      Hmix=Hmix+(heat+I_0)*dti
      if(depmix.eq.0.0) depmix=h 
      call open_water(nilay,I_0,Sice_bulk,Hmix,Tice,depmix,TSS,Fh,heat,precip,precip_i)
   else
! reset I_0 so it can be recalculated for ice
      I_0 = I_0/(1-alb) ! Here alb is water albedo from GOTM

!  do surface energy budget and heat conduction calculations via N-R 
!  iteration on surface temperature
!  This will update the ice temperature at each layer

      call nr_iterate(hum_method,longwave_radiation_method,fluxes_method,nilay,&
                      airt,rh,cloud,I_0,Told,Tice,Pari,&
                      Sice_bulk,ice_hi,ice_hs,dzi,Cond,rhoCp,zi,Sint,&
                      lat,u10,v10,precip,airp,evap,alb)

!  construct depth array by adding up dzi's
      zi(1)= 0
      do k=2,nilay+1
         zi(k)=zi(k-1)+dzi(k-1)
      enddo

!   add short wave radiation that penetrates snow/ice slab to mixed layer 
!nsnote, check if ok that if snow_dist can be removed and calculation done with PAR for all cases
!if (snow_dist) then 
       Hmix=Hmix+Pari(size(Pari))*dti
!else
!       Hmix=Hmix+PenSW*exp(-swkappa*zi(nilay+1))*dti
!endif

!add mixed layer heat storage to oceanic heat flux -

       ohflux=Hmix/dti
       Hmix=0.D+00

!...Calculate sea-ice turbulent heat flux and input it to the ohflux
       Fh=ohflux !NSnote check
       if(ice_hi .gt. 0) then  
!          heat=-Fh !NSnote, don't think so, if, ohflux is used for ice growth
!          heat=0.D00 
       end if ! NSnote: What if no ice?, all melted..what heat then? 
!         I_0=PenSW  
         I_0=Pari(nilay+1) !NScheck

!...........calculate growth or melt amount at top and bottom surface, 
!...........converting snow to ice if ice surface is submerged

!  initialize top, interior and bottom melt and growth  
!
    


      call growthtb(rhowater,nilay,rhoCp,dzi,Tice,TopGrowth,TerMelt,TopMelt,&
                       BotGrowth,BotMelt,Hmix,ohflux)
!NSnote, Hmix contains extra heat left from melting, 
!this should be transferred back into ocean to heat the water
       heat=Hmix/dti

    if (ice_salt) &
      call saltice_prof_simple(dti,nilay,SSS,simasso,snmasso,meltmasso, &
                    rhoice,rhosnow,rhowater,rhowaterfresh,zi,Tice,Sice_bulk,Ff,Fs, &
                     TopGrowth,BotGrowth,TopMelt,BotMelt,TerMelt)


   endif

!...Update ice and snow thickness from mass and density so accurate for end of timestep
!
       ice_hs=snmass/rhosnow
       ice_hi=simass/rhoice
       ice_hm=meltmass/rhowaterfresh

!  Change units from m (per timestep) to m/s for potential transfer to FABM 
!  (ice algae, since FABM is not aware of the timestep) 
                    
       TopMelt=TopMelt/dti
       BotMelt=BotMelt/dti
       TerMelt=TerMelt/dti
       TopGrowth=TopGrowth/dti
       BotGrowth=BotGrowth/dti

! set area fractions of ice so they can be written out
       Aice_i=Aice
       Asnow_i=Asnow
       Amelt_i=Amelt
       sfall_i=sfall

   return

!#endif

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


 !-----------------------------------------------------------------------



 !      LEVEL1'therm1d'

!
!...Calculate ice and snow thickness from mass and density
!...(minimum thickness to avoid dividing by zero elsewhere in code)
!
   ice_hs=max(hzero*0.1D+00,snmass/rhosnow)
   ice_hi=max(hzero,simass/rhoice)

!...snow_dist: max. snow height within a year to get the ratio of current snow
!    depth to maximum snow depth for the Weibull-distributed snow 
!      if(snow_dist .and. snmass .eq. 0.D+00) then
   if(snmass .eq. 0.D+00) then
      hsmax=hsmin
!      else if (snow_dist .and. snmass .gt. 0.D+00) then
   else if (snmass .gt. 0.D+00) then
      hsmax=max(ice_hs,hsmax)
   end if
!...end snow_dist
!
!...Calculate array of layer thicknesses from ice and snow thickness
!...and calculate conductivity and heat capacity arrays. Parameterizations
!...based on Ebert and Curry (1993)
!
!...Total slab thickness
!
   hlay=ice_hi+ice_hs
!
!...Check if sufficient snow to warrant snow layers
!
   if(ice_hs.gt.hsmin) then
!
!.......estimate layer thickness (i.e. if evenly-spaced)
      thicklay=hlay/float(nilay)
!
!.......number of snow layers
!         nslay=max(1,nint(ice_hs/thicklay))
!         nslay=min(nslay,nilay-1)

!***Temporary fix***

! --- ONLY ONE SNOW LAYER USED
! --- Do not change this until a conservative method for reapportioning
! --- temperature when number of snow layers changes is developed
! --- G. Flato 14/Sept/94

   nslay=1

!***Temporary fix***
!
!.......thickness of snow layers
      snlaythick=ice_hs/float(nslay)
!
!.......loop over snow layers
      do k=1,nslay
        dzi(k)=snlaythick
!.........conductivity of snow, including parameterization of water vapor
!.........diffusion
        Ticeav=(Tice(k)+Tice(k+1))/2.D+00
        Cond(k)=(2.845D-06*rhosnow**2)+(2.7D-04*2.D+00**((Ticeav-233.D+00)/5.D+00))
!....snow_dist
        if (snow_dist .and. distr_type .eq. 0 .and. ice_hs .ge. hsmax) then
            Cond(k)=Cond(k)*pi/4.D+00
        end if

!end snow_dist
!.........heat capacity of snow 
        rhoCp(k)=rhosnow*(92.88D+00+7.364D+00*Ticeav)
      enddo
!.......reset 'hlay' to be thickness spanned by ice layers
      hlay=ice_hi
   else
!....if not, number of snow layers equals zero
      nslay=0
   endif
!
!...Now do ice layers
!
   do k=nslay+1,nilay
!.......ice layer thickness
      dzi(k)=hlay/float(nilay-nslay)
!.......salinity of ice
 if (ice_salt) then 
        Sice_bulk=Sice_bulk 
      else
        if(ice_hi.ge.0.57D+00) then
          Sice_bulk=3.2
        else
          Sice_bulk=(14.24D+00-19.39D+00*ice_hi)
        endif
      endif
      Ticeav=min(Tmelti-0.1D+00,(Tice(k)+Tice(k+1))/2.D+00)
!.......conductivity of ice
      Cond(k)=Condfi+0.1172D+00*Sice_bulk/(Ticeav-Tmelti)
!.......minimum value to avoid zero or negative conductivities
      Cond(k)=max(Cond(k),Condfi/5.D+00)
!.......heat capacity of ice
      rhoCp(k)=rhoCpfi+1.715D+07*Sice_bulk/(Ticeav-Tmelti)**2
!.......limit value to 100 times fresh water value to avoid
!.......extreme values when temperature nears melting
      rhoCp(k)=min(rhoCpfi*100.D+00,rhoCp(k))
   enddo
!
!...Calculate internal sources due to penetrating shortwave radiation
!...Note, internal sources specified at layer interfaces
!
   zi(1)=0.D+00
   zi(2)=zi(1)+dzi(1)/2.D+00
   zi(3)=zi(2)+dzi(1)/2.D+00+dzi(2)/2.D+00
   Sint(1)=0.D+00            !internal sources must be zero at
   Sint(nilay+1)=0.D+00       !top and bottom surfaces
   Pari(1)=PenSW !Pari at top of snow set to PenSW
!....if no ice, put all short wave radiation into surface flux
! NSnote: not sure why this is here, since this routine should not 
! be called if no ice

   if(hlay.lt.hzero) then
     fluxt=fluxt+PenSW
   else
!....snow_dist
     if (snow_dist .and. distr_type .eq. 0 .and. ice_hs .ge. hsmax) then 
! absorption at surface, use PenSW from albedo_ice since Anow=1.0
        fluxt=fluxt+PenSW*swkappas(Tsav)*(1.D+00-swkappas(Tsav)*ice_hs *erfc(swkappas(Tsav)&
              *ice_hs/sqrt(pi))*exp((swkappas(Tsav)*ice_hs)**2.D+00/pi))
        Pari(2)=PenSW*(1-swkappas(Tsav)*ice_hs*erfc(swkappas(Tsav)*ice_hs/sqrt(pi))&
                *exp((swkappas(Tsav)*ice_hs)**2.D+00/pi))
     else if(snow_dist .and. distr_type .eq. 0 .and. ice_hs .gt. 0.D+00 &
             .and. ice_hs .lt. hsmax .and. Asnow.gt.0.D+00) then
! absortion at surface redone here to accurately account for different surface types
        fluxt=fluxt+Asnow*I_0*(1.D+00-albsnow)*transs(Tsav)*swkappas(Tsav) & 
              *(exp(-(inverfc(ice_hs/hsmax))**2.D+00-swkappas(Tsav)*2.D+00 &
              *hsmax/sqrt(pi)*inverfc(ice_hs/hsmax))&
              -swkappas(Tsav)*hsmax*erfc(swkappas(Tsav)*hsmax/sqrt(pi) &
              +inverfc(ice_hs/hsmax))*exp((swkappas(Tsav)*hsmax)**2.D+00/pi)) &
              +Aice*I_0*(1.D+00-albice)*transi(Tsav)

        Pari(2)=Asnow*I_0*(1.D+00-albsnow)*transs(Tsav) &
                *(exp(-(inverfc(ice_hs/hsmax))**2.D+00-swkappas(Tsav)*2.D+00 &
                *hsmax/sqrt(pi)*inverfc(ice_hs/hsmax))&
                -swkappas(Tsav)*hsmax*erfc(swkappas(Tsav)*hsmax/sqrt(pi) &
                +inverfc(ice_hs/hsmax))*exp((swkappas(Tsav)*hsmax)**2.D+00/pi)) &
                +Aice*I_0*(1.D+00-albice)*transi(Tsav)
     else if(snow_dist .and. distr_type .eq. 1 .and. ice_hs .ge. hsmax) then
! absorption at surface, use PenSW from albedo_ice since Anow=1.0
        fluxt=fluxt+PenSW*swkappas(Tsav)*4.D+00/ice_hs**2.D+00*1.D+00 &
              /(2.D+00/ice_hs+swkappas(Tsav))**2
        Pari(2)=PenSW*4.D+00/ice_hs**2.D+00*1.D+00/(2.D+00/ice_hs+swkappas(Tsav))**2
     else if(snow_dist .and. distr_type .eq. 1 .and. ice_hs .gt. 0.D+00 .and. ice_hs .lt. hsmax) then
! absortion at surface redone, here to accurately account for different surface types
        fluxt=fluxt+Asnow*I_0*(1.D+00-albsnow)*transs(Tsav)*swkappas(Tsav) &
              *4.D+00/hsmax**2.D+00*1.D+00/(2.D+00/hsmax+swkappas(Tsav)) &
              *exp(hsmax/2.D+00*(W(-2.D+00*ice_hs/hsmax*exp(-2.D+00))+2.D+00) &
              *(2.D+00/hsmax+swkappas(Tsav)))*(hsmax/2.D+00*(-W(-2.D+00*ice_hs/hsmax &
              *exp(-2.D+00))-2.D+00)+1.D+00/(2.D+00/hsmax+swkappas(Tsav))) &
              + Aice*I_0*(1.D+00-albice)*transi(Tsav)       
        Pari(2)=Asnow*I_0*(1.D+00-albsnow)*transs(Tsav) &
              *4.D+00/hsmax**2.D+00*1.D+00/(2.D+00/hsmax+swkappas(Tsav)) &
              *exp(hsmax/2.D+00*(W(-2.D+00*ice_hs/hsmax*exp(-2.D+00))+2.D+00) &
              *(2.D+00/hsmax+swkappas(Tsav)))*(hsmax/2.D+00*(-W(-2.D+00*ice_hs/hsmax &
              *exp(-2.D+00))-2.D+00)+1.D+00/(2.D+00/hsmax+swkappas(Tsav))) &
              + Aice*I_0*(1.D+00-albice)*transi(Tsav)
     else if (.not. snow_dist .and. ice_hs .gt. 1.D-02) then
        fluxt=fluxt+PenSW*swkappas(Tsav)*(exp(-swkappas(Tsav)*zi(1)) &
              -exp(-swkappas(Tsav)*zi(2)))
        Pari(2)=PenSW*exp(-swkappas(Tsav)*ice_hs)
     else
        fluxt=fluxt+PenSW*swkappai(Tsav)*(exp(-swkappai(Tsav)*zi(1)) &
              -exp(-swkappai(Tsav)*zi(2)))
        Pari(2)=PenSW*exp(-swkappai(Tsav)*zi(2))
     end if
    
     if (meltpond) then
       !note for ice_hs ge hsmax Asnow is 1.0 and Amelt is 0.0, so no doublecounting
        fluxt=fluxt+Amelt*I_0*(1.D+00-albmelt)*transm
       Pari(2)=Pari(2)+Amelt*I_0*(1.D+00-albmelt)*transm
     end if

!......put short wave radiation absorbed in upper half of top layer
!......into surface warming
       Sint(2)=Pari(1)*swkappai(Tice(2))*(exp(-swkappai(Tice(2))*zi(2))- &
               exp(-swkappai(Tice(2))*zi(3)))/(0.5D+00*(dzi(1)+dzi(2)))

!......do ice layers 
    do k=3,nilay-1
       zi(k+1)=zi(k)+dzi(k-1)/2.D+00+dzi(k)/2.D+00
       
       Pari(k)=Pari(k-1)*exp(-swkappai(Tice(k))*dzi(k-1))
       Sint(k)=Pari(1)*swkappai(Tice(k))*(exp(-swkappai(Tice(k))*zi(k))- &
               exp(-swkappai(Tice(k))*zi(k+1)))/(0.5D+00*(dzi(k-1)+dzi(k)))
     enddo
     
     zi(nilay+1)=zi(nilay)+dzi(nilay-1)/2.D+00+dzi(nilay)
!        zi(nilay+1)=zi(nilay)+dzi(nilay-1)/2.D+00+dzi(nilay)/2.D+00 !NScheck
!......all of short wave radiation absorbed in bottom-most layer
!......applied at interface between bottom two layers for simplicity
     Pari(nilay)=Pari(nilay-1)*exp(-swkappai(Tice(nilay))*dzi(nilay-1))
     Sint(nilay)=Pari(1)*swkappai(Tice(k))*(exp(-swkappai(Tice(k))*zi(nilay)) &
                -exp(-swkappai(Tice(k))*zi(nilay+1)))/(0.5D+00*dzi(nilay-1)+dzi(nilay))
     Pari(nilay+1)=Pari(nilay)*exp(-swkappai(Tice(nilay+1))*dzi(nilay)) !NScheck
   endif
!...end snow_dist


!     
!...Set up boundary conditions for 1-D heat conduction equation
!
!....if upper surface is melting, use temp. boundary condition
   if(fluxt.ge.0.D+00.and.Tice(1).ge.Tmelts) then
     bctype(1)=-1.D+00       !indicates temp. boundary condition
     bcs(1)=Tmelts           !Note, under melting conditions, surface
!                                temperature is taken to be melting point
!                                of snow (fresh water)
   else
!
!....if upper surface is cooling, use flux boundary condition
     bctype(1)=1.D+00        !indicates flux boundary condition
     bcs(1)=fluxt            !Note, under cooling conditions, fluxt is
!                                the outgoing flux calculated in SEBUDGET
   endif
!
!....bottom boundary condition always freezing temp. of sea water
   bctype(2)=-1.D+00      !indicates lower boundary has temp specified
   bcs(2)=Tfreezi

!
!...Solve 1-D heat conduction equation to get new temperature profile
!...and conductive flux at top and bottom
!     
   if(hlay.ge.hlaymin) then
      call cndiffus(bctype,bcs,nilay,dzi,rhoCp,Cond,Sint,Tice)

!...Check if layer is thinner than HLAYMIN but greater than 0. If so,
!...assume linear temperature profile and calculate conduction from
!...constant conductivity rather than solving transient diffusion
!...equation. If ice and snow have vanished, calculate ocean surface
!...temperature assuming specified mixed-layer depth.
!
   elseif(hlay.lt.hlaymin .and. hlay.gt.0.D+00)then
!.......thin layer of ice: assume linear temp. profile and calculate
!.......new surface temperature
      Tbot=Tfreezi
      Ttop=(hlay/Condfi)*fluxt+Tbot
      if(Ttop.gt.Tmelts) then
         Ttop=Tmelts
         do k=1,nilay
            Tice(k)=Ttop
         enddo
         Tice(nilay+1)=Tbot
         Iceflux(1)=0.D+00
         Iceflux(2)=0.D+00
      else
         do k=1,nilay+1
            Tice(k)=Ttop+(Tbot-Ttop)*float(k-1)/float(nilay)
         enddo
         Iceflux(1)=fluxt
         Iceflux(2)=fluxt
      endif
   endif
!
!
!    LEVEL1'end therm1d'


return




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


 !-----------------------------------------------------------------------

!      LEVEL1'cndiffus'
!
!...check if nilay is greater than maximum allowed
!
   if(nilay.gt.nlmax) then
      print*,'!!!! FATAL ERROR IN CNDIFFUS - NILAY > NLMAX'
      stop
   endif
!
!...check if nilay is less than minimum allowed
   if(nilay.lt.2) then
      print*,'!!!! FATAL ERROR IN CNDIFFUS - NILAY < 2'
      stop
   endif
!
!...calculate average volumetric heat capacity at temperature grid points
!...from layer values
!
   do j=1,nilay-1
      rCpav(j+1)=(rhoCp(j)*dzi(j)+rhoCp(j+1)*dzi(j+1))/(dzi(j)+dzi(j+1))
   enddo
   rCpav(1)=rhoCp(1)
   rCpav(nilay+1)=rhoCp(nilay)
!
!...set up tridiagonal matrix coefficients for Crank-Nicholson scheme
!
   do l=2,nilay
      C(l,1)=-theta*(2.D+00*Cond(l-1)*dti/dzi(l-1))/(dzi(l)+dzi(l-1))/  &
            rCpav(l)
      C(l,2)=1.D+00+(theta*2.D+00*Cond(l-1)*dti/dzi(l-1)+theta*2.D+00*  &
            Cond(l)*dti/dzi(l))/(dzi(l)+dzi(l-1))/rCpav(l)
      C(l,3)=-theta*(2.D+00*Cond(l)*dti/dzi(l))/(dzi(l)+dzi(l-1))/      &
            rCpav(l)
   enddo
!
!...set up right hand side vector
!
   do l=2,nilay
      R(l)=(1.D+00-theta)*(((2.D+00*Cond(l)*dti/dzi(l))/(dzi(l)+ &
          dzi(l-1)))*(Tice(l+1)-Tice(l))-((2.D+00*Cond(l-1)*dti/dzi(l-1))/ &
          (dzi(l)+dzi(l-1)))*(Tice(l)-Tice(l-1)))/rCpav(l)+Tice(l)
!......add internal heat sources (if any)
      R(l)=R(l)+Sint(l)*dti/rCpav(l)
   enddo
!
!...add temperature boundary conditions (if any) to RHS 
!
!.....upper boundary
   if(bctype(1).lt.0.D+00) then
      C(1,1)=0.D+00
      C(1,2)=1.D+00
      C(1,3)=0.D+00
      R(1)=bcs(1)
   endif
!.....lower boundary
   if(bctype(2).lt.0.D+00) then
      C(nilay+1,1)=0.D+00
      C(nilay+1,2)=1.D+00
      C(nilay+1,3)=0.D+00
      R(nilay+1)=bcs(2)
   endif
!
!...add flux boundary conditions (if any) to RHS
!
!.....upper boundary
   if(bctype(1).gt.0.D+00) then
      C(1,1)=0.D+00
      C(1,2)=1.D+00+(2.D+00*dti*Cond(1)*theta/dzi(1)**2)/rCpav(1)
      C(1,3)=-(2.D+00*dti*Cond(1)*theta/dzi(1)**2)/rCpav(1)
      R(1)=Tice(1)+(Tice(2)-Tice(1))*(2.D+00*dti*Cond(1)*(1.D+00-theta)/ &
          dzi(1)**2)/rCpav(1)+2.D+00*bcs(1)*dti/rCpav(1)/dzi(1)
   endif
!.....lower boundary
   if(bctype(2).gt.0.D+00) then
      C(nilay+1,1)=-(2.D+00*dti*Cond(nilay)*theta/dzi(nilay)**2)/ &
                 rCpav(nilay+1)
      C(nilay+1,2)=1.D+00+(2.D+00*dti*Cond(nilay)*theta/dzi(nilay)**2)/ &
                 rCpav(nilay+1)
      C(nilay+1,3)=0.D+00
      R(nilay+1)=Tice(nilay+1)-((Tice(nilay+1)-Tice(nilay))*2.D+00*dti*Cond(nilay)* &
               (1.D+00-theta)/dzi(nilay)**2)/rCpav(nilay+1)-(2.D+00*   &
               bcs(2)*dti/rCpav(nilay+1)/dzi(nilay))
   endif
!
!...save temperatures from previous time step for use in flux
!...calculation later
!
   do l=1,nilay+1
      To(l)=Tice(l)
   enddo
!
!...solve system of equations using efficient tridiagonal
!...matrix solver. Subroutine returns new values of temperature
!...in array 'Tice' - overwriting previous values.
!
   call trisol(C,R,nilay+1,Tice)
!
!...calculate fluxes at boundaries where temperature boundary 
!...conditions were specified. Flux returned in array 'flux' -
!...overwriting previous value.
!
!.....upper boundary
   Iceflux(1)=-Cond(1)*((1.D+00-theta)*(To(2)-To(1))+theta* &
           (Tice(2)-Tice(1)))/dzi(1)
!.....lower boundary
   Iceflux(2)=-Cond(nilay)*((1.D+00-theta)*(To(nilay+1)-To(nilay)) &
          +theta*(Tice(nilay+1)-Tice(nilay)))/dzi(nilay)
!

!
!    LEVEL1'end cndiffus'

return

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

   !-----------------------------------------------------------------------

!      LEVEL1'trisol'

!...dimension arrays
!
!       dimension Tice(nilay+1),C(nmax,3),R(nmax),temp(nmax)
!
!...check that 'n' is not greater than nmax
!
   if(n.gt.nmax) then
      print*,'!!!! FATAL ERROR IN TRISOL - N > NMAX'
      stop
   endif
!
!...check that first diagonal element is not too small
!
   if(C(1,2).lt.eps) then
      print*,'!!!! FATAL ERROR IN TRISOL '
      print*,'        - FIRST DIAGONAL ELEMENT < EPS'
      stop
   endif
!
!..."decomposition" phase of solution
!
   D1=C(1,2)
   Tice(1)=R(1)/D1
   do j=2,n
      temp(j)=C(j-1,3)/D1
      D1=C(j,2)-C(j,1)*temp(j)
      if(D1.lt.eps) then
         print*,'!!!! FATAL ERROR IN TRISOL '
         print*,'        - DIVISOR, D1 < EPS'
         print*,'        - Try making hlaymin larger'
         stop
      endif
      Tice(j)=(R(j)-C(j,1)*Tice(j-1))/D1
   enddo           
!
!..."back-substitution" phase of solution
!
   do j=n-1,1,-1
      Tice(j)=Tice(j)-temp(j+1)*Tice(j+1)
   enddo
!
!
!    LEVEL1'end trisol'

return


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


 !-----------------------------------------------------------------------

!      LEVEL1'growthtb'
!
!...Calculate surface melt, if any
!
   call surfmelt(Tice(1),TopMelt)

   !
   !...Calculate change in thickness due to bottom growth or melt
   !
   
         dhib=-(Iceflux(2)+ohflux)*dti/(rhoice*Hfi)
         simass=simass+dhib*rhoice
         if(dhib.gt.0.D+00) BotGrowth=BotGrowth+dhib
         if(dhib.lt.0.D+00) BotMelt=BotMelt-dhib
   
   !
   !...Check if sufficient snow to depress ice surface below sea level, if
   !...so, convert some snow to ice ("flooding"). 
   !...Add latent heat released by sea water
   !...freezing in snow pores to upper ice layer. (Note, if this heat raises 
   !...temperature above freezing, some ice will subsequently be melted by 
   !...the calculations which follow
   !
   
   
         snsub=snmass-(rhowater-rhoice)*max(0.D+00,simass)/rhoice
         if(snsub.gt.0.D+00) then
   !......limit submerged snow to that available
           if(snmass-snsub.lt.0.D+00) then
             snsub=snmass
             snmass=0.D+00
           else
             snmass=snmass-snsub
           endif
   
   !......ice mass increase enhanced due to sea water filling pores in snow
   !......and freezing
           simass=simass+snsub*(1.D+00+(rhoice-rhosnow)/rhoice)
           TopGrowth=TopGrowth+snsub*(1.D+00+(rhoice-rhosnow)/ &
                                                   rhoice)/rhoice
   !......add latent heat released by freezing sea water to ice (this may
   !......subsequently cause some melt if temp. is already near melting)
           porewat=snsub*(rhoice-rhosnow)/rhoice
           Tice(nslay+1)=Tice(nslay+1)+porewat*Hfi/(0.5D+00*(rhoCp(nslay)* &
                                 dzi(nslay)+rhoCp(nslay+1)*dzi(nslay+1)))
           
         endif
   !
   !...Calculate melt if temperature has risen above melting point). Note: 
   !...don't need to do bottom temperatures because it is always at the 
   !...freezing temp. of sea water. However, have to do surface temperature 
   !...in case conductive flux causes melting even when surface flux is 
   !...outward.
   !
   !....if there is snow cover, check for snow melt
         if(nslay.gt.0) then
   !......first do surface layer
           Ts=Tice(1) !nsnote this is not in CA code?
                 if(Ts.gt.Tmelts) then
   !...snow_dist
                     if(snow_dist) then
                     dhst=Asnow*(Ts-Tmelts)*rhoCp(1)*dzi(1)*0.5D+00/(Hfw*rhosnow)
                     dhi=(Aice+Amelt)*(Ts-Tmelts)*rhoCp(1)*dzi(2)*0.5D+00/(Hfi*rhoice)
                     simass=simass-dhi*rhoice 
                     Topmelt=Topmelt+dhi
                     else
                     dhst=(Ts-Tmelts)*rhoCp(1)*dzi(1)*0.5D+00/(Hfw*rhosnow)
                     end if
   !...end snow_dist
                     Ts=Tmelts
   !........reduce snow amount accordingly and replace heat that went unused
                     snmass=snmass-dhst*rhosnow
                     if(snmass.lt.0.D+00) then
                        Ts=Ts-snmass*Hfw/(rhoCp(1)*dzi(1)*0.5D+00)
                        snmass=0.D+00 
                        nslay=0
                     endif
                     Tice(1)=Ts
                 endif
   !
   !
   !......now check internal points 
           if(nslay.gt.1) then     
             do k=2,nslay
               Tl=Tice(k)
               if(Tl.gt.Tmelts) then
                 dhs=(Tl-Tmelts)*0.5D+00*(rhoCp(k-1)*dzi(k-1)+rhoCp(k) &
                                                *dzi(k))/(Hfw*rhosnow)
                 Tl=Tmelts
   !..........reduce snow thickness and replace heat that went unused
                 snmass=snmass-dhs*rhosnow
                 TerMelt=TerMelt+dhs
                 if(snmass.lt.0.D+00) then
                   TerMelt=TerMelt+snmass/rhosnow
                   Tl=Tl-snmass*Hfw/(0.5D+00*(rhoCp(k-1)*dzi(k-1) &
                                                        +rhoCp(k)*dzi(k)))
                   snmass=0.D+00
                 endif
                 Tice(k)=Tl
               endif
             enddo
           endif 
         endif
   !
   !....do top level if no snow (NOTE: surface temp. taken as Tmelts)
   
         if(nslay.eq.0) then
           Tl=Tice(1)
           if(Tl.gt.Tmelts) then
             dhi=(Tl-Tmelts)*rhoCp(1)*dzi(1)*0.5D+00/(Hfi*rhoice)
             Tl=Tmelts
   !........first use heat to melt remaining snow, if any
             snmass=snmass-dhi*rhoice
             if(snmass.lt.0.D+00) then
               dhi=-snmass/rhoice
               snmass=0.D+00
             else
               dhi=0.D+00
             endif
   !........reduce ice mass accordingly and replace heat that went unused
             simass=simass-dhi*rhoice
             TopMelt=TopMelt+dhi
             Tice(1)=Tl
           endif
         endif     
   !
   !....now check internal points (NOTE: internal points melt at Tmelti)
         do k=max(2,nslay+1),nilay
           Tl=Tice(k)
           if(Tl.gt.Tmelti) then
             dhi=(Tl-Tmelti)*0.5D+00*(rhoCp(k-1)*dzi(k-1)+rhoCp(k) &
                                                *dzi(k))/(Hfi*rhoice)
             Tl=Tmelti
   !........reduce ice thickness accordingly and replace heat that went unused
             simass=simass-dhi*rhoice
             TerMelt=TerMelt+dhi
             Tice(k)=Tl
           endif
         enddo
   !
   !...Check if negative ice has been created. If so, put heat into 
   !...mixed layer (Note, negative ice is used to account for additional
   !...heat which may remain after all ice is melted)
   !
         if(simass.le.0.D+00) then
   !......add negative ice to mixed layer heat !NSnote, this is the dtemp comp from Winton? => should go into heat, is this done???
           Hmix=Hmix-min(0.D+00,simass)*Hfi
           simass=0.D+00
         endif
   !
   !    LEVEL1'end growthtb'
   
      return  

end subroutine growthtb

!-----------------------------------------------------------------------

subroutine sebudget(hum_method,longwave_radiation_method,fluxes_method,&
                     TTss,airt,rh,cloud,ice_hi,ice_hs,&
                     lat,u10,v10,precip,airp,evap)
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
      integer, intent(in)       :: longwave_radiation_method ! method for LW
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
      call longwave_radiation(longwave_radiation_method, &
                          lat,TTss,airt+kelvin,cloud,qb)     ! subroutine longwave_radiation(method,dlat,tw,ta,cloud,ql)
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

!-----------------------------------------------------------------------

!      LEVEL1'surfmelt'
  
!
!...Calculate flux divergence at surface. Note: Iceflux(1) is the downward
!...positive heat flux defined at the centre of the uppermost layer
!
   Fdiv=Iceflux(1)-fluxt
   !
   !...If no melt is implied, return; otherwise, calculate surface melt
   ! still apply constant drainrate for meltponds
   !
         if(Fdiv.ge.0.D+00 .or. fluxt.lt.0.D+00 .or. Ts.lt.Tmelts) then
             if (meltpond) then 
                meltmass=meltmass-drainrate*dti*Amelt
               if (meltmass .lt. 0.D+00) then
                 meltmass= 0.D+00
                 Amelt=0.D+00
                 Aice=1-Asnow
               end if 
             end if
           return
         else
   !
   !......use heat to melt available snow, then ice
           dmsw=-Fdiv*dti/Hfw
           fluxt=0.D+00
           snmass=snmass-Asnow*dmsw
           simass=simass-(Aice+Amelt)*dmsw*Hfw/Hfi 
   !nsnote again double snow check! rm one?
             if (meltpond .and. snmass .gt. 0.D+00 .and. snmass/rhosnow .gt. hsmin) then
             meltmass=meltmass+dmsw*rhosnow/rhowaterfresh-drainrate*dti*Amelt
                   
                   if (meltmass .lt. 0.D+00) then
                       meltmass=0.D+00
                       Amelt=0.D+00
                       Aice=1-Asnow
                   end if  
             else if (meltpond .and. snmass .lt. 0.D+00) then
             meltmass=meltmass+dmsw*rhoice/rhowaterfresh-drainrate*dti*Amelt
                   if (meltmass .lt. 0.D+00) then
                       meltmass=0.D+00
                       Amelt=0.D+00
                       Aice=1-Asnow
                   end if
             end if
   
           if(snmass.lt.0.D+00) then
   !
   !........if all of snow melted, use remaining heat to melt ice
             if (snow_dist) then
                hsmax=hsmin
             end if
             dmi=-snmass*Hfw/Hfi
             snmass=0.D+00
             simass=simass-dmi
             TopMelt=TopMelt+dmi/rhoice
             if(simass.lt.0.D+00) then
   !
   !..........if all ice is melted, restore remaining surface flux to 
   !..........be added to mixed layer heat storage
               TopMelt=TopMelt+simass/rhoice
               fluxt=-simass*Hfi/dti
               simass=0.D+00
             endif
           endif
         endif
   
   !
   !    LEVEL1'end surfmelt'
   
      return


end subroutine surfmelt

!-----------------------------------------------------------------------



subroutine nr_iterate(hum_method,longwave_radiation_method,fluxes_method,&
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
   integer, intent(in)       :: longwave_radiation_method ! method for LW
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
   integer                   ::    nrit,ksearch,kb_uvic,l     !jpnote renamed kb 
   real(rk),parameter        ::    toler=1.D-02 

!-----------------------------------------------------------------------


#if 0
!      LEVEL1'nr_iterate'
!
!
!...Save initial temperature profile
!
   do l=1,nilay+1
      Told(l)=Tice(l)
   enddo
!
!...First guess at surface temperature
!
   Ts=Tice(1)

!
!-----------------------------------------------------------------------------
!--- Start initial Newton-Raphson iteration - do a maximum of 5 iterations ---
!-----------------------------------------------------------------------------
!
   do nrit=1,5
!      
!......calculate surface energy budget terms
      call sebudget(hum_method,longwave_radiation_method,fluxes_method,&
                    Ts,airt,rh,cloud,ice_hi,ice_hs,&
                    lat,u10,v10,precip,airp,evap)
      call  albedo_ice_uvic(fluxt,I_0,PenSW,alb,ice_hs,ice_hi,Ts)
!
!
!......restore intial temperature profile as therm1d takes a forward time step
      do l=1,nilay+1
         Tice(l)=Told(l)
      enddo

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
      call sebudget(hum_method,longwave_radiation_method,fluxes_method,&
                    Tsp,airt,rh,cloud,ice_hi,ice_hs,&
                    lat,u10,v10,precip,airp,evap)
      call  albedo_ice_uvic(fluxt,I_0,PenSW,alb,ice_hs,ice_hi,Tsp)
!
      do l=1,nilay+1
         Tice(l)=Told(l)
      enddo
      call therm1d(nilay,Sice_bulk,ice_hi,ice_hs,dzi, &
         Cond,rhoCp,zi,Sint,Pari,Tice,I_0)
      fTsdT1=Tsp-Tice(1)
      Tsp=Ts-dTemp
      call sebudget(hum_method,longwave_radiation_method,fluxes_method,&
                    Tsp,airt,rh,cloud,ice_hi,ice_hs,&
                    lat,u10,v10,precip,airp,evap)
      call  albedo_ice_uvic(fluxt,I_0,PenSW,alb,ice_hs,ice_hi,Tsp)
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
!
!----------------------------------------------------------------------
!--- If first try at Newton-Raphson iteration fails, try method of  ---
!--- bisection which is less efficient, but more robust             ---
!----------------------------------------------------------------------
!
!...To start, find where Ts-Tice(1) changes sign
!
!....first, do initial guess again
   Ts1=Told(1)
   call sebudget(hum_method,longwave_radiation_method,fluxes_method,&
                 Ts1,airt,rh,cloud,ice_hi,ice_hs,&
                    lat,u10,v10,precip,airp,evap)
      call  albedo_ice_uvic(fluxt,I_0,PenSW,alb,ice_hs,ice_hi,Ts1)
!
   do l=1,nilay+1
      Tice(l)=Told(l)
   enddo
   call therm1d(nilay,Sice_bulk,ice_hi,ice_hs,dzi, &
         Cond,rhoCp,zi,Sint,Pari,Tice,I_0)
   fTs1=Ts1-Tice(1)
!
!....increase and decrease Ts in 1 degree intervals until zero-crossing
!....is bracketed
   dTs=1.
   do ksearch=1,10
      Tsu=Ts1+dTs*float(ksearch-1)
      call sebudget(hum_method,longwave_radiation_method,fluxes_method,&
                    Tsu,airt,rh,cloud,ice_hi,ice_hs,&
                    lat,u10,v10,precip,airp,evap)
      call  albedo_ice_uvic(fluxt,I_0,PenSW,alb,ice_hs,ice_hi,Tsu)
!
      do l=1,nilay+1
         Tice(l)=Told(l)
      enddo
      call therm1d(nilay,Sice_bulk,ice_hi,ice_hs,dzi, &
           Cond,rhoCp,zi,Sint,Pari,Tice,I_0)
      fTsu=Tsu-Tice(1)
!......check if crossing is between Ts1 and Tsu>Ts1
      if(fTsu.gt.0.D+00.and.fTs1.le.0.D+00 &
           .or.                            &
        fTsu.lt.0.D+00.and.fTs1.ge.0.D+00) then
            Ts2=Tsu
            go to 887
      endif
!
      Tsl=Ts1-dTs*float(ksearch-1)
      call sebudget(hum_method,longwave_radiation_method,fluxes_method,&
                    Tsl,airt,rh,cloud,ice_hi,ice_hs,&
                    lat,u10,v10,precip,airp,evap)
      call  albedo_ice_uvic(fluxt,I_0,PenSW,alb,ice_hs,ice_hi,Tsl)
!
      do l=1,nilay+1
         Tice(l)=Told(l)
      enddo
      call therm1d(nilay,Sice_bulk,ice_hi,ice_hs,dzi, &
         Cond,rhoCp,zi,Sint,Pari,Tice,I_0)
      fTsl=Tsl-Tice(1)
!......check if crossing is between Ts1 and Tsl<Ts1
      if(fTsl.gt.0.D+00.and.fTs1.le.0.D+00 &
                       .or.                &
        fTsl.lt.0.D+00.and.fTs1.ge.0.D+00) then
            Ts2=Tsl
            go to 887
      endif
!
!.....if crossing not found, go back expand search
   enddo
!
!...If no zero-crossing is found, have to resort to a forward time
!...step with surface temperature from last time. First print message.
!  
   print*,'***Search for zero-crossing failed in nr_iterate!'
   go to 900
!
!...If zero-crossing found, come here to start method of Bisection
!
 887 continue
!
   do kb_uvic=1,20
!......one end of interval

      call sebudget(hum_method,longwave_radiation_method,fluxes_method,&
                    Ts1,airt,rh,cloud,ice_hi,ice_hs,&
                    lat,u10,v10,precip,airp,evap)
      call  albedo_ice_uvic(fluxt,I_0,PenSW,alb,ice_hs,ice_hi,Ts1)
!
      do l=1,nilay+1
         Tice(l)=Told(l)
      enddo
      call therm1d(nilay,Sice_bulk,ice_hi,ice_hs,dzi, &
        Cond,rhoCp,zi,Sint,Pari,Tice,I_0)
      fTs1=Ts1-Tice(1)
!......middle of interval
      Tsm=(Ts1+Ts2)/2.D+00
      call sebudget(hum_method,longwave_radiation_method,fluxes_method,&
                    Tsm,airt,rh,cloud,ice_hi,ice_hs,&
                    lat,u10,v10,precip,airp,evap)
      call  albedo_ice_uvic(fluxt,I_0,PenSW,alb,ice_hs,ice_hi,Tsm)
!
      do l=1,nilay+1
         Tice(l)=Told(l)
      enddo
      call therm1d(nilay,Sice_bulk,ice_hi,ice_hs,dzi, &
        Cond,rhoCp,zi,Sint,Pari,Tice,I_0)
      fTsm=Tsm-Tice(1)
!......check if solution found
      if(abs(fTsm).lt.toler) then
        print*,'Met tolerance in Bisection Method'
        go to 901
      endif
!......if not, find which half interval and try again
      if(fTs1*fTsm.gt.0.D+00) then
         Ts1=Tsm
      else
         Ts2=Tsm
      endif
!......check if interval is smaller than tolerance
      if(abs(Ts1-Ts2).lt.toler) then
        print*,'Interval < tolerance in Bisection Method'
        go to 901
      endif
   enddo
!
!...If method of Bisection fails, print message and do a forward
!...time step
!
   print*,'***Bisection scheme failed: Ts1 =',Ts1,'Ts2 =',Ts2
   go to 900
!
!...If method of Bisection is successful, take its solution, Tsm
!...as the surface temperature and do a forward time step
!
 901 continue
   Ts=Tsm
   go to 902
!
!------------------------------------------------------------------------
!--- If all else fails, print message and do a forward time step      ---
!--- using last surface temperature. This will allow model to proceed ---
!--- but may leave a spurious spasm in the surface temperature.       ---
!------------------------------------------------------------------------
!
 900 continue
   print*,'***!!!! Doing a forward time step anyway !!!***'
   Ts=Told(1)
 902 continue
   call sebudget(hum_method,longwave_radiation_method,fluxes_method,&
                 Ts,airt,rh,cloud,ice_hi,ice_hs,&
                    lat,u10,v10,precip,airp,evap)
      call  albedo_ice_uvic(fluxt,I_0,PenSW,alb,ice_hs,ice_hi,Ts)
!
   do l=1,nilay+1
      Tice(l)=Told(l)
   enddo
   call therm1d(nilay,Sice_bulk,ice_hi,ice_hs,dzi, &
        Cond,rhoCp,zi,Sint,Pari,Tice,I_0)
!
 987 continue
!
!
!    LEVEL1'end nr_iterate'

return

#endif

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

!-----------------------------------------------------------------------
!BOC
      ssi=0
      !get old time step bulk salinity
      !
         So=Sice_bulk
      
      !...convert rate of mass to thickness
      !
         ratei=(simass-simasso)/dti
         V = ratei/rhoice          !ice thickness change in m
      
      !Ice salinity conservation equation terms:
      !...brine entrapment during ice growth
      !...get the fractionation coefficient (Cox and Weeks, 1988)
      !    
           
              f=0.12D+00
         if(higb .gt. 0) then
           if(V .gt. 3.6e10-7) then
              f=0.26D+00/(0.26D+00+0.74D+00*exp(-724300D+00*V)) !V in m
           end if
         
           if(V .lt. 3.6e-7 .and. V .gt. 2e-8) then
              f=0.8925D+00+0.0568D+00*log(100.D+00*V)
           end if
      
      
           x1=(higb)*(f*SSS-So)/(simass/rhoice)     !NSnote fixed units 
         else
           x1=0
         end if
      
      !...brine entrapment during snow ice formation
      !...calculate if snow ice has been created
      !
         if(higs .gt. 0.D+00 .and. simass .gt. 0) then
           Ssi=(rhoice-rhosnow)*SSS/rhoice
           x2=(max(0.D+00,higs))*((Ssi-So)/(simass/rhoice)) 
         else
           x2=0
         end if
      
      !...gravity drainage term
         if(Tice(nilay+1) .gt. Tice(2)) then  !only if dT/dz < 0
            x3=-((So-Sg)/Tg)*dti       !decreases Sice_bulk
         else
            x3=0
         endif
      
      !...flushing term
         if(hims .gt. 0.0D+00) then
           x4=-((So-Sfl)/Tfl)*dti      !decreases Sice_bulk
         else
           x4=0
         end if
      
      !Get time-step ice bulk salinity!
      !
         Sice_bulk=So+x1+x2+x3+x4 !NSnote unit inconsistency x1,x2 - fixed
      !Get time-step salinity profile, which depends on the bulk salinity
      !
      !...get the ice salinity linear profile to be used in the cases of
      !...constrained from 0 ppt at the surface to the bulk salinity (Sice_bulk)
      !
         if(simass .gt. 0) then
           do k=1,nilay+1
             Szero(k) = Sice_bulk*zi(k)/(simass/rhoice)
           end do
         end if
      
         
         
      !...for a low bulk salinity (Sice_bulk < S1 = 3.5 ppt) >> linear profile
              if(Sice_bulk .lt. S1) alpha=1.D+00
      !...for a high bulk salinity (Sice_bulk > S2 = 4.5 ppt) >> const profile
              if(Sice_bulk .gt. S2) alpha=0.D+00
      !...for intermediate bulk salinity >> intermediate values between previous cases
              if(Sice_bulk .gt. S1 .and. Sice_bulk .lt. S2) then
                 alpha=(Sice_bulk-S2)/(S1-S2)
              end if
      
         do k=1,nilay+1
            Sice(k) = alpha*Szero(k)+(1.D+00-alpha)*Sice_bulk
         end do
      
      
      ! Ice-ocean salt flux (Tartinville et al. 2001)
      ! Terms representing salt rejection of the sea ice and
      ! salt input to the ocean.
      !
      ! Freshwater flux
      
      ! Flux due to snow melt
         if (meltpond.and.Amelt.ne.0.0) then
          rates=meltmass-meltmasso
         else  
          rates=snmass-snmasso
          rates=rates*(rhosnow/rhof)
         endif
      ! to be summed to ppt,evap and river inflow in GOTM
        Ff =  (-1.D+00)* min(0.D+00,rates)/dti    ![kg m-2 s-1] !sign positive if any freshwater goes into the ocean
        Ff = Ff/rhow   ![m s-1]
        rates=0
      
      ! Flux due to brine drainage
      !
         Fb=simass*(x3+x4)/dti
         
      ! Basal accretion
      !
        xx2=-rhoice*(SSS-Sice(nilay+1))*max(0.D+00,higb)/dti !Vancoppenole multiplies it by the ice fractionation at the mth categorie,
      !check how it would work properly using Hibler scheme (check 'g' in the black book)
      
      ! Snow ice formation
      !
        xx3=-rhoice*(SSS-Ssi)*(max(0.D+00,higs))/dti ![ppt kg m-2 s-1]
      
      ! Melt of saline ice
      !
        xx4=rhoice*(SSS-Sice_bulk)*(hims+himb+himi)/dti
      
        Feq=xx2+xx3+xx4 
      
        Fs = (Fb + Feq)/rhow         ! [ppt m s-1] - sign with respect to ice ! neg => S transfer from ice to water
                                     ! pos => S transfer from water to ice ! need to invert sign in S equation  
      
      
         return      

end subroutine saltice_prof_simple

!-----------------------------------------------------------------------

! call  albedo_ice_uvic(fluxt,I_0,PenSW,alb,ice_hs,ice_hi,Ts)
subroutine albedo_ice_uvic(fluxt,I_0,PenSW,alb,ice_hs,ice_hi,TTss) !??? renamed from albedo_ice to albedo_ice_uvic to avoid conflict with winton variable 
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
!snow_dist
!
!EOP
!-----------------------------------------------------------------------------
!BOP
!-----------------------------------------------------------------------------
!  LEVEL1 'albedo_ice'

!...snow_dist
      if (ice_hs .ge. hsmax.and.ice_hs.ne.0.D+00) then
      !nsnote isn't hsmax always gt 0.0 ?
               Asnow=1.D+00
               Aice=0.D+00
               Amelt=0.D+00
               meltmass=0.D+00
            else if (ice_hs .lt. hsmin .and. ice_hi .ge. hsmin) then
                  if (meltpond .and. meltmass .gt. 0.D+00) then
                     Amelt=Ameltmax
                     Aice=1.D+00-Amelt
                  else
                     Amelt=0.D+00
                     Aice=1.D+00
                  end if
               Asnow=0.D+00
            else if (ice_hi .lt. hsmin) then
               Aice=0.D+00
               Asnow=0.D+00
               Amelt=1.D+00
               meltmass=0.D+00
            else  
               if(snow_dist .and. distr_type .eq. 0) then
                     Asnow=exp(-1.D+00*(inverfc(ice_hs/hsmax)**2))
                  else if(snow_dist .and. distr_type .eq. 1) then
                     Asnow=exp(W(-2.D+00*ice_hs/hsmax*exp(-2.D+00))+2.D+00) &
                           *(-W(-2.D+00*ice_hs/hsmax*exp(-2.D+00))-2.D+00+1.D+00)
                  else!Nfix
                     Asnow=1.D+00 !Nfix
                  end if
                  if (meltpond .and. meltmass .gt. 0.D+00) then
                  if(snow_dist) then
                     Amelt=min(Ameltmax,1.D+00-Asnow)
                     Aice=1.D+00-Asnow-Amelt
                  else
                     Amelt=Ameltmax !uniform snow with meltponds
                     Asnow=1.0D+00-Amelt !uniform snow
                     Aice=0.0D+00 !Nfix
                  end if
                  else
                     Aice=1.D+00-Asnow
                     Amelt=0.D+00
                  end if
            endif
      !H!begin
      !  ice and snow albedo
         if (albice_method.eq.1) then
            albice=albice_f
            albsnow=albsnow_f
         else if (albice_method.eq.2) then
            albice=(0.44D+00*ice_hi**0.28D+00)+0.08D+00
            albice=max(albice,albmelt)
            albsnow=albsnow_f
         else if (albice_method.eq.3) then
            if (airtk.lt.273.05) then
            albice=max(albmelt,(0.44D+00*ice_hi**0.28D+00)+0.08D+00)
            else
            albice=min(albice_m,(0.075D+00*ice_hi**2.D+00)+albmelt)
            end if
            if (airtk.lt.273.15) then
            albsnow=albsnow_f
            else
            albsnow=albsnow_m
            end if
         else if (albice_method.eq.4) then
            albice=((albice_m+albice_f)+(albice_m-albice_f)*tanh(TTss-273.15D+00))/2.D+00
            albsnow=((albsnow_m+albsnow_f)+(albsnow_m-albsnow_f)*tanh(TTss-273.15D+00))/2.D+00
         end if
      !H!end
      !...Calculate penetrating SW flux
      !
      !
      
      
            PenSW=Asnow*I_0*(1.D+00-albsnow)*transs(TTss)+ &
                  Aice*I_0*(1.D+00-albice)*transi(TTss)+ &
                  Amelt*I_0*(1.D+00-albmelt)*transm   
      
      
      !   Calculate short wave flux absorbed at the surface
      
      
            I_0_tmp=Asnow*I_0*(1.D+00-albsnow)*(1.D+00-transs(TTss))+ &
                  Aice*I_0*(1.D+00-albice)*(1.D+00-transi(TTss))+ &
                  Amelt*I_0*(1.D+00-albmelt)*(1.D+00-transm)
      
      
      !....Add short wave flux absorbed at surface to total flux
               fluxt=fluxt+I_0_tmp
      
               alb=Asnow*albsnow+Aice*albice+Amelt*albmelt
      !
      !...end snow_dist
      !
      
         return

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
   !LEVEL1 'clean_ice_uvic'  !commented for now jpnote

   !LEVEL3 'closing ice.nc file...'  !commented for now  jpnote
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

!BOF
! !IROUTINE: snow extinction
!
! !INTERFACE:

real function swkappas(Tin)
!DESCRIPTION:

!USES:
             IMPLICIT NONE
!INPUT PARAMETERS:
             real(rk), intent(in) :: Tin
! !REVISION HISTORY:
!  Original author(s): 
!
!  See log for the icemodel module
!            
!-----------------------------------------------------------------------
!BOC
             swkappas=((swkappasm+swkappasf)+(swkappasm-swkappasf)* &
                     tanh(Tin-273.15D+00))/2.D+00

!test
!             swkappas=1.5               
 
end function swkappas
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
