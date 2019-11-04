!-----------------------------------------------------------------------
!BOP
!
! !MODULE: 'lebedev' ice module
!
! !INTERFACE:
   module stim_lebedev
!
! !DESCRIPTION:
!  https://agupubs.onlinelibrary.wiley.com/doi/10.1002/2016JC012199
!
! !USES:
   use stim_variables, only: transmissivity, albedo_ice
   use stim_variables, only: rk, Hice, dHis, dHib, Tf, fdd
   use stim_variables, only: rho_ice, L_ice, ocean_ice_flux
   IMPLICIT NONE
!  Default all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public init_stim_lebedev, do_stim_lebedev, clean_stim_lebedev
!
! !PUBLIC DATA MEMBERS:
!
! !PRIVATE DATA MEMBERS:
!  Constants from Adolf Stips
   real(rk)           :: lebedev_fac=1.33_rk
   real(rk)           :: damp_leb_swr = -1.6_rk
   real(rk)           :: damp_leb_wind = -1.6_rk
   real(rk)           :: damp_leb_shf = -1.6_rk
   real(rk)           :: freeze_fac = -0.0575_rk
   real(rk)           :: lebedev_albedo = 0.545
!  REALTYPE :: sens_ice_water = 46.9d0 ! sensible heta transfer
! basal heat flux fixed at 20 W/m^2 according to Holland 1998
!   REALTYPE :: sens_ice_water = 20.d0 ! sensible heta transfer
!   REALTYPE :: dt_ice=3600.
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialisation of the ice variables
!
! !INTERFACE:
   subroutine init_stim_lebedev(ice_cover)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
   integer, intent(inout)  :: ice_cover
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  See log for the ice module
!
! !LOCAL VARIABLES:
   integer                   :: rc
!EOP
!-----------------------------------------------------------------------
!BOC
   if (Hice .gt. 0._rk) then
      fdd = (100._rk*Hice/lebedev_fac)**(1._rk/0.58_rk)
      ice_cover = 2
   else
      fdd = 0._rk
      ice_cover = 0
   end if
   return
   end subroutine init_stim_lebedev
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: do the 'lebedev' ice calculations
!
! !INTERFACE:
   subroutine do_stim_lebedev(ice_cover,dt,Tw,S,Ta,precip)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   real(rk), intent(in)    :: dt,Ta,S,precip
!
! !INPUT/OUTPUT PARAMETERS:
   integer, intent(inout)  :: ice_cover
   real(rk), intent(inout) :: Tw
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  See log for the ice module
!
! !LOCAL VARIABLES:
   real(rk) :: x
!EOP
!-----------------------------------------------------------------------
!BOC
   ! calculate ice thickness according to Lebedev 1938 
   Tf = freeze_fac*S
   if (ice_cover .eq. 2) then
      ! calculate the cumulative freezing days
      fdd = fdd+(Tf-Ta)*dt/86400._rk
   else ! ice_hi was 0._rk
      if (Ta .lt. Tf .and. Tw .lt. Tf) then
         fdd = fdd+(Tf-Ta)*dt/86400._rk
      end if 
   end if

   ! have at least 1 cm of ice,.. - the melting is not done correct - KB
   if (fdd .gt. 1._rk) then
      x = Hice
      Hice = 0.01_rk*lebedev_fac*fdd**0.58_rk
      dHis = Hice - x
#if 0
      if (dHis .lt. 0._rk) then
         dHib = (-dt*sensible_ice_water)/(rho_ice*L_ice)
      else
         dHib = 0._rk
      end if
#endif
      Tw = Tf
      ice_cover = 2
      albedo_ice = lebedev_albedo
      transmissivity = exp(Hice/damp_leb_swr)
!      transmissivity = exp(Hice*damp_leb_swr)
      if (Hice .lt. 0._rk) then
         albedo_ice = 0._rk
         transmissivity = 1._rk
         Hice = 0._rk
         ice_cover = 0
      end if
   else
      fdd = 0._rk
      Hice = 0._rk
      dHis = 0._rk
      ice_cover = 0
   end if

#if 0
   ! if we have ice change the heatfluxes !
   where (ice_hi .gt. 0._rk)
      ! damp the shortwave radiation
      swr = swr * exp(damp_leb_swr*ice_hi ) 
      albedo = albedo_ice
      ! damp the wind stress
      tausx = tausx * exp(damp_leb_wind*ice_hi ) 
      tausy = tausy * exp(damp_leb_wind*ice_hi )
      ! damp the heatflux ­ this is of course tricky as latent heatflux and
      ! longwave radiation behave different from sensible heatflux
      ! flux routine must be rewritten to treat them all separate
      ! heat = heat * exp(damp_leb_shf * ice_hi)
      ! estimate sensible heat flux between ice and ocean water shf = sens_ice_water * ( freezing_fac * sss ­ sst )
   end where
#endif

   return
   end subroutine do_stim_lebedev
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Cleaning up the 'lebedev' ice variables
!
! !INTERFACE:
   subroutine clean_stim_lebedev()
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  See log for the ice module
!
!EOP
!-----------------------------------------------------------------------
!BOC
   !LEVEL2 'clean_stim_lebedev'

   !LEVEL2 'de-allocation ice memory ...'
   !LEVEL2 'done.'

   return
   end subroutine clean_stim_lebedev
!EOC

!-----------------------------------------------------------------------

   end module stim_lebedev

!-----------------------------------------------------------------------
! Copyright by the STIM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
