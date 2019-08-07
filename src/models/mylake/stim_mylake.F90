!-----------------------------------------------------------------------
!BOP
!
! !MODULE: 'mylake' ice module
!
! !INTERFACE:
   module stim_mylake
!
! !DESCRIPTION:
!
! !USES:
   use stim_variables, only: rk
   use stim_variables, only: albedo_ice, attenuation_ice
   use stim_variables, only: Cw,K_ice,L_ice,rho_ice
   use stim_variables, only: Tf,Tice_surface
   use stim_variables, only: Hice,Hfrazil,dHis,dHib
   IMPLICIT NONE
!  Default all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public init_stim_mylake, do_stim_mylake, clean_stim_mylake
!
! !PUBLIC DATA MEMBERS:
!
! !PRIVATE DATA MEMBERS:
   real(rk), pointer :: Tice, albedo, attenuation
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
   subroutine init_stim_mylake()
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
#if 0
   integer, intent(in)                      :: namlst
   character(len=*), intent(in)             :: fn
#endif
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
   !LEVEL2 'init_stim_mylake'

!  Read namelist from file.
#if 0
   open(namlst,file=fn,status='old',action='read',err=80)
   !LEVEL2 'reading ice namelists..'
   read(namlst,nml=ice,err=81)
   close (namlst)
   !LEVEL2 'done.'
#endif

   Tice => Tice_surface
   albedo => albedo_ice
   attenuation => attenuation_ice

   !LEVEL2 'done.'

   return
#if 0
80 FATAL 'I could not open: ',trim(fn)
   stop 'init_ice'
81 FATAL 'I could not read "ice" namelist'
   stop 'init_ice'
#endif

   end subroutine init_stim_mylake
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: do the 'mylake' ice calculations
!
! !INTERFACE:
   subroutine do_stim_mylake(ice_cover,dz,dt,Tw,S,Ta,precip,Qsw,Qfluxes)
!
! !DESCRIPTION:
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
!  See log for the ice module
!
! !LOCAL VARIABLES:
   real(rk)                  :: max_frazil=0.03_rk
   real(rk)                  :: alpha
   real(rk)                  :: ice_energy
   real(rk)                  :: qh,qe,qb,Qflux
!EOP
!-----------------------------------------------------------------------
!BOC

   Tf = -0.0575*S

   ice_energy = (Tw-Tf)*dz*Cw
!   dHib = (-ice_energy+sensible_ice_water)/(rho_ice*L_ice)
   dHis = 0._rk
   dHib = (-ice_energy)/(rho_ice*L_ice)

   if (ice_cover .eq. 0) then ! No ice
      if (dHib .gt. 0._rk) then
         Hfrazil = Hfrazil + dHib + dt*precip
      end if
      if (Hfrazil .ge. max_frazil) then      !initial freezing when frazil ice is above threshold value
         ice_cover = 2
         Hice = Hfrazil
         Hfrazil = 0._rk
      end if
      if (Hfrazil .lt. 0.) then  ! excess of melting energy returned to water temp
         Tw = -Hfrazil*rho_ice*L_ice/(dz*Cw)+Tf
         Hfrazil = 0._rk
      endif

   else ! Ice-cover - frazil or solid

      Tw=Tf

!     surface of ice
      if (Ta .lt. Tf) then ! top ice growth when air temperature is below freezing pt
         dHis = Hice
         alpha = 1._rk/(10.*Hice)
         Tice = (alpha*Tf+Ta)/(1._rk+alpha) ! mylake - but why?
         Hice = sqrt(Hice**2+2._rk*K_ice/(rho_ice*L_ice)*dt*(Tf-Tice)) !Stefan's law
         dHis = Hice - dHis
      else
         Tice = Tf         ! top ice melting due to solar radiation and heat fluxes
         Tice = 0._rk     ! top ice melting due to solar radiation and heat fluxes
         call Qfluxes(Tice,qh,qe,qb)
         Qflux = qh+qe+qb
         dHis = -dt*(Qsw+Qflux)/(rho_ice*L_ice)
         Hice = Hice+min(0._rk,dHis)
      endif
      Hice = Hice+dt*precip

!     bottom of ice - melting or freezing depending in flux direction
      Hice = Hice + dHib

      if (Hice .le. 0.) then  ! excess of melting energy returned to water temp
          Tw = -Hice*rho_ice*L_ice/(dz*Cw) + Tf
          ice_cover = 0       ! no ice
          Hice = 0._rk
          attenuation = 0._rk
      else
         albedo = 0.3
      endif
   end if

   return
   end subroutine do_stim_mylake
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Cleaning up the 'mylake' ice variables
!
! !INTERFACE:
   subroutine clean_stim_mylake()
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
   !LEVEL2 'clean_stim_mylake'

   return
   end subroutine clean_stim_mylake
!EOC

!-----------------------------------------------------------------------

   end module stim_mylake

!-----------------------------------------------------------------------
! Copyright by the STIM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
