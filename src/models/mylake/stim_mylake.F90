!> The MyLAKE ice model
!>
!> authors: Karsten Bolding (after ?? and ??)

   MODULE stim_mylake

   use stim_variables, only: rk
   use stim_variables, only: albedo_ice, attenuation_ice
   use stim_variables, only: Cw,K_ice,L_ice,rho_ice
   use stim_variables, only: Tf,Tice_surface
   use stim_variables, only: Hice,Hfrazil,dHis,dHib
   use stim_variables, only: transmissivity
   use stim_variables, only: surface_ice_energy, bottom_ice_energy
   use stim_variables, only: ocean_ice_flux

   IMPLICIT NONE

   private

   public init_stim_mylake, do_stim_mylake, clean_stim_mylake

   real(rk), pointer :: Tice
   real(rk), pointer :: ice_energy
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
   SUBROUTINE init_stim_mylake()
!-----------------------------------------------------------------------
   Tice => Tice_surface
   ice_energy => bottom_ice_energy

!KB   albedo => albedo_ice
!KB   attenuation => attenuation_ice

   ! https://github.com/biogeochemistry/MyLake_public/blob/master/v12/v12_1/VAN_para_v12_1b.xls
   ! Phys_par(13) = lambda_i=5
   attenuation_ice = 0.7_rk
   attenuation_ice = 5.0_rk
   transmissivity = exp(-Hice*attenuation_ice)
   END SUBROUTINE init_stim_mylake

!-----------------------------------------------------------------------

   SUBROUTINE do_stim_mylake(ice_cover,dz,dt,Tw,S,Ta,precip,Qsw,Qfluxes)

   real(rk), intent(in) :: dz,dt,Ta,S,precip,Qsw
   integer, intent(inout)  :: ice_cover
   real(rk), intent(inout) :: Tw

   interface
      SUBROUTINE Qfluxes(T,qh,qe,qb)
         integer, parameter :: rk = kind(1.d0)
         real(rk), intent(in) :: T
         real(rk), intent(out) :: qh,qe,qb
      END SUBROUTINE
   end interface

   real(rk) :: max_frazil=0.03_rk
   real(rk) :: alpha
   real(rk) :: qh,qe,qb,Qflux
!-----------------------------------------------------------------------
   Tf = -0.0575*S

   ice_energy = (Tw-Tf)*dz*Cw + ocean_ice_flux*dt
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
      end if

   else ! Ice-cover - frazil or solid

      Tw=Tf

      ! surface of ice
      if (Ta .lt. Tf) then ! top ice growth when air temperature is below freezing pt
         dHis = Hice
         alpha = 1._rk/(10.*Hice)
         Tice = (alpha*Tf+Ta)/(1._rk+alpha) ! mylake - but why?
         Hice = sqrt(Hice**2+2._rk*K_ice/(rho_ice*L_ice)*dt*(Tf-Tice)) !Stefan's law
         dHis = Hice - dHis
      else
         Tice = Tf        ! top ice melting due to solar radiation and heat fluxes
         Tice = 0._rk     ! top ice melting due to solar radiation and heat fluxes
         call Qfluxes(Tice,qh,qe,qb)
         Qflux = qh+qe+qb
         dHis = -dt*(Qsw+Qflux)/(rho_ice*L_ice)
         Hice = Hice+min(0._rk,dHis)
      end if
      Hice = Hice+dt*precip

      ! bottom of ice - melting or freezing depending in flux direction
      Hice = Hice + dHib

      if (Hice .le. 0.) then  ! excess of melting energy returned to water temp
         Tw = -Hice*rho_ice*L_ice/(dz*Cw) + Tf
         ice_cover = 0       ! no ice
         Hice = 0._rk
         attenuation_ice = 0._rk
         transmissivity = 1._rk
      else
         albedo_ice = 0.3
         transmissivity = exp(-Hice*attenuation_ice)
      end if
   end if
   END SUBROUTINE do_stim_mylake

!-----------------------------------------------------------------------
! !INTERFACE:
   SUBROUTINE clean_stim_mylake()
!-----------------------------------------------------------------------
!BOC
   !LEVEL2 'clean_stim_mylake'
   END SUBROUTINE clean_stim_mylake

!-----------------------------------------------------------------------

   END MODULE stim_mylake

!-----------------------------------------------------------------------
! Copyright by the STIM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
