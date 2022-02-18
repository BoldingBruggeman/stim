!> The basal_melt ice model
!>
!> authors: Hans Burchard and Karsten Bolding

   module stim_basal_melt

   use stim_variables, only: rk
   use stim_variables, only: Hice
   use stim_variables, only: z0 => z0i
   use stim_variables, only: rhoi => rho_ice
   use stim_variables, only: Li => L_ice
   use stim_variables, only: vm => melt_rate
   use stim_variables, only: Tb => T_melt
   use stim_variables, only: Sb => S_melt
   use stim_variables, only: fH => ocean_ice_heat_flux
   use stim_variables, only: fS => ocean_ice_salt_flux

   IMPLICIT NONE

   private

   public init_stim_basal_melt, do_stim_basal_melt, clean_stim_basal_melt

!-----------------------------------------------------------------------

   contains

   SUBROUTINE init_stim_basal_melt(ice_cover)
     !! Initialize the basal_melt model with an ice-cover

   IMPLICIT NONE

   integer, intent(inout)  :: ice_cover
   integer                   :: rc
!-----------------------------------------------------------------------

   END SUBROUTINE init_stim_basal_melt

!-----------------------------------------------------------------------

   SUBROUTINE do_stim_basal_melt(h,ustar,T,S)

   IMPLICIT NONE

   real(rk), intent(in) :: h      
     !! upper layer thickness [m]
   real(rk), intent(inout) :: ustar  
     !! friction velocity ice/water interface [m/s]
   real(rk), intent(in) :: T
     !! upper layer temperature [C]
   real(rk), intent(in) :: S      
     !! upper layer salinity [g/kg] 

#if 0
   ! param. for freezing in-situ temperature eq. 
   real(rk) :: l1 =  -0.0575    
   real(rk) :: l2 =   0.0901
   real(rk) :: l3 =  7.61e-4
#else
   ! param. for freezing potential temperature eq. 
   real(rk), parameter :: l1 = -5.6705121472405570E-002 
   real(rk), parameter :: l2 =  7.5436448744204881E-002
   real(rk), parameter :: l3 =  7.6828512903539831E-004
#endif

   real(rk) :: nu     = 1.95e-6 
     !! molecular viscosity [m2/s]
   real(rk) :: kappa  = 0.4     
     !! von Karman constant [-]
   real(rk) :: Bs     = 8.5     
     !! constant for tracer roughness [-]
   real(rk) :: Pr_turb= 0.7     
     !! turbulent Prandtl number [-]
   real(rk) :: Pr_temp= 13.8    
     !! molecular Prandtl for temp [-]
   real(rk) :: Pr_salt= 2432.   
     !! turbulent Prandtl for salt [-]
   real(rk) :: c      = 4180.   
     !! heat capacity of sea water [J/(kg K)]
   real(rk) :: ci     = 1995.   
     !! heat capacity of glacial ice [J/(kg K)]
   real(rk) :: Ti     = -20.    
     !! ice core temperature [deg C]
!KB   real(rk) :: Li     = 3.33e5  
!KB     !! latent heat of fusion [J/kg]
   real(rk) :: rho0   = 1030.   
     !! reference density of sea water [kg/m3]
!KB   real(rk) :: rhoi   =  920.   
!KB     !! reference density of glacial ice [kg/m3]

   real(rk) ::  log_zk_p_z0_z0t,log_zk_p_z0_z0s
   real(rk) ::  betaT,betaS
   real(rk) ::  A1T,A1S,LL
   real(rk) ::  s1,s2,s3,pp,qq

   real(rk) :: X,beta
!-----------------------------------------------------------------------
#if 0
   X=0.55*exp(0.5*kappa*Bs)*sqrt(z0*ustar/nu)

   beta=X*(Pr_temp**(2./3.)-0.2)-Pr_turb*Bs+9.5
   log_zk_p_z0_z0t=log((0.5*h+z0)/z0)+kappa/Pr_turb*beta

   beta=X*(Pr_salt**(2./3.)-0.2)-Pr_turb*Bs+9.5
   log_zk_p_z0_z0s=log((0.5*h+z0)/z0)+kappa/Pr_turb*beta
#else
   betaT=0.55*exp(0.5*kappa*Bs)*sqrt(z0*ustar/nu)
   betaT=betaT*(Pr_temp**(2./3.)-0.2)
   betaT=betaT-Pr_turb*Bs+9.5

   betaS=0.55*exp(0.5*kappa*Bs)*sqrt(z0*ustar/nu)
   betaS=betaS*(Pr_salt**(2./3.)-0.2)
   betaS=betaS-Pr_turb*Bs+9.5

   log_zk_p_z0_z0t=log((0.5*h+z0)/z0)+kappa/Pr_turb*betaT
   log_zk_p_z0_z0s=log((0.5*h+z0)/z0)+kappa/Pr_turb*betaS
#endif

   A1T=kappa*ustar/(Pr_turb*log_zk_p_z0_z0t)
   A1S=kappa*ustar/(Pr_turb*log_zk_p_z0_z0s)
!KB   LL=l2+l3*zb
   LL=l2+l3*Hice*rhoi/rho0

   s1=l1*(A1T-A1S*ci/c)
   s2=A1S*S*l1*ci/c-A1S/c*(ci*(LL-Ti)+Li)-A1T*(T-LL)
   s3=A1S*S/c*(ci*(LL-Ti)+Li)
   pp=s2/s1
   qq=s3/s1

   Sb = -0.5*pp+sqrt(0.25*pp**2-qq)     
     !! melt layer salinity
   Tb = l1*Sb+LL                        
     !! melt layer temperature

   vm = A1T*(T-Tb)/(ci/c*(Tb-Ti)+Li/c)  
     !! melt rate in m/s (>0 for melting)

   fH = (A1T*(T-Tb)-vm*Tb)*c*rho0       
     !! upward heat flux in W/m2 (>0 for melting)

   fS = A1S*(S-Sb)-vm*Sb                
     !! upward salinity flux (should be zero for glacial ice)

   if (ustar < 1e-10) then
      fH = 0.0
      fS = 0.0
      vm = 0.0
   end if
   END SUBROUTINE do_stim_basal_melt

!-----------------------------------------------------------------------
!
   SUBROUTINE clean_stim_basal_melt()
!
   IMPLICIT NONE
!-----------------------------------------------------------------------

   END SUBROUTINE clean_stim_basal_melt
!EOC

!-----------------------------------------------------------------------

   end module stim_basal_melt

!-----------------------------------------------------------------------
! Copyright by the STIM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
