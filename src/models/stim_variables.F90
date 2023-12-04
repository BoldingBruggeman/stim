!> Variables common to all ice models.
!>
!> Not all variables are used by all models.
!>
!> authors: Karsten Bolding

   MODULE stim_variables

   use, intrinsic :: iso_fortran_env

   IMPLICIT NONE

   integer, parameter, public :: rk=real64
   real(rk), public :: ice_cord = 1._rk
   real(rk), public :: Cw = 4.18e+6_rk
     !! volumetric heat capacity of water (J K-1 m-3)
   real(rk), public :: L_ice = 333500._rk
     !! latent heat of freezing (J kg-1)
   real(rk), public :: K_ice = 2.1_rk
     !! ice heat conduction coefficient (W m-1 K-1)
   real(rk), public :: lambda_ice = 5._rk !
   real(rk), public :: rho_ice = 910._rk
     !! kg m-3

   real(rk), public :: init_ice_energy = 0._rk    !
   real(rk), target, public :: surface_ice_energy = 0._rk    !
   real(rk), target, public :: bottom_ice_energy = 0._rk    !

   real(rk), target, public :: Tf = 0._rk
     !! freezing temperature
   real(rk), target, public :: Tice_surface = 0._rk
   real(rk), allocatable, target, public :: Tice(:)
   real(rk), target, public :: Hfrazil = 0._rk
     !! Total frazil ice thickness
   real(rk), target, public :: Hsnow = 0._rk
     !! Total snow thickness
   real(rk), target, public :: Hice = 0._rk
     !! Total ice thickness
   real(rk), target, public :: dHis = 0._rk
     !! surface ice growth
   real(rk), target, public :: dHib = 0._rk
     !! bottom ice growth
   real(rk), target, public :: albedo_ice = 0._rk
   real(rk), target, public :: attenuation_ice = 0._rk
   real(rk), target, public :: melt_rate = 0._rk
      !! ice melt rate [m/y] !KB why?????
   real(rk), target, public :: T_melt = 0._rk
     !! melt layer temperature
   real(rk), target, public :: S_melt = 0._rk
     !! melt layer salinity
   real(rk), target, public :: ocean_ice_flux = 0._rk
     !! downward water flux [kg/m2]
   real(rk), target, public :: ocean_ice_heat_flux = 0._rk
     !! upward heat flux [W/m2]
   real(rk), target, public :: ocean_ice_salt_flux = 0._rk
     !! upward salinity flux [(m/s)*(g/kg)]
   real(rk), target, public :: transmissivity = 1._rk
   real(rk), target, public :: z0i = 0.01_rk

   ! Lebedev
   real(rk), target, public :: fdd
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
   SUBROUTINE init_stim_variables(ice_model)

   integer, intent(in)     :: ice_model
!-------------------------------------------------------------------------

   END SUBROUTINE init_stim_variables

!-----------------------------------------------------------------------
   SUBROUTINE clean_stim_variables()
!-----------------------------------------------------------------------

   END SUBROUTINE clean_stim_variables

!-----------------------------------------------------------------------

   END MODULE stim_variables

!-----------------------------------------------------------------------
! Copyright (C) 2019 - Karsten Bolding (BB)                            !
!-----------------------------------------------------------------------
