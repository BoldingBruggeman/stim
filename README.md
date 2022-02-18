# Simple Thermodynamic Ice Models (STIM)

STIM provides a common interface of ice models to hydrodynamic models.

The aim is to separate the actual ice algoorithms from the hydrodynamic model such that only a well-defined small interface is needed in order to test and use a variaty of different ice models.

As a proof of concept two different models are included in this initial release. Also included is a driver for [GOTM](https://www.gotm.net).

### Code documentation (alpha - and not always up to date)

Source code [documentation](https://bolding-bruggeman.com/portfolio/eat/ford/)

### API change

The API for _do\_stim_ in the GOTM-driver has changed between the legacy and master branch:

Legacy
```
subroutine do_stim(dz,dt,Tw,S,Ta,precip,Qsw,Qfluxes)
```

Master
```
  SUBROUTINE do_stim(dz,dt,ustar,Tw,S,Ta,precip,Qsw,Qfluxes)
 
    !! Arguments
   real(rk), intent(inout)    :: dz
      !! layer thickness [m]
   real(rk), intent(inout)    :: dt
      !! time step [s]
   real(rk), intent(inout)    :: ustar
      !! surface friction velocity [m/s]
   real(rk), intent(inout) :: Tw
      !! water temperature [C]
   real(rk), intent(inout)    :: Ta
      !! air temperature [C]
   real(rk), intent(inout)    :: S
      !! salinity [g/kg]
   real(rk), intent(inout)    :: precip
      !! precipitation [mm?]
   real(rk), intent(inout)    :: Qsw
      !! short wave radiation [W/m^2]
      
   interface
      SUBROUTINE Qfluxes(T,qh,qe,qb)
         use, intrinsic :: iso_fortran_env
         integer, parameter :: rk=real64

         real(rk), intent(in) :: T
           !! temperature [C]
         real(rk), intent(out) :: qh
           !! latent heat [W/m^2]
         real(rk), intent(out) :: qe
           !! sensible heat [W/m^2]
         real(rk), intent(out) :: qb
           !! net longwave radiation [W/m^2]
      END SUBROUTINE
   end interface
```
