real*8 function asknewet(e,Tm,lambda)

!! This subroutine determines the zenith wet delay based on the
!! equation 22 by Aske and Nordius (1987) 
!!
!! c Reference:
!! Askne and Nordius, Estimation of tropospheric delay for microwaves from
!! surface weather data, Radio Science, Vol 22(3): 379-386, 1987.
!!
!! input parameters:
!!
!! e:      water vapor pressure in hPa 
!! Tm:     mean temperature in Kelvin
!! lambda: water vapor lapse rate (see definition in Askne and Nordius 1987)
!! 
!! output parameters:
!!
!! zwd:  zenith wet delay in meter 
!!
!! Example 1 :
!!
!! e =  10.9621 hPa
!! Tm = 273.8720
!! lambda = 2.8071
!!
!! output:
!! zwd = 0.1176 m
!!
!! Johannes Boehm, 3 August 2013
!! Johannes Boehm, 24 December 2014, converted to Fortran
!! ---

implicit double precision (a-h,o-z)

double precision k1,k2,k2p,k3,lambda

!% coefficients
k1  = 77.604d0                     !% K/hPa
k2 = 64.79d0                       !% K/hPa
k2p = k2 - k1*18.0152d0/28.9644d0  !% K/hPa
k3  = 377600.d0                    !% KK/hPa

!% mean gravity in m/s**2
gm = 9.80665d0
!% molar mass of dry air in kg/mol
dMtr = 28.965d-3
!% universal gas constant in J/K/mol
R = 8.3143d0

!% specific gas constant for dry consituents
Rd = R/dMtr 

asknewet = 1.0d-6*(k2p + k3/Tm)*Rd/(lambda + 1.d0)/gm*e
return
end
