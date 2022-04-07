!!   INPUT:    p      Pressure, mbar
!!             t      Temperature, deg C
!!             wetvar Water vapor variable, defined below
!!             h2otyp H2OTYP = 'R' --> WETVAR is relative humidity (0-1)
!!                    H2OTYP = 'D' --> WETVAR is dew point temperature (deg C)
!!             phi    Geocentric latitude (radians)
!!             h      Elevation above ellipsoid, km
!!             elev   Elevation angle, radians

real*8 function cfa(pr, temp, wvv, rdtyp, phi, h, elev)
!
implicit none
real*8 pr, temp, wvv, phi, h, elev, wp
real*8 a0, a1, a2, a3,a4, a5, b0, b1, b2, b3,b4, b5, c0
real*8 a, b, c, sinel, tanel
!
real*8 wpress
!
character*1 rdtyp, upper_string
!
a0=0.001185d+00
a1=0.6071d-04
a2=-0.1471d-03
a3=0.3072d-02
a4=0.1965d-01
a5=-0.5645d-02
b0=0.001144d+00
b1=0.1164d-04
b2=0.2795d-03
b3=0.3109d-02
b4=0.3038d-01
b5=-0.1217d-01
c0=-0.009d+00
!
if (rdtyp .eq. upper_string('R')) then
  wp = wpress(wvv, temp)
else
  wp = wpress(1.0d0, wvv)
endif
!
if (elev .eq. 90.d0) then
  cfa = 1.d0
  goto 100
endif

sinel = sin(elev)
if (dabs(sinel - 1.d0) .lt. 1.d-10) then
  cfa = 1.d0
  goto 100
endif
tanel = sinel/cos(elev)
!
a = a0*(1.0d+00 + a1*(pr - 1000.0) + a2*wp + a3*(temp - 20.0) + a5*(-h))
b = b0*(1.0d+00 + b1*(pr - 1000.0) + b2*wp + b3*(temp - 20.0) + b5*(-h))
c = c0
!
cfa = 1.0d0/(sinel + a/(tanel + b/(sinel + c)))
!
100 return
end
