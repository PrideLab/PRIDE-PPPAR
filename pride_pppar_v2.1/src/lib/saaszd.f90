!!   INPUT  : pr    (pressure(mbar))
!!            temp  (temperature)
!!            wvv   (Water vapor variable)
!!            rdtyp (R:relative humdity, D:dew point temperature)
!!            phi   (Geocentric latitude)
!!            elev  (Elevation above the geoid)

real*8 function saaszd(pr, temp, wvv, rdtyp, phi, elev)
implicit none
real*8 pr, temp, wvv, phi, elev, wp
character*1 rdtyp
!! function
real*8 wpress,ffun
character*1 upper_string
!
if (rdtyp .eq. upper_string('R')) then
  wp = wpress(wvv, temp)
else
  wp = wpress(1.d0, wvv)
endif
!! zenith delay
saaszd = 0.2277d-2*((0.1255d4/(temp + 273.15d0) + 0.5d-1)*wp + pr)/ffun(phi, elev)
return
end
