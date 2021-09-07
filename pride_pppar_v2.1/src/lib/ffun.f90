!!   INPUT  :
!!            PHI    Geocentric latitude, radians
!!            H      Elevation of site above geoid, km

real*8 function ffun(phi, valu)
  implicit none
  real*8 phi, valu
  ffun = 1.d0 - 0.266d-2*dcos(2.d0*phi) - 0.28d-3*valu
  return
end
