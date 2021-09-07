subroutine vmf1_ht(ah, aw, dmjd, dlat, ht, zd, vmf1h, vmf1w)

!     !!! This is the version with height correction !!!
!     !!! It has to be used with the grid !!!
!
!     This subroutine determines the VMF1 (Vienna Mapping Functions 1)
!     Reference: Boehm, J., B. Werl, H. Schuh (2006),
!     Troposphere mapping functions for GPS and very long baseline interferometry
!     from European Centre for Medium-Range Weather Forecasts operational analysis data,
!     J. Geoph. Res., Vol. 111, B02406, doi:10.1029/2005JB003629.
!
!     input data
!     ----------
!     ah:   hydrostatic coefficient a (www.hg.tuwien.ac.at/~ecmwf1)
!     aw:   wet coefficient a         (www.hg.tuwien.ac.at/~ecmwf1)
!     dmjd: modified julian date
!     dlat: latitude in radians
!     ht:   ellipsoidal height in meter
!     zd:   zenith distance in radians
!
!     output data
!     -----------
!     vmf1h: hydrostatic mapping function
!     vmf1w: wet mapping function
!
!     Johannes Boehm, 2005 October 2
!

  implicit double precision(a - h, o - z)

  pi = 3.14159265359d0

!     reference day is 28 January
!     this is taken from Niell (1996) to be consistent
  doy = dmjd - 44239.d0 + 1 - 28

  bh = 0.0029d0
  c0h = 0.062d0
  if (dlat .lt. 0.d0) then   ! southern hemisphere
    phh = pi
    c11h = 0.007d0
    c10h = 0.002d0
  else                     ! northern hemisphere
    phh = 0.d0
    c11h = 0.005d0
    c10h = 0.001d0
  end if
  ch = c0h + ((dcos(doy/365.25d0*2.d0*pi + phh) + 1.d0)*c11h/2.d0 + c10h)*(1.d0 - dcos(dlat))

  sine = dsin(pi/2.d0 - zd)
  beta = bh/(sine + ch)
  gamma = ah/(sine + beta)
  topcon = (1.d0 + ah/(1.d0 + bh/(1.d0 + ch)))
  vmf1h = topcon/(sine + gamma)

!  height correction [Niell, 1996]
  a_ht = 2.53e-5
  b_ht = 5.49e-3
  c_ht = 1.14e-3
  hs_km = ht/1000.d0
  beta = b_ht/(sine + c_ht)
  gamma = a_ht/(sine + beta)
  topcon = (1.d0 + a_ht/(1.d0 + b_ht/(1.d0 + c_ht)))
  ht_corr_coef = 1.d0/sine - topcon/(sine + gamma)
  ht_corr = ht_corr_coef*hs_km
  vmf1h = vmf1h + ht_corr

  bw = 0.00146d0
  cw = 0.04391d0
  beta = bw/(sine + cw)
  gamma = aw/(sine + beta)
  topcon = (1.d0 + aw/(1.d0 + bw/(1.d0 + cw)))
  vmf1w = topcon/(sine + gamma)
  return
end subroutine
