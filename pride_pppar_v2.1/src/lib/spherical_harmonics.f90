!
subroutine spherical_harmonics(dlat, dlon, V, W, undu)
!
!     This subroutine determines Global Spherical Harmonics up to degree and order 9
!
!     input data
!     ----------
!     dlat: latitude in radians
!     dlon: longitude in radians
!
!     output data
!     -----------
!     V,W : spherical harmonics up to 9x9
!     undu: Geoid undulation in m (from a 9x9 EGM based model)
!
!     Johannes Boehm, 2006 June 12
!     rev 2006 June 16: geoid undulation is accounted for
!     ref 2006 Aug. 14: recursions for Legendre polynomials (O. Montenbruck)
!
  implicit none

  REAL*8 dlat, dlon, undu
  REAL*8 V(10, 10), W(10, 10)

  REAL*8 a_geoid(55), b_geoid(55)

  INTEGER*4 i, n, m, nmax, mmax
  REAL*8 x, y, z

  data(a_geoid(i), i=1, 55)/ &
    -5.6195d-001, -6.0794d-002, -2.0125d-001, -6.4180d-002, -3.6997d-002, &
    +1.0098d+001, +1.6436d+001, +1.4065d+001, +1.9881d+000, +6.4414d-001, &
    -4.7482d+000, -3.2290d+000, +5.0652d-001, +3.8279d-001, -2.6646d-002, &
    +1.7224d+000, -2.7970d-001, +6.8177d-001, -9.6658d-002, -1.5113d-002, &
    +2.9206d-003, -3.4621d+000, -3.8198d-001, +3.2306d-002, +6.9915d-003, &
    -2.3068d-003, -1.3548d-003, +4.7324d-006, +2.3527d+000, +1.2985d+000, &
    +2.1232d-001, +2.2571d-002, -3.7855d-003, +2.9449d-005, -1.6265d-004, &
    +1.1711d-007, +1.6732d+000, +1.9858d-001, +2.3975d-002, -9.0013d-004, &
    -2.2475d-003, -3.3095d-005, -1.2040d-005, +2.2010d-006, -1.0083d-006, &
    +8.6297d-001, +5.8231d-001, +2.0545d-002, -7.8110d-003, -1.4085d-004, &
    -8.8459d-006, +5.7256d-006, -1.5068d-006, +4.0095d-007, -2.4185d-008/

  data(b_geoid(i), i=1, 55)/ &
    +0.0000d+000, +0.0000d+000, -6.5993d-002, +0.0000d+000, +6.5364d-002, &
    -5.8320d+000, +0.0000d+000, +1.6961d+000, -1.3557d+000, +1.2694d+000, &
    +0.0000d+000, -2.9310d+000, +9.4805d-001, -7.6243d-002, +4.1076d-002, &
    +0.0000d+000, -5.1808d-001, -3.4583d-001, -4.3632d-002, +2.2101d-003, &
    -1.0663d-002, +0.0000d+000, +1.0927d-001, -2.9463d-001, +1.4371d-003, &
    -1.1452d-002, -2.8156d-003, -3.5330d-004, +0.0000d+000, +4.4049d-001, &
    +5.5653d-002, -2.0396d-002, -1.7312d-003, +3.5805d-005, +7.2682d-005, &
    +2.2535d-006, +0.0000d+000, +1.9502d-002, +2.7919d-002, -8.1812d-003, &
    +4.4540d-004, +8.8663d-005, +5.5596d-005, +2.4826d-006, +1.0279d-006, &
    +0.0000d+000, +6.0529d-002, -3.5824d-002, -5.1367d-003, +3.0119d-005, &
    -2.9911d-005, +1.9844d-005, -1.2349d-006, -7.6756d-009, +5.0100d-008/

! degree n and order m
  nmax = 9
  mmax = 9

! unit vector
  x = DCOS(dlat)*DCOS(dlon)
  y = DCOS(dlat)*DSIN(dlon)
  z = DSIN(dlat)

! Legendre polynomials
  V(1, 1) = 1.0D0
  W(1, 1) = 0.0D0; 
  V(2, 1) = z*V(1, 1); 
  W(2, 1) = 0.0; 
  DO n = 2, nmax
    V(n + 1, 1) = ((2*n - 1)*z*V(n, 1) - (n - 1)*V(n - 1, 1))/n
    W(n + 1, 1) = 0.0D0
  ENDDO

  DO m = 1, nmax
    V(m + 1, m + 1) = (2*m - 1)*(x*V(m, m) - y*W(m, m))
    W(m + 1, m + 1) = (2*m - 1)*(x*W(m, m) + y*V(m, m))
    IF (m < nmax) THEN
      V(m + 2, m + 1) = (2*m + 1)*z*V(m + 1, m + 1)
      W(m + 2, m + 1) = (2*m + 1)*z*W(m + 1, m + 1)
    ENDIF
    DO n = m + 2, nmax
      V(n + 1, m + 1) = ((2*n - 1)*z*V(n, m + 1) - (n + m - 1)*V(n - 1, m + 1))/(n - m)
      W(n + 1, m + 1) = ((2*n - 1)*z*W(n, m + 1) - (n + m - 1)*W(n - 1, m + 1))/(n - m)
    ENDDO
  ENDDO

! Geoidal height
  undu = 0.d0
  i = 0
  DO n = 0, nmax
    DO m = 0, n
      i = i + 1
      undu = undu + (a_geoid(i)*V(n + 1, m + 1) + b_geoid(i)*W(n + 1, m + 1))
    ENDDO
  ENDDO

  return
END SUBROUTINE
