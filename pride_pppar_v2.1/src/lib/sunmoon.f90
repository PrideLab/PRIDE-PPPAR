!
!! sunmoon.f90
!!
!!    Copyright (C) 2021 by Wuhan University
!!
!!    This program belongs to PRIDE PPP-AR which is an open source software:
!!    you can redistribute it and/or modify it under the terms of the GNU
!!    General Public License (version 3) as published by the Free Software Foundation.
!!
!!    This program is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License (version 3) for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with this program.  If not, see <https://www.gnu.org/licenses/>.
!!
!! Contributor: Ran Zeng
!! 
!!
!!
subroutine getghar(jd, sod, ghar)
  implicit none
  include "../header/const.h"
  real*8 sod, ghar, D2R, tsecgps, tsecutc, fmjdutc, d, ghad
  integer*4 i, jd
  
  D2R=PI/180.d0
  !/* need UT to get sidereal time */
  tsecgps = sod  !/* GPS time (sec of day)           */
  tsecutc = sod-18  !/* UTC time (sec of day)           */
  fmjdutc = tsecutc / 86400.0  !/* UTC time (fract. day)           */
  d = (jd - 51544) + (fmjdutc - 0.50)  !/* days since J2000                */
  ghad = 280.46061837504 + 360.9856473662862 * d  !/* corrn.   (+digits)         */
  i = int((ghad / 360.0))
  ghar = (ghad - i * 360.0) * D2R
  
  do while(ghar >= 2 * PI)
    ghar = ghar - 2 * PI
  enddo
  do while(ghar < 0.0)
    ghar = ghar + 2 * PI;
  enddo
end

subroutine rot3(theta, x, y, z, u, v, w)
  implicit none
  real*8 theta, x, y, z, u, v, w, s, c
  
  s = dsin(theta)
  c = dcos(theta)
  u = c * x + s * y
  v = c * y - s * x
  w = z
end

subroutine sunxyz(jd, sod, rs)
  implicit none
  include "../header/const.h"
  real*8 sod, rs(3), D2R
  real*8 obe, sobe, cobe, opod
  real*8 tsecgps, tsectt, fmjdtt
  real*8 tjdtt, t, emdeg, em, em2
  real*8 r, slond
  real*8 slon, sslon, cslon, rs1, rs2, rs3
  real*8 ghar
  integer*4 jd
  
  D2R=PI/180.0
  !/* mean elements for year 2000, sun ecliptic orbit wrt. Earth */
  obe     = 23.43929111 * D2R  !/* obliquity of the J2000 ecliptic */
  sobe    = dsin(obe)
  cobe    = dcos(obe)
  opod    = 282.94
  
  !/* use TT for solar ephemerides */
  tsecgps = sod  !/* GPS time (sec of day)           */
  tsectt  = sod+19.0+32.184  !/* TT  time (sec of day)           */
  fmjdtt  = tsectt / 86400.0
  
  !/* julian centuries since 1.5 january 2000 (J2000) */
  !/*    (note also low precision use of mjd --> tjd) */
  tjdtt   = jd + fmjdtt + 2400000.5  !/* Julian Date, TT                 */
  t       = (tjdtt - 2451545.0) / 36525.0  !/* Julian centuries, TT            */
  emdeg   = 357.5256 + 35999.049 * t  !/* degrees                         */
  em      = emdeg * D2R  !/* radians                         */
  em2     = em + em
  
  !/* series expansions in mean anomaly, em */
  r       = (149.619 - 2.499 * cos(em) - 0.021 * cos(em2)) * 1.0e9  !/* m.                              */
  slond   = opod + emdeg + (6892.0 * sin(em) + 72.0 * sin(em2)) / 3600.0

  !/* precession of equinox */
  slond         = slond + 1.3972 * t    !/* degrees                         */

  !/* position vector of sun (mean equinox & ecliptic of J2000) (EME2000, ICRF) */
  !/*                                    (plus long. advance due to precession) */
  slon    = slond * D2R  !/* radians                         */
  sslon   = sin(slon)
  cslon   = cos(slon)

  rs1     = r * cslon  !/* meters                          */
  rs2     = r * sslon * cobe  !/* meters                          */
  rs3     = r * sslon * sobe  !/* meters                          */

  !/* convert position vector of sun to ECEF (ignore polar motion / LOD) */
  ghar = 0.0
  call getghar(jd, sod, ghar)
  call rot3(ghar, rs1, rs2, rs3, rs(1), rs(2), rs(3))
end

subroutine moonxyz(jd, sod, rm)
  implicit none
  include "../header/const.h"
  real*8 sod, rm(3)
  real*8 tsecgps, tsectt, fmjdtt, D2R
  real*8 tjdtt, t
  real*8 el0, el, elp, f, d
  real*8 selond
  real*8 q, selatd
  real*8 rse
  real*8 oblir, sselat, cselat, sselon, cselon, t1, t2, t3, rm1, rm2, rm3
  real*8 ghar
  integer*4 jd
  
  !/* use TT for lunar ephemerides */
  tsecgps = sod  !/* GPS time (sec of day)           */
  tsectt  = sod+19.0+32.184  !/* TT  time (sec of day)           */
  fmjdtt  = tsectt / 86400.0  !/* TT  time (fract. day)           */
  D2R=PI/180.0
  
  !/* julian centuries since 1.5 january 2000 (J2000) */
  !/*    (note also low precision use of mjd --> tjd) */
  tjdtt   = jd + fmjdtt + 2400000.5  !/* Julian Date, TT                 */
  t       = (tjdtt - 2451545.0) / 36525.0
  
  !/* el0 -- mean longitude of Moon (deg)                             */
  !/* el  -- mean anomaly of Moon (deg)                               */
  !/* elp -- mean anomaly of Sun (sun)                                */
  !/* f   -- mean angular distance of Moon from ascending node (deg)  */
  !/* d   -- difference between mean longitudes of Sun and Moon (deg) */
  el0     = 218.31617 + 481267.88088 * t - 1.3972 * t
  el      = 134.96292 + 477198.86753 * t
  elp     = 357.52543 +  35999.04944 * t
  f       =  93.27283 + 483202.01873 * t
  d       = 297.85027 + 445267.11135 * t
  
  !/* longitude w.r.t. equinox and ecliptic of year 2000 */
  selond  = el0 &
     + ((22640.0 * dsin((el              ) * D2R)) &
     +  (  769.0 * dsin((el + el         ) * D2R)) &
     -  ( 4586.0 * dsin((el - d - d      ) * D2R)) &
     +  ( 2370.0 * dsin((d + d           ) * D2R)) &
     -  (  668.0 * dsin((elp             ) * D2R)) &
     -  (  412.0 * dsin((f + f           ) * D2R)) &
     -  (  212.0 * dsin((el + el - d - d ) * D2R)) &
     -  (  206.0 * dsin((el + elp - d - d) * D2R)) &
     +  (  192.0 * dsin((el + d + d      ) * D2R)) &
     -  (  165.0 * dsin((elp - d - d     ) * D2R)) &
     +  (  148.0 * dsin((el - elp        ) * D2R)) &
     -  (  125.0 * dsin((d               ) * D2R)) &
     -  (  110.0 * dsin((el + elp        ) * D2R)) &
     -  (   55.0 * dsin((f + f - d - d   ) * D2R))) / 3600.0
     
  !/* latitude w.r.t. equinox and ecliptic of year 2000 */
  q      = (412.0 * dsin((f + f) * D2R) + 541.0 * dsin((elp) * D2R)) / 3600.0
  selatd = &
    + ((18520.0 * dsin((f + selond - el0 + q) * D2R)) &
    - (   526.0 * dsin((f - d - d           ) * D2R)) &
    + (    44.0 * dsin((el + f - d - d      ) * D2R)) &
    - (    31.0 * dsin((-el + f - d - d     ) * D2R)) &
    - (    25.0 * dsin((-el - el + f        ) * D2R)) &
    - (    23.0 * dsin((elp + f - d - d     ) * D2R)) &
    + (    21.0 * dsin((-el + f             ) * D2R)) &
    + (    11.0 * dsin((-elp + f - d - d    ) * D2R))) / 3600.0

  !/* distance from Earth center to Moon (m) */
  rse    = (385000.0 &
    - ( 20905.0 * dcos((el              ) * D2R)) &
    - (  3699.0 * dcos((d + d - el      ) * D2R)) &
    - (  2956.0 * dcos((d + d           ) * D2R)) &
    - (   570.0 * dcos((el + el         ) * D2R)) &
    + (   246.0 * dcos((el + el - d - d ) * D2R)) &
    - (   205.0 * dcos((elp - d - d     ) * D2R)) &
    - (   171.0 * dcos((el + d + d      ) * D2R)) &
    - (   152.0 * dcos((el + elp - d - d) * D2R))) * 1000.0

  !/* convert spherical ecliptic coordinates to equatorial cartesian */
  !/* precession of equinox wrt. J2000 */
  selond       = selond + 1.3972 * t  !/* degrees                         */

  !/* position vector of moon (mean equinox & ecliptic of J2000) (EME2000, ICRF) */
  !/*                                     (plus long. advance due to precession) */
  oblir  = 23.43929111 * D2R  !/* obliquity of the J2000 ecliptic */

  sselat = dsin(selatd * D2R)
  cselat = dcos(selatd * D2R)
  sselon = dsin(selond * D2R)
  cselon = dcos(selond * D2R)

  t1     = rse * cselon * cselat  !/* meters                         */
  t2     = rse * sselon * cselat  !/* meters                         */
  t3     = rse *          sselat  !/* meters                         */
  
  rm1    = 0.0
  rm2    = 0.0
  rm3    = 0.0
  call rot1(-oblir, t1, t2, t3, rm1, rm2, rm3)

  !/* convert position vector of moon to ECEF (ignore polar motion / LOD) */
  ghar = 0.0
  call getghar(jd, sod, ghar)
  call rot3(ghar, rm1, rm2, rm3, rm(1), rm(2), rm(3))
end

subroutine rot1(theta, x, y, z, u, v, w)
  implicit none
  real*8 theta, x, y, z, u, v, w, s, c
  
  s = dsin(theta)
  c = dcos(theta)
  u = x
  v = c*y+s*z
  w = c*z-s*y
end
