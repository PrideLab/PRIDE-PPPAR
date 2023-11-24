!
!! precession_angle.f90
!!
!!    Copyright (C) 2023 by Wuhan University
!!
!!    This program belongs to PRIDE PPP-AR which is an open source software:
!!    you can redistribute it and/or modify it under the terms of the GNU
!!    General Public License (version 3) as published by the Free Software Foundation.
!!
!!    This program is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!!    GNU General Public License (version 3) for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with this program. If not, see <https://www.gnu.org/licenses/>.
!!
!! Contributor: Maorong Ge
!! 
!!
!!
!! sub. name : subroutine precession_angle(jd_tdt,sod_tdt,p1,p2,p3,epsilon0)
!!
!! purpose   : compute precesion quatities(4 angles) at epoch t
!!
!! parameters: jd_tdt,sod_tdt, julian day and second of day of epoch t in TDT system
!!             p1,p2,p3,epsilon0, 4 angles in the sequence of the formulas in Ref. 1,
!!                unit is radians.
!!
!!

subroutine precession_angle(jd_tdt, sod_tdt, p1, p2, p3, epsilon0)
  implicit none
  integer*4 jd_tdt
  real*8 sod_tdt, p1, p2, p3, epsilon0
!
!! local
  integer*4 i, j
  logical*1 first
  real*8 coeff(3, 4), coeff40, x(4), sec2rad, tdt_centuries

! data from IERS standard (1992) page 31
  data coeff/2306.2181d0, 0.30188d0, 0.017998d0, &
    2306.2181d0, 1.09468d0, 0.018203d0, &
    2004.3109d0, -0.42665d0, -0.041833d0, &
    -46.8150d0, -0.00059d0, 0.001813d0/
  data coeff40/84381.448d0/
  data first/.true./

  save coeff, coeff40, first
!
!! convert coefficients to match unit rad.
  if (first) then
!! angle second to radian
    sec2rad = datan(1.d0)/162000.d0
    do i = 1, 4
      do j = 1, 3
        coeff(j, i) = coeff(j, i)*sec2rad
      enddo
    enddo
    coeff40 = coeff40*sec2rad
    first = .false.
  endif
!
!! calculate the angles
  tdt_centuries = (jd_tdt + sod_tdt/86400.d0 - 51544.5d0)/36525d0
  do i = 1, 4
    x(i) = coeff(3, i)*tdt_centuries
    do j = 2, 1, -1
      x(i) = (x(i) + coeff(j, i))*tdt_centuries
    enddo
  enddo
  p1 = x(1)
  p2 = x(2)
  p3 = x(3)
  epsilon0 = coeff40 + x(4)

  return
end
