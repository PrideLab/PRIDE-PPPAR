!
!! everett_coeff.f90
!!
!!    Copyright (C) 2021 by Wuhan University
!!
!!    This program is an open source software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License (version 3) as
!!    published by the Free Software Foundation.
!!
!!    This program is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License (version 3) for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with this program.  If not, see <https://www.gnu.org/licenses/>.
!!
!! Contributor: Maorong Ge
!! 
!!
!!
!! purpose  : compute coefficients of the everett interpolation
!!
!! parameter: n  -- degree of the interpolation, with 2*n+1 points
!!            ec -- everett coefficients
!!
!
subroutine everett_coeff(ndim, n, ec)
  implicit none
  integer*4 n, ndim
  real*8 ec(0:ndim, 0:ndim)
!
!! local
  integer*4 i, j, k
  real*8 a(0:n, 0:n), factor(0:2*n + 1), sign

  factor(0) = 1.d0
  do i = 1, 2*n + 1
    factor(i) = factor(i - 1)*i
  enddo

!
!! coefficients of general bi-?
  a(0, 0) = 1.d0
  do i = 1, n
    a(i, 0) = -a(i - 1, 0)*i*i
    do j = 1, i - 1
      a(i, j) = a(i - 1, j - 1) - i*i*a(i - 1, j)
    enddo
    a(i, i) = 1.d0
  enddo
!
!! everett coefficient
  do k = 0, n
    do j = 0, n
      ec(k, j) = 0.d0
      do i = max(k, j), n
        if (mod(i + j, 2) .eq. 1) then
          sign = -1.d0
        else
          sign = 1.d0
        endif
        ec(k, j) = ec(k, j) + sign*a(i, k)/(2.d0*i + 1.d0)/factor(i + j)/factor(i - j)
      enddo
      if (j .eq. 0) ec(k, j) = ec(k, j)/2.d0
    enddo
  enddo

  return
end
