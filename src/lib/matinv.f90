!!
!! matinv.f90
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
!! purpose    :  inversion of a matrix by means of gauss....
!!
!! parameters : a : matrix to be inverted
!!           ndim : dimension of a matrix
!!              n : dimension of a to be invert
!!              d : value of determinant a
!!

subroutine matinv(a, ndim, n, d)

  implicit none
  integer*4 ndim, n
  real*8 a(ndim, ndim), d
!
!! local
  integer*4 i, j, k, l(ndim), m(ndim)
  real*8 max_a, swap

  d = 1.d0
!
!! loop over all row
  do k = 1, n
    l(k) = k
    m(k) = k
    max_a = a(k, k)
!
!! find out the diag. element with max. abs. value.
    do i = k, n
      do j = k, n
        if (dabs(max_a) - dabs(a(i, j)) .lt. 0.d0) then
          max_a = a(i, j)
          m(k) = i
          l(k) = j
        endif
      enddo
    enddo
!
!! rank defection or bad condition
    if (dabs(max_a) .lt. 1d-13) then
      d = 0.d0
      return
    endif
!
!! swap k and l(k) column
    if (l(k) .gt. k) then
      do i = 1, n
        swap = -a(i, k)
        a(i, k) = a(i, l(k))
        a(i, l(k)) = swap
      enddo
    endif
!
!! swap k and m(k) raw
    if (m(k) .gt. k) then
      do j = 1, n
        swap = -a(k, j)
        a(k, j) = a(m(k), j)
        a(m(k), j) = swap
      enddo
    endif
!
!! elemination
    do i = 1, n
      if (i .ne. k) a(k, i) = -a(k, i)/max_a
    enddo

    do i = 1, n
      do j = 1, n
        if (i .ne. k .and. j .ne. k) a(j, i) = a(j, i) + a(k, i)*a(j, k)
      enddo
    enddo
!
!! save invert part
    do j = 1, n
      if (j .ne. k) a(j, k) = a(j, k)/max_a
    enddo
    a(k, k) = 1/max_a
  enddo
!
!! re-range raw and column as input
  do k = n, 1, -1
    if (l(k) .gt. k) then
      do j = 1, n
        swap = a(k, j)
        a(k, j) = -a(l(k), j)
        a(l(k), j) = swap
      enddo
    endif

    if (m(k) .gt. k) then
      do i = 1, n
        swap = a(i, k)
        a(i, k) = -a(i, m(k))
        a(i, m(k)) = swap
      enddo
    endif
  enddo
  return
end
