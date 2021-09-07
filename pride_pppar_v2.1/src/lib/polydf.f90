!
!! polydf.f90
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
!! Contributor: Maorong Ge
!! 
!!
!!
!! purpose  : approximation of a set of discrete functional values
!!            by a polynomial and vdetermining the integral of the
!!            polynomial using differences of consecutive values
!! parameter:
!!   input  :     x --- vector with independent arguments
!!                y --- vector with functional values
!!                l --- l(i)=0: value # i is used for approx.
!!                          =1: value # i is not used for approx.
!!                              value # i is not used at all
!!                n --- # of values (=dimension of x,y,v)
!!               ip --- degree of polynomial
!!  output  :
!!                c --- vector with coefficients of polynomial
!!                v --- vector with residuals (diff approx-y)
!!            dx,ft --- normalisation of ind. argument
!!                                 xnorm=(x-dx)*ft
!!             ierr --- flag indicating not enough obs if != 0
!! restict  :  maximum degree== 14
subroutine polydf(x, y, l, n, ip, c, v, rms, ipt, dx, ft, ierr)
  implicit none

  integer*4 l(1:*), ipt(1:*), l_save(n)
  real*8 x(1:*), y(1:*), c(1:*), v(1:*), rms, dx, ft
  integer*4 n, ip, ierr
!! local
  real*8 xmin, xmax, ys, xnorm, xnorj, det, vmax, rms1
  integer*4 i, j, k, m, neff, ik, ki
  real*8 q(ip, ip), f(ip), a(ip)

  ierr = 0
!! normalisation of independent argument
  xmin = 1.d20
  xmax = -1.d20
!! searching the Max.  and  min.  argument values
  do i = 1, n
    l_save(i) = l(i)
    ipt(i) = 0
    if (l(i) .le. 1) then
      if (x(i) .gt. xmax) xmax = x(i)
      if (x(i) .lt. xmin) xmin = x(i)
      v(i) = 0.d0
    endif
  enddo

  if (xmin .eq. xmax) then
    ierr = 1
    return
  endif
  dx = (xmax + xmin)/2.d0
  ft = 2.d0/(xmax - xmin)
!! initial importent parameter
  do i = 1, ip
    do j = 1, ip
      q(i, j) = 0.d0
    enddo
    f(i) = 0.d0
  enddo

  neff = 0
!! loop over all functional values
  j = n
  do while (j .ge. 2)
!!  marked values
    if (l(j) .eq. 0) then

!!  searching next usable functional value to form difference
      m = j - 1
      do while (m .ge. 1 .and. l(m) .gt. 1)
        m = m - 1
      enddo
      if (m .ne. 0) then
        xnorj = (x(j) - dx)*ft
        xnorm = (x(m) - dx)*ft
        neff = neff + 1

! computaion of a-matrix (row number j)
        do k = 1, ip
          a(k) = xnorm**k - xnorj**k
        enddo

! contribution to q-matrix (upper triangle)
        do i = 1, ip
          do k = i, ip
            q(i, k) = q(i, k) + a(i)*a(k)
          enddo
        enddo

! contribution to right hand side
        do i = 1, ip
          f(i) = f(i) + a(i)*(y(m) - y(j))
        enddo
      endif
    endif
! next difference
    j = j - 1
  enddo

  if (neff .lt. ip) then
!  write(*,*) neff,ip
    ierr = 2
    return
  endif

! lower triangle of q matrix
  do i = 1, ip
    do k = 1, i
      q(i, k) = q(k, i)
    enddo
  enddo

! inversion of q matrix
  call matinv(q, ip, ip, det)
  if (det .eq. 0) then
    ierr = 3
    return
  endif

! solution vector (polynomial coefficients)
  do i = 1, ip
    c(i) = 0.d0
    do k = 1, ip
      c(i) = c(i) + q(i, k)*f(k)
    enddo
  enddo

!  residuals, rms
  rms = 0.d0
  vmax = 0.d0
  j = n
  do while (j .ge. 2)
    if (l(j) .le. 1) then
      m = j - 1
      do while (m .ge. 1 .and. l(m) .gt. 1)
        m = m - 1
      enddo
      if (m .ne. 0) then
        xnorj = (x(j) - dx)*ft
        xnorm = (x(m) - dx)*ft
        ys = 0.d0
        do k = 1, ip
          ys = ys + c(k)*(xnorm**k - xnorj**k)
        enddo
        v(j) = ys - (y(m) - y(j))
        if (l(j) .eq. 0) then
          if (dabs(v(j)) .gt. vmax) then
            vmax = dabs(v(j))
          endif
          rms = rms + v(j)*v(j)
        endif
        ipt(j) = m
      endif
    endif
    j = j - 1
  enddo

  if (neff - ip .ne. 0 .and. rms .gt. 0.d0) then
    vmax = rms - vmax*vmax
    rms = dsqrt(rms/(neff - ip))
    if (neff - ip - 1 .ne. 0 .and. vmax .gt. 0.d0) then
      rms1 = dsqrt(vmax/(neff - ip - 1))
      if (rms/rms1 .gt. 2.d0) rms = rms1
    endif
  else
!  write(*,*) ' BAD rms ',rms,neff,ip,n
!  do j=2,n
!    write(*,*) n,j,l(j),l_save(j),ipt(j),y(j),v(j)
!  enddo
    ierr = 4
  endif

  return
end
