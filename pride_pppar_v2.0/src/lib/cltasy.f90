!
!! cltasy.f90
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
!! Contributor: Maorong Ge, Jianghui Geng, Songfeng Yang
!! 
!!
!!
!! purpose  : correlation analysis
!! parameter:
!!    input : nobs -- # of observations
!!            nukn -- # of unknown parameters
!!            amat -- design matrix
!!            invn -- inversed normal matrix
!!            resi -- reisudal vector
!!    output: iptx -- pointer to deleted observation
!
subroutine cltasy(nobs, nukn, amat, invn, resi, iptx)
  implicit none
  include '../header/const.h'

  integer*4 iptx, nobs, nukn
  real*8 amat(MAXSAT, nukn), invn(nukn, nukn), resi(MAXSAT)
!
!! local
  integer*4 i, j, k, imax
  real*8 rmtx(MAXSAT, MAXSAT), help(MAXSAT, nukn), uresi(MAXSAT), rhelp(MAXSAT), corr, dummy
!
!! function called
  real*8 dot
!
!! initialization
  iptx = 0
  do j = 1, nobs
    do i = 1, nobs
      rmtx(i, j) = 0.d0
      if (j .le. nukn) help(i, j) = 0.d0
    enddo
  enddo
!
!! formation of R matrix
!.. B * N-1
  do i = 1, nobs
    do j = 1, nukn
      do k = 1, nukn
        help(i, j) = help(i, j) + amat(i, k)*invn(k, j)
      enddo
    enddo
  enddo
!.. B * N-1 * Bt - E
  do i = 1, nobs
    do j = 1, nobs
      do k = 1, nukn
        rmtx(i, j) = rmtx(i, j) + help(i, k)*amat(j, k)
      enddo
      if (i .eq. j) rmtx(i, j) = rmtx(i, j) - 1.d0
    enddo
  enddo
!
!! correlation coefficients
  call unit_vector(nobs, resi, uresi, dummy)
  imax = 0
  corr = -10.d0
  do j = 1, nobs
    call unit_vector(nobs, rmtx(1, j), rhelp, dummy)
    dummy = dabs(dot(nobs, rhelp, uresi))
    if (dummy .gt. corr) then
      imax = j
      corr = dummy
    endif
  enddo
!
!! Hypothesis test needed. Future work
  if (corr .gt. 0.5d0) then
    iptx = imax
  endif

  return
end
