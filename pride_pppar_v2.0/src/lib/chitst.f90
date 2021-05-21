!
!! chitst.f90
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
!! Contributor: Jianghui Geng, Songfeng Yang
!! 
!!
!!
!! purpose  : chi-square test
!! parameter:
!!    input : flag -- 0: dual-side test
!!                   -1: left side test, H0 nsig^2 < osig^2   H1 nsig^2 >= osig^2
!!                    1: right side test, H0 nsig^2 > osig^2   H1 nsig^2 <= osig^2
!!            ndof -- degree of freedom
!!            osig -- reference sigma
!!            nsig -- new sigma
!!            zone -- belief zone 0.95, 0.99 ...
!!    output: chitst -- true, then reject H0, otherwise false, accept H0
!! NOTE     : go to page 206 of Chinese book of statistics (Edition of Gaojiao)
!
logical*1 function chitst(flag, ndof, osig, nsig, zone)
  implicit none

  integer*4 flag, ndof
  real*8 osig, nsig, zone
!
!! local
  integer*4 ierr
  real*8 chi2l, chi2r, afa, stat
!
!! chi test
  stat = ndof*(nsig/osig)**2
  chitst = .false.
  if (flag .eq. 0) then
    afa = (1.d0 - zone)/2.d0
    call pchi2(ndof, afa, 1, chi2r)
    afa = 1.d0 - afa
    call pchi2(ndof, afa, 1, chi2l)
    if (stat .gt. chi2l .and. stat .lt. chi2r) chitst = .true.
  else if (flag .gt. 0) then
    afa = 1.d0 - zone
    call pchi2(ndof, afa, 1, chi2r)
    if (stat .lt. chi2r) chitst = .true.
  else
    afa = zone
    call pchi2(ndof, afa, 1, chi2l)
    if (stat .gt. chi2l) chitst = .true.
  endif

  return
end
