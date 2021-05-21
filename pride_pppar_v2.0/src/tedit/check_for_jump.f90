!
!! check_for_jump.f90
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
!! purpose  : Check jump using polynomial fitting
!!
!! input :
!!        string --------  ' LG ' or ' LGC '
!!        lfnout --------  -1 or 6 or 0
!!        nepo   --------  the number of epoch in check arc
!!        ti     --------  first epoch time
!!        li     --------  geometry-free obs value
!!        flg    --------  10 : initialization
!!                          1 : "ok"
!!                          0 : "good"
!!        ndgr   --------   2 : default value
!!        niter  --------   3 : default value
!!        bias0  --------  2.0: default value
!!        sig0   --------  0.4: default value
!! output:
!!        a0     --------
!!        rms    --------
!!        v      --------
!!        ipt    --------
!
subroutine check_for_jump(string, lfnout, nepo, ti, li, flg, ndgr, niter, bias0, sig0, a0, rms, v, ipt, ierr, interval)
  implicit none
  integer*4 nepo, niter, ndgr, flg(1:*), ierr, lfnout, ipt(1:*)
  real*8 ti(1:*), li(1:*), sig0, bias0, interval, a0, rms, v(1:*)
  character*(*) string
!
!! local
  integer*4 i, k, iter, nflg, iepo0, idgr
  real*8 coeff(ndgr*2), dx, ft, prerms
  logical*1 again, lwrite

  lwrite = .true.
  if (lfnout .eq. 0) lfnout = 6
  if (lfnout .lt. 0) lwrite = .false.

  iepo0 = nint(dmod(ti(1), 86400.d0)/interval)
  ierr = 0
  iter = 0
  idgr = 0
  prerms = 1.d10
  again = .true.
  a0 = 0.d0
  do while (again)
! polynomial fit
    call polydf(ti, li, flg, nepo, idgr, coeff, v, rms, ipt, dx, ft, ierr)
    if (ierr .ne. 0) then
      if (lwrite) write (lfnout, *) ' ploydf error ', ierr
      return
    endif
    a0 = coeff(1)*ft
    if (lwrite) write (lfnout, '(a,4i6,f12.2)') string//' fitrms,ndgr,iter,', iepo0 + 1, iepo0 + nepo, idgr, iter, rms
    nflg = 0
    do k = 1, nepo
      if (ipt(k) .ne. 0) then
        if (dabs(v(k)) .gt. max(3*rms, bias0)) then
          if (flg(k) .eq. 0) then
            if (lwrite) write (lfnout, '(a,i6,f15.3,4f10.3)') string//' flagged(b):', k + iepo0, ti(k), v(k), bias0, rms
            nflg = nflg + 1
            flg(k) = 1
          endif
        endif
      endif
    enddo
    iter = iter + 1
    again = .false.
    if (nflg .ne. 0) then
      again = .true.
    else if (rms .gt. sig0) then
      again = .true.
      if (rms/prerms .gt. 0.95d0) then
        idgr = idgr + 1      ! increase degrees gradually
        iter = 1
      endif
      if (idgr .gt. ndgr) again = .false.
    endif
    prerms = rms
    if (again .and. iter .le. niter) again = .true.
  enddo
  if (lwrite) write (lfnout, '(a,4i6,f12.2)') string//' fitrms,ndgr,iter,', iepo0 + 1, iepo0 + nepo, idgr, iter, rms

  return
end
