!
!! lsq_rmv_prmt.f90
!!
!!    Copyright (C) 2018 by Wuhan University
!!
!!    This program is free software: you can redistribute it and/or modify
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
!! Contributor: Maorong Ge, Jianghui Geng
!! 
!!
!!
!! purpose   : pre-eliminate one selected parameter
!! parameter :
!!    input  : lcompact -- compact ntxal matrix or not
!!             lfnrem   -- tmp rem file for recovery
!!             ipar     -- index of selected parameter in PAR table
!!             mp,pw    -- for process parameters
!!             NM,PM    -- normal matrix & PAR table
!!             ntx      -- one-dimensional normal matrix
!
subroutine lsq_rmv_prmt(lcompact, lfncid, lfnrem, ipar, NM, PM, ntx)
  implicit none
  include '../header/const.h'
  include 'lsq.h'

  logical*1 lcompact
  integer*4 lfncid, lfnrem, ipar
  real*8 ntx(1:*)
  type(norm) NM
  type(prmt) PM(1:*)
!
!! local
  integer*4 i, j, k, iptn, nelem, ipt(0:MAXPAR + 1), ipx(0:MAXPAR)
  real*8 mp, pw, elem(MAXPAR + 1), elem0(MAXPAR + 1)
!
!! initialization
  mp = 0.d0; pw = 0.d0
  if (.not. lcompact) then
    mp = PM(ipar)%map
    pw = PM(ipar)%rw**2
  endif
  nelem = 0
  do i = 1, NM%imtx + 1
    elem(i) = 0.d0
    elem0(i) = 0.d0
  enddo
!
!! get row and column to be removed
  iptn = PM(ipar)%ipt   ! index to normal matrix column of par.
  ipt(0) = 0            ! index to normal matrix column of removing elements
  ipx(0) = 0            ! index to par. array of each removing elements
  do i = 1, NM%imtx + 1
    if (i .le. iptn) then
      k = NM%iptx(iptn) + i
    else
      k = NM%iptx(i) + iptn
    endif
    if (ntx(k) .eq. 0.d0 .and. i .ne. NM%imtx + 1) cycle
    nelem = nelem + 1
    ipt(nelem) = i
    if (i .ne. NM%imtx + 1) ipx(nelem) = NM%iptp(i)
    if (i .eq. iptn) then
      ipt(0) = nelem
      ipx(0) = nelem
    endif
    elem(nelem) = ntx(k)
    ntx(k) = 0.d0
  enddo
  if (ipt(0) .eq. 0) then
    write (*, '(2a,2i4)') '***ERROR(lsq_rmv_prmt): zero pivot ', PM(ipar)%pname, PM(ipar)%pcode(1:2)
    call exit(1)
  endif
!
!! check pivot only for process parameters
  if (PM(ipar)%iobs .lt. 1) then
    if (.not. lcompact) then
      ntx(NM%iptx(iptn) + iptn) = pw
      return
    endif
  endif
!
!! pivot adding
  elem(ipt(0)) = elem(ipt(0)) + mp*pw*mp
!
!! ltpl
  NM%ltpl = NM%ltpl - elem(nelem)*elem(nelem)/elem(ipt(0))
!
!! write removed par. into rem file
  if (lfncid .ne. 0 .and. lfnrem .ne. 0) then
    if (PM(ipar)%ptype .eq. 'P') then
      write (lfncid) 'pc'
      write (lfnrem) ipar, nelem, ipx(0), (ipx(i), i=1, nelem - 1), &
        (elem(i), i=1, nelem), PM(ipar)%iobs, PM(ipar)%xini, PM(ipar)%ptime(1:2)
    else if (PM(ipar)%ptype .eq. 'S') then
      write (lfncid) 'am'
      write (lfnrem) ipar, nelem, ipx(0), (ipx(i), i=1, nelem - 1), (elem(i), i=1, nelem)
    endif
  endif
!
!! remove related row and column
  do j = 1, nelem
    elem0(j) = elem(j)/elem(ipt(0))
    if (ipt(j) .eq. iptn) cycle
    do i = 1, j
      if (ipt(i) .eq. iptn .or. ipt(i) .eq. NM%imtx + 1) cycle
      ntx(NM%iptx(ipt(j)) + ipt(i)) = ntx(NM%iptx(ipt(j)) + ipt(i)) - elem0(j)*elem(i)
    enddo
  enddo
!
!! add state transition
  if (.not. lcompact) then
    if (mp .eq. 0.d0) then
      ntx(NM%iptx(iptn) + iptn) = pw
    else
      do i = 1, nelem
        if (ipt(i) .lt. iptn) then
          ntx(NM%iptx(iptn) + ipt(i)) = mp*pw*elem0(i)
        else if (ipt(i) .gt. iptn) then
          ntx(NM%iptx(ipt(i)) + iptn) = mp*pw*elem0(i)
        else
          ntx(NM%iptx(iptn) + iptn) = pw - mp*pw/elem(ipt(0))*mp*pw
        endif
      enddo
      NM%nobs = NM%nobs + 1
    endif
    NM%nuk = NM%nuk + 1
  endif
!
!! compact normal matrix
  if (lcompact) then
    do j = iptn + 1, NM%imtx + 1
      k = 0
      do i = 1, j
        if (i .eq. iptn .or. i .eq. NM%imtx + 1) cycle
        k = k + 1
        elem(k) = ntx(NM%iptx(j) + i)
      enddo
      do i = 1, j - 1
        if (i .eq. NM%imtx) then
          ntx(NM%iptx(j - 1) + i) = 0.d0
        else
          ntx(NM%iptx(j - 1) + i) = elem(i)
        endif
      enddo
    enddo
    do i = 1, NM%imtx
      ntx(NM%iptx(NM%imtx + 1) + i) = 0.d0
    enddo
    if (PM(ipar)%ptype .eq. 'P' .and. PM(ipar)%iobs .gt. 0) NM%nuk = NM%nuk + 1
  endif

  return
end
