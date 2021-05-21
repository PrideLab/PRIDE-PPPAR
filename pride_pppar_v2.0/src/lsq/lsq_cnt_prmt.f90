!
!! lsq_cnt_prmt.f90
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
!!! purpose   : count number of the three types of parameters
!!             process, state and determinated parameters
!!
!
subroutine lsq_cnt_prmt(LCF, SITE, NM)
  implicit none
  include '../header/const.h'
  include '../header/station.h'
  include 'lsqcfg.h'
  include 'lsq.h'

  type(station) SITE
  type(lsqcfg) LCF
  type(norm) NM
!
!! local
  integer*4 i

  NM%nc = 0
  NM%np = 0
  NM%ns = 0
!
!! STATION PARAMETERS
!! We only estimate position, troposphere and clock parameters for stations
!
!! STATIC means only estimating 3 positions as constant parameters
  if (SITE%skd(1:1) .eq. 'S') then
    NM%nc = NM%nc + 3
!! K means more state paramters to be estimated including x as p parameters
  else if (SITE%skd(1:1) .eq. 'K') then
    NM%np = NM%np + 3
  else if (SITE%skd(1:1) .ne. 'F') then
    write (*, '(2a,a2)') '***ERROR(lsq_cnt_prmt): unknown site type ', SITE%name, SITE%skd
    call exit(1)
  endif
!
!! atmospheric parameters is process parameters
  if (LCF%ztdmod(1:3) .eq. 'PWC' .or. LCF%ztdmod(1:3) .eq. 'STO') then
    NM%np = NM%np + 1
  else if (LCF%ztdmod(1:3) .ne. 'NON' .and. LCF%ztdmod(1:3) .ne. 'FIX') then
    write (*, '(2a)') '***ERROR(lsq_cnt_prmt): ztd mode ', LCF%ztdmod
    call exit(1)
  endif
  if (LCF%htgmod(1:3) .eq. 'PWC' .or. LCF%htgmod(1:3) .eq. 'STO') then
    NM%np = NM%np + 2
  else if (LCF%htgmod(1:3) .ne. 'NON' .and. LCF%htgmod(1:3) .ne. 'FIX') then
    write (*, '(2a)') '***ERROR(lsq_cnt_prmt): htg mode ', LCF%htgmod
    call exit(1)
  endif
!
!! receiver clock offset
  do i=1,LCF%sysnum
    NM%np = NM%np + 1
  enddo

  NM%imtx = NM%nc + NM%np
!
!! check consistence
  if (NM%imtx .gt. MAXPAR) then
    write (*, '(a,i8)') '***ERROR(lsq_cnt_prmt): too many parameters ', NM%imtx
    call exit(1)
  endif

  return
end
