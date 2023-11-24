!
!! lsq_cnt_prmt.f90
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
!! Contributor: Maorong Ge, Jianghui Geng, Songfeng Yang, Jihang Lin, Jing Zeng
!!
!!
!!
!!! purpose   : count number of the three types of parameters
!!              process, state and determinated parameters
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
  if (SITE%skd(1:1) .eq. 'S' .or. SITE%skd(1:1) .eq. 'F') then
    NM%nc = NM%nc + 3
  else if (SITE%skd(1:1) .eq. 'P') then
    NM%np = NM%np + 3
  else if (SITE%skd(1:1) .eq. 'K') then
    NM%np = NM%np + 3
  else if (SITE%skd(1:1) .eq. 'L') then
    NM%np = NM%np + 3
  else
    write (*, '(2a,1x,a)') '***ERROR(lsq_cnt_prmt): unknown positioning mode ', SITE%name, SITE%skd
    call exit(1)
  end if
!
!! atmospheric parameters is process parameters
  if (LCF%ztdmod(1:3) .eq. 'PWC' .or. LCF%ztdmod(1:3) .eq. 'STO') then
    NM%np = NM%np + 1
  else if (LCF%ztdmod(1:3) .ne. 'NON') then
    write (*, '(2a)') '***ERROR(lsq_cnt_prmt): ztd mode ', LCF%ztdmod
    call exit(1)
  end if
  if (LCF%htgmod(1:3) .eq. 'PWC' .or. LCF%htgmod(1:3) .eq. 'STO') then
    NM%np = NM%np + 2
  else if (LCF%htgmod(1:3) .ne. 'NON') then
    write (*, '(2a)') '***ERROR(lsq_cnt_prmt): htg mode ', LCF%htgmod
    call exit(1)
  end if
!
!! receiver clock offset
  NM%np = NM%np + len_trim(LCF%sys)

  NM%imtx = NM%nc + NM%np
!
!! check consistence
  if (NM%imtx .gt. MAXPAR) then
    write (*, '(a,i8)') '***ERROR(lsq_cnt_prmt): too many parameters ', NM%imtx
    call exit(1)
  end if

  return
end subroutine
