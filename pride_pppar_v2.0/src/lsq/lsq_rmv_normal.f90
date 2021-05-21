!
!! lsq_rmv_normal.f90
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
!! purpose  : eliminate non-active ambiguities in the LSQ estimator
!! parameter: lfncid,lfnrem -- tmp rem file for recovery
!!            ncol,icol     -- index of parameters to be removed
!!            NM,PM         -- normal matrix & PAR table
!
subroutine lsq_rmv_normal(lfncid, lfnrem, ncol, icol, NM, PM)
  implicit none
  include '../header/const.h'
  include 'lsq.h'

  integer*4 lfncid, lfnrem, ncol, icol(1:*)
  type(norm) NM
  type(prmt) PM(1:*)
!
!! local
  integer*4 ipar, ib, i

  do ib = 1, ncol
    ipar = icol(ib)
    if (ipar .ne. NM%iptp(PM(ipar)%ipt)) then
      write (*, '(a)') '***ERROR(lsq_rmv_normal): index not consistent '
      call exit(1)
    endif
    call lsq_rmv_prmt(.true., lfncid, lfnrem, ipar, NM, PM, NM%norx)
    do i = 1, NM%ipm
      if (PM(i)%ipt .gt. PM(ipar)%ipt) PM(i)%ipt = PM(i)%ipt - 1
    enddo
    do i = PM(ipar)%ipt + 1, NM%imtx
      NM%iptp(i - 1) = NM%iptp(i)
    enddo
    if (PM(ipar)%ptype .eq. 'P') then
      PM(ipar)%ipt = -1        ! act as a flag
      NM%np = NM%np - 1
    else
      PM(ipar)%ipt = 0
      NM%ns = NM%ns - 1
    endif
    NM%iptp(NM%imtx) = 0
    NM%imtx = NM%imtx - 1
  enddo

  return
end
