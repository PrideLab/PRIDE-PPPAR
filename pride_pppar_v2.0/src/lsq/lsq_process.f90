!
!! lsq_process.f90
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
!! purpose   : time update of process parameter in LSQ
!! parameters: lfncid, lfnrem -- lsq structure
!!             jd,sod,dintv   -- time tag
!!             NM,PM          -- normal matrix & PAR table
!! Reference : Ge Maorong  PhD Dissertation
!!       v= p(i+1)-M(i+1)*p(i)+w(i+1), Pw(i+1)
!
subroutine lsq_process(lfncid, lfnrem, jd, sod, dintv, NM, PM)
  implicit none
  include '../header/const.h'
  include 'lsq.h'

  integer*4 lfncid, lfnrem, jd
  real*8 sod, dintv
  type(norm) NM
  type(prmt) PM(1:*)
!
!! local variables
  integer*4 minut, ipar, k
!
!! one each time
  do ipar = NM%nc + 1, NM%nc + NM%np
    if (PM(ipar)%ptype .ne. 'P') cycle
    k = index(PM(ipar)%pname, ':')
    if (k .ne. 0) then
      read (PM(ipar)%pname(k + 1:), *) minut
      if ((jd - PM(ipar)%ptime(2))*86400.d0 + sod .lt. -MAXWND) cycle
    endif
!
!! remove process parameter
    if (PM(ipar)%iobs .gt. 0) then
      call lsq_rmv_prmt(.false., lfncid, lfnrem, ipar, NM, PM, NM%norx)
      PM(ipar)%iobs = 0
    endif
    if (k .ne. 0) then
      PM(ipar)%ptime(1) = PM(ipar)%ptime(1) + minut/1440.d0
      PM(ipar)%ptime(2) = PM(ipar)%ptime(2) + minut/1440.d0
    else
      PM(ipar)%ptime(1) = jd + (sod + dintv)/86400.d0
      PM(ipar)%ptime(2) = jd + (sod + dintv)/86400.d0
    endif
!! next process parameter
  enddo

  return
end
