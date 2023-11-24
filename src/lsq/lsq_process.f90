!
!! lsq_process.f90
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
!! purpose   : time update of process parameter in LSQ
!! parameters: lfncid, lfnrem -- lsq structure
!!             jd,sod,dintv   -- time tag
!!             NM,PM          -- normal matrix & PAR table
!! Reference : Ge Maorong  PhD Dissertation
!!       v= p(i+1)-M(i+1)*p(i)+w(i+1), Pw(i+1)
!
subroutine lsq_process(lfncid, lfnrem, jd, sod, LCF, NM, PM, pdop)
  implicit none
  include '../header/const.h'
  include 'lsqcfg.h'
  include 'lsq.h'

  integer*4 lfncid, lfnrem, jd
  real*8 sod, pdop, pdop_sav
  type(lsqcfg) LCF
  type(norm) NM
  type(prmt) PM(1:*)
!
!! local variables
  integer*4 sec, ipar, k

  save pdop_sav

!
!! one each time
  do ipar = NM%nc + 1, NM%nc + NM%np
    if (PM(ipar)%ptype .ne. 'P') cycle
    k = index(PM(ipar)%pname, ':')
    if (k .ne. 0) then
      read (PM(ipar)%pname(k + 1:), *) sec 
      if ((jd*864.d2 + sod + LCF%dintv .lt. PM(ipar)%ptime(2)*864.d2  - MAXWND) .and. &
          (jd*864.d2 + sod             .lt. LCF%jd1*864.d2 + LCF%sod1 - MAXWND)) cycle
    end if
!
!! remove process parameter
    if (PM(ipar)%iobs .gt. 0) then
      call lsq_rmv_prmt(.false., lfncid, lfnrem, ipar, NM, PM, NM%norx)
      if (PM(ipar)%pname(1:5) .eq. 'STAPX') then
        write (lfncid) 'dp'
        write (lfnrem) (pdop - pdop_sav)/PM(ipar)%iepo 
        pdop_sav = pdop
      end if
      PM(ipar)%iepo = 0
      PM(ipar)%iobs = 0
      PM(ipar)%iobs_G = 0
      PM(ipar)%iobs_R = 0
      PM(ipar)%iobs_E = 0
      PM(ipar)%iobs_C = 0
      PM(ipar)%iobs_3 = 0
      PM(ipar)%iobs_J = 0
    end if
    if (k .ne. 0) then
      PM(ipar)%ptime(1) = PM(ipar)%ptime(1) + sec/864.d2
      PM(ipar)%ptime(2) = PM(ipar)%ptime(2) + sec/864.d2
    else
      PM(ipar)%ptime(1) = jd + (sod + LCF%dintv)/864.d2
      PM(ipar)%ptime(2) = jd + (sod + LCF%dintv)/864.d2
    end if
!! next process parameter
  end do

  return
end subroutine
