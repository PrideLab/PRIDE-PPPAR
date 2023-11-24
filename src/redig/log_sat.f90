!
!! log_sat.f90
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
!! Contributor: Jianghui Geng
!
subroutine log_sat(RCF,nepo)
implicit none
include '../header/const.h'
include 'rescfg.h'

type(rescfg) RCF
integer*4 nepo
!
!! local
integer*4 iy, imon, id, ih, im, jdb, iepo, isat, ierr
real*8 sec,tb,intv
character*256 line
character*3 iprn
!
!! function
integer*4 modified_julday, pointer_string
real*8 timdif,cutemod

do while (.true.)
  read (RCF%lfnrhd, '(a)', iostat=ierr) line
  if (ierr .ne. 0) exit
  if (line(1:3) .ne. 'TIM') cycle
  read (line, '(3x,i5,4i3,f11.7)') iy, imon, id, ih, im, sec
  jdb = modified_julday(id, imon, iy)
  tb  = ih*3600.d0 + im*60.d0 + sec
  intv=timdif(jdb,tb,RCF%jd0,RCF%sod0)
  if(cutemod(intv,RCF%dintv).ne.0.d0) cycle
  iepo = nint(intv/RCF%dintv) + 1
  if (iepo .gt. nepo) then
     write (*, '(a,i5,f10.2)') '%%%MESSAGE(log_sat): log file truncated at ',jdb,tb
     exit
  endif
  do while (.true.)
    read (RCF%lfnrhd, '(a)', iostat=ierr) line
    if (ierr .ne. 0) exit
    if (line(1:3) .eq. 'TIM') then
      backspace RCF%lfnrhd
      exit
    endif
    if (line(61:63) .eq. 'AMB') cycle
    read (line, *) iprn
    isat = pointer_string(RCF%nprn, RCF%prn, iprn)
    if (isat .eq. 0) then
      RCF%nprn=RCF%nprn+1
      RCF%prn(RCF%nprn)=iprn
    endif
  enddo
enddo
rewind RCF%lfnrhd

return
end 
