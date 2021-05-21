!
!! rdglschn.f90
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
!! Read Glonass frequency channel number according to slot number

integer*4 function rdglschn(iprn,jd,sod)
implicit none

character*3  iprn
integer*4    jd
real*8       sod
!
!! local
logical*1    lfirst
integer*4    lfn,ierr
character*80 line
!
!! function used
integer*4    get_valid_unit
real*8       timdif

data lfirst,lfn/.true.,0/
save lfirst,lfn

if(lfirst) then
  lfirst=.false.
  lfn   =get_valid_unit(10)
  open(lfn,file='glonass_chn',status='old',iostat=ierr)
  if(ierr.ne.0) then
    write(*,'(a)') '***ERROR(rdglschn): glonass_chn isn''t found'
    call exit(1)
  endif
endif
rdglschn=999
rewind lfn
do while(.true.)
  read(lfn,'(a)',end=100) line
  if(line(1:1).eq.'#') cycle
  if(line(1:3).eq.iprn) then
    read(line(5:),*) rdglschn
    exit
  endif
enddo
100 return
end
