!
!! istrue.f90
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
!! purpose  : check if a data flag is or belongs to the specified
!!            type.(GOOD,NODATA,OK etc)
!! parameter:
!!   j -- flag code
!!   c -- flag type description
!!   istrue -- return, if flag description matchs the code.
!!
!
logical*1 function istrue(i, c)
  include 'data_flag.h'

  integer*4 i
  character*(*) c

  istrue = .false.
!
!! data will be used in estimating
  if (c(1:2) .eq. 'ok' .or. c(1:2) .eq. 'OK') then
    if (i .eq. GOOD .or. i .eq. NEWAMB .or. i .eq. OLDAMB) istrue = .true.
!
!! data will be used in estimating but without amb. flag
  else if (c(1:4) .eq. 'good' .or. c(1:4) .eq. 'GOOD') then
    if (i .eq. GOOD) istrue = .true.

  else if (c(1:3) .eq. 'amb' .or. c(1:3) .eq. 'AMB') then
    if (i .eq. NEWAMB .or. i .eq. OLDAMB) istrue = .true.

  else if (c(1:6) .eq. 'nodata' .or. c(1:6) .eq. 'NODATA') then
    if (i .eq. NODATA) istrue = .true.

  else if (c(1:6) .eq. 'newamb' .or. c(1:6) .eq. 'NEWAMB') then
    if (i .eq. NEWAMB) istrue = .true.

  else if (c(1:6) .eq. 'oldamb' .or. c(1:6) .eq. 'OLDAMB') then
    if (i .eq. OLDAMB) istrue = .true.

  else if (c(1:6) .eq. 'delrms' .or. c(1:6) .eq. 'DELRMS') then
    if (i .eq. DELRMS) istrue = .true.

  else if (c(1:6) .eq. 'delbad' .or. c(1:6) .eq. 'DELBAD') then
    if (i .eq. DELBAD) istrue = .true.

  else if (c(1:6) .eq. 'delsht' .or. c(1:6) .eq. 'DELSHT') then
    if (i .eq. DELSHT) istrue = .true.

  else if (c(1:6) .eq. 'newbad' .or. c(1:6) .eq. 'NEWBAD') then
    if (i .eq. NEWBAD) istrue = .true.
  else
    write (*, '(2a)') '***ERROR(istrue): unknow ', trim(c)
    call exit(1)
  endif

  return
end
