!
!! mjd2date.f90
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
!! Contributor: Maorong Ge, Jianghui Geng, Songfeng Yang
!! 
!!
!!
!! purpose  : transform modified julday to date
!! parameter:
!!    input : jd,sod -- modified julday
!!    output: iy,imon,id,ih,im,is -- date
!
subroutine mjd2date(jd, sod, iy, imon, id, ih, im, is)
  implicit none

  integer*4 jd, iy, imon, id, ih, im
  real*8 sod, is
!
!! local
  integer*4 mjd, doy
  real*8 msod
!
!! check type of modified julday
  if (jd .ne. 0) then
    mjd = jd
    msod = sod
  else
    mjd = int(sod)
    msod = (sod - mjd)*86400.d0
  endif
!
!! transformation
  call mjd2doy(mjd, iy, doy)
  call yeardoy2monthday(iy, doy, imon, id)
  call sod2hms(msod, ih, im, is)

  return
end
