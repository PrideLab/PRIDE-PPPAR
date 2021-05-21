!
!! interp_ionex.f90
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
!! Contributor: Jianghui Geng
!! 
!!
!!
!! ionex interpolation
subroutine interp_ionex(IM,ind,lat,lonost,iondelay,ionrms)
implicit none
include '../header/const.h'
include '../header/ionex.h'

type(ionex) IM
integer*4   ind
real*8      lat,lonost,iondelay,ionrms ! in TECU
!
!! local
logical*1   lfirst
integer*4   i,j,ie,je
real*8      rad2deg,lon,dump1,dump2

data lfirst/.true./
save lfirst,rad2deg

if(lfirst) then
  lfirst=.false.
  rad2deg=180.d0/pi
endif
iondelay=0.d0
ionrms=0.d0
! transform longitude to [-pi, pi]
if(lonost.gt.2.d0*pi) then
  lon=lonost-2.d0*pi
else if(lonost.lt.0.d0) then
  lon=lonost+2.d0*pi
else
  lon=lonost
endif
if(lon.gt.pi) lon=lon-2.d0*pi
!
!! locate grid
if(lat*rad2deg.ge.IM%lat1) then
  i =1
  ie=1
else if(lat*rad2deg.le.IM%lat2) then
  i =IM%blat
  ie=IM%blat
else
  i =int((lat*rad2deg-IM%lat1)/IM%dlat)+1
  ie=i+1
endif
if(lon*rad2deg.le.IM%lon1) then
  j =1
  je=1
else if(lon*rad2deg.ge.IM%lon2) then
  j =IM%blon
  je=IM%blon
else
  j =int((lon*rad2deg-IM%lon1)/IM%dlon)+1
  je=j+1
endif
!
!! bi-linear interpolation (in TECU)
if(IM%tecmp(j, i,ind).eq.9999.or.IM%tecmp(je, i,ind).eq.9999.or.&
   IM%tecmp(j,ie,ind).eq.9999.or.IM%tecmp(je,ie,ind).eq.9999) return
dump1=IM%tecmp(j, i,ind)+(IM%tecmp(je, i,ind)-IM%tecmp(j, i,ind))*(lon*rad2deg-IM%lon1-(j-1)*IM%dlon)/IM%dlon
dump2=IM%tecmp(j,ie,ind)+(IM%tecmp(je,ie,ind)-IM%tecmp(j,ie,ind))*(lon*rad2deg-IM%lon1-(j-1)*IM%dlon)/IM%dlon
iondelay=dump1+(dump2-dump1)*(lat*rad2deg-IM%lat1-(i-1)*IM%dlat)/IM%dlat

if(IM%rmsmp(j, i,ind).eq.9999.or.IM%rmsmp(je, i,ind).eq.9999.or.&
   IM%rmsmp(j,ie,ind).eq.9999.or.IM%rmsmp(je,ie,ind).eq.9999) return
dump1=IM%rmsmp(j, i,ind)+(IM%rmsmp(je, i,ind)-IM%rmsmp(j, i,ind))*(lon*rad2deg-IM%lon1-(j-1)*IM%dlon)/IM%dlon
dump2=IM%rmsmp(j,ie,ind)+(IM%rmsmp(je,ie,ind)-IM%rmsmp(j,ie,ind))*(lon*rad2deg-IM%lon1-(j-1)*IM%dlon)/IM%dlon
ionrms=dump1+(dump2-dump1)*(lat*rad2deg-IM%lat1-(i-1)*IM%dlat)/IM%dlat

return
end
