!
!! ipp.f90
!!
!!    Copyright (C) 2021 by Wuhan University
!!
!!    This program belongs to PRIDE PPP-AR which is an open source software:
!!    you can redistribute it and/or modify it under the terms of the GNU
!!    General Public License (version 3) as published by the Free Software Foundation.
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
! Calculate ionosphere pierce point (IPP)
subroutine ipp(x,elev,azim,bradius,height,latipp,lonipp,zipp)
implicit none

real*8 x(1:*)          ! receiver coordinates in km
real*8 elev,azim       ! elevation and azimuth of satellite in radian
real*8 bradius,height  ! Earth radius and layer height in km
real*8 latipp,lonipp   ! IPP coordicantes in radian
real*8 zipp            ! zenith distance in radians at IPP
!
!! local
logical*1 lfirst
real*8 pi,alpha,lat,lon

data lfirst/.true./
save lfirst,pi
!
!! pi
if(lfirst) then
  lfirst=.false.
  pi=4.d0*datan(1.d0)
endif
!
!! geocentric cooridinates
lat=datan(x(3)/dsqrt(x(1)**2+x(2)**2))
lon=datan2(x(2),x(1))
!write(*,*)x(1),x(2),x(3)
!
!! zenith distance from satellite to IPP
zipp=dasin(bradius/(bradius+height)*dcos(elev))
!
!! geocentric angle between receiver and satellite
alpha=pi/2-elev-zipp
!
!! latitutde of IPP
if(1.d0 .lt. dsin(lat)*dcos(alpha)+dcos(lat)*dsin(alpha)*dcos(azim))then
  latipp=pi/2
elseif(-1.d0 .gt. dsin(lat)*dcos(alpha)+dcos(lat)*dsin(alpha)*dcos(azim))then
  latipp=-pi/2;
else
  latipp=dasin(dsin(lat)*dcos(alpha)+dcos(lat)*dsin(alpha)*dcos(azim))
endif
!
!! longitude of IPP
if(1.d0 .lt. dsin(alpha)*dsin(azim)/dcos(latipp))then
  lonipp=lon+pi/2
elseif(-1.d0 .gt. dsin(alpha)*dsin(azim)/dcos(latipp))then
  lonipp=lon+(-pi/2)
else
  lonipp=lon+dasin(dsin(alpha)*dsin(azim)/dcos(latipp))
endif
if(lonipp.gt.2.d0*pi) lonipp=lonipp-2.d0*pi

return
end
