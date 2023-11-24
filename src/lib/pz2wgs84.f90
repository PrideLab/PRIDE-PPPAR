!
!! pz2wgs84.f90
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
!! Contributor: Jianghui Geng, Songfeng Yang
!! refer to RTKLIB
!! 
!!
!!
! transform PZ cooridnates to WGS84
subroutine pz2wgs84(jd,sod,xpz,xwgs)
implicit none

integer*4 jd
real*8    sod,xpz(1:*),xwgs(1:*)
!
!! local
integer*4 jd02,jd11
real*8    mas2rad,sod02,sod11
real*8    bias(3),srot(3,3)
!
!! function used
real*8    timdif

data mas2rad/4.84813681d-9/
data jd02,jd11,sod02,sod11/54363,56657,64814.d0,54016.d0/
save mas2rad,jd02,jd11,sod02,sod11

if(timdif(jd,sod,jd02,sod02).lt.0.d0) then
! PZ90
  bias(1)=0.07d0
  bias(2)=0.d0
  bias(3)=-0.77d0
  srot(1,1)=-3.d-9
  srot(2,1)=353.d0*mas2rad
  srot(3,1)=4.d0*mas2rad
  srot(1,2)=-353.d0*mas2rad
  srot(2,2)=-3.d-9
  srot(3,2)=-19.d0*mas2rad
  srot(1,3)=-4.d0*mas2rad
  srot(2,3)=19.d0*mas2rad
  srot(3,3)=-3.d-9
else if(timdif(jd,sod,jd11,sod11).lt.0.d0) then
! PZ90.02
  bias(1)=-0.36d0
  bias(2)=0.08d0
  bias(3)=0.18d0
  srot=0.d0
else
! PZ90.11
  bias(1)=0.003d0
  bias(2)=0.001d0
  bias(3)=0.001d0
  srot=0.d0
endif
xwgs(1:6)=0.d0
call matmpy(srot,xpz,xwgs,3,3,1)
xwgs(1:3)=xwgs(1:3)+bias(1:3)
call matmpy(srot,xpz(4),xwgs(4),3,3,1)
xwgs(1:6)=xpz(1:6)+xwgs(1:6)

return
end
