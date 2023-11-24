!
!! tide_displace.f90
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
!! Contributor: Maorong Ge, Jianghui Geng, Shuyin Mao, Songfeng Yang
!! 
!!
!!
!! purpose   : compute all tide related station position deformation (IERS2010)
!! parameters: 
!!        iers -- which terms are included
!!        jd, t  -- epoch time
!!        xsite, xsun xlun -- station-, solar- and lunar-positions
!!        rot_f2j, rot_l2f -- rotation matrix from earth-fixed to inertial
!!                and from station to earth-fixed system
!!        xpole,ypole -- sideral time, x and y pole positions
!!        olc  -- ocean loading coefficients
!!        dx  -- position correction
!
subroutine tide_displace(tide,jd,t,xsit_j,xsit_f,xsun,xlun,rot_f2j,rot_l2f,&
                         lat,lon,sidtm,xpole,ypole,olc,disp)
implicit none
include '../header/const.h'

character*(*) tide
integer*4 jd
real*8    xsit_j(1:*),xsit_f(1:*),xsun(1:*),xlun(1:*),disp(1:*)
real*8    t,lat,lon,xpole,ypole,sidtm,rot_f2j(3,3),rot_l2f(3,3),olc(11,6)
!
!! local
integer*4 i,j,jdutc,jdtdt,iy,imon,id,ih,im
real*8    tsec,fhr,xpm,ypm,tutc,ttdt,dxi(3),colat,xs(3),xl(3)
!
!! function called
real*8    dot,taiutc
!
!! initialization
do i=1,3
  disp(i)=0.d0
enddo
colat=2.d0*datan(1.d0)-lat
call timinc(jd,t,19.d0-taiutc(jd),jdutc,tutc)
call mjd2date(jdutc,tutc,iy,imon,id,ih,im,tsec)
fhr=ih+im/60.d0+tsec/3600.d0
call timinc(jd,t,GPSTDT,jdtdt,ttdt)
if(index(tide,'SOLID').ne.0) then
  !
  !! 1. Displacement due to frequency-independent solid-Earth tide(in J2000)
  call matmpy(xsun,rot_f2j,xs,1,3,3)
  call matmpy(xlun,rot_f2j,xl,1,3,3)
  call DEHANTTIDEINEL(xsit_f(1:3)*1.d3,iy,imon,id,fhr,xs*1.d3,xl*1.d3,dxi)
  call matmpy(rot_f2j,dxi,dxi,3,3,1)
  do i=1,3
    disp(i)=disp(i)+dxi(i)*1.d-3
  enddo
endif
  !
  !! 2. Displacement due to the pole tide 
  !! displacement in east, south and radial in mm
  !! xpole,ypole in seconds of arc
  !! dxi(3) east-north-radial
if(index(tide,'POLE').ne.0) then
  call mean_pole(jdtdt+ttdt/86400.d0,xpm,ypm)
  xpm=xpole-xpm
  ypm=-(ypole-ypm)
  dxi(1)=  9.d0*dcos(colat)     *(xpm*dsin(lon)-ypm*dcos(lon))
  dxi(2)=  9.d0*dcos(2.d0*colat)*(xpm*dcos(lon)+ypm*dsin(lon))  ! to north
  dxi(3)=-33.d0*dsin(2.d0*colat)*(xpm*dcos(lon)+ypm*dsin(lon))
  !
  !! rotation matrix from east-north-radial to x-y-z, then to J2000
  call matmpy(rot_l2f,dxi,dxi,3,3,1)  
  call matmpy(rot_f2j,dxi,dxi,3,3,1)
  do i=1,3
    disp(i)=disp(i)+dxi(i)*1.d-6
  enddo
endif
  !
  !! 3. Displacement due to ocean-loading
if(index(tide,'OCEAN').ne.0) then
  dxi(1:3)=0.d0
  call hardisp(jdutc,tutc,olc,1,0.d0,dxi)
  call matmpy(rot_l2f,dxi,dxi,3,3,1) 
  call matmpy(rot_f2j,dxi,dxi,3,3,1)
  do i=1,3
    disp(i)=disp(i)+dxi(i)*1.d-3
  enddo
endif

return
end
