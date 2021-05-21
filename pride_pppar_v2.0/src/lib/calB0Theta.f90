!
!! calB0Theta.f90
!!
!!    Copyright (C) 2021 by Wuhan University
!!
!!    This program is an open source software: you can redistribute it and/or
!modify
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
!! Contributor: Shaoming Xin
!! 
!!
!!
!
subroutine calB0Theta(rf,DATE,ALT,XLT,XLN,b0,theta) 
  implicit none
  real*8 DATE,ALT,XLT,XLN,b0,theta !pp's lat lon
  real*8 D,S,IH,IX,IY,IZ,NF
  real*8 Drad,Srad,aaa
  real*8 pi
  real*8 ppenu(3),dump0(3,1),dumpm(3,1),dump(3)
  real*8 rf(3),rfu(3) ! vector in fixed system
  real*8 rot_l2f(3,3),d0d(3,3)
  pi=4.d0*datan(1.d0)


  call IGRF12(DATE,ALT,XLT,XLN,D,S,IH,IX,IY,IZ,NF)
  b0=dabs(NF) !unit nt

  !---cal theta
  Drad=D/180.d0*pi
  Srad=S/180.d0*pi

  dump=0.d0
  ppenu(1)=IY  !1.d0*dcos(Drad)
  ppenu(2)=IX  !*ALT !1.d0*dsin(Drad)
  ppenu(3)=-IZ !*(-ALT*dsin(pi/2.d0-XLT/180.d0*pi)) !1.d0*dtan(Srad)
  !
  call rot_enu2xyz(XLT/180.d0*pi,XLN/180.d0*pi,rot_l2f)
  call matmpy(rot_l2f,ppenu,dump,3,3,1)  !dump mag xyz vecto


  !
  ! sat- > sit
  rfu(1:3)=rf(1:3)/dsqrt(rf(1)**2+rf(2)**2+rf(3)**2)
 
  theta=dacos((rfu(1)*dump(1)+rfu(2)*dump(2)+rfu(3)*dump(3))/&
        (dsqrt(rfu(1)**2+rfu(2)**2+rfu(3)**2)*dsqrt(dump(1)**2+dump(2)**2+dump(3)**2)))
  return
end 

!subroutine calB0Theta4(zipp,xs,DATE,ALT,XLT,XLN,b0,theta)
!  implicit none
!  real*8 zipp,DATE,ALT,XLT,XLN,b0,theta !pp's lat lon
!  real*8 D,S,IH,IX,IY,IZ,NF
!  real*8 Drad,Srad
!  real*8 pi
!  real*8 rfu(3),dump(3),ixyz(3),xs(3),azi
!  
!  pi=4.d0*datan(1.d0)
!
!  call IGRF12(DATE,ALT,XLT,XLN,D,S,IH,IX,IY,IZ,NF)
!  b0=dabs(NF) !unit nt
!
!  !ppenu
!  dump(1)=IY 
!  dump(2)=IX 
!  dump(3)=-IZ
!  dump(1:3)=dump(1:3)/dsqrt(dump(1)**2+dump(2)**2+dump(3)**2)
!
!  ixyz(1)=ALT*dsin(XLT/180.d0*pi)*dsin(XLN/180.d0*pi)*1.d3
!  ixyz(2)=ALT*dsin(XLT/180.d0*pi)*dcos(XLN/180.d0*pi)*1.d3
!  ixyz(3)=ALT*dcos(XLT/180.d0*pi)*1.d3
!  azi=datan2((xs(1)-ixyz(1)),(xs(2)-ixyz(2)))
!  rfu(2)=-dsin(zipp)*dcos(azi)
!  rfu(1)=-dsin(zipp)*dsin(azi)
!  rfu(3)=dcos(zipp)
!  !rfu=-rfu
!
!  theta=dacos((rfu(1)*dump(1)+rfu(2)*dump(2)+rfu(3)*dump(3))/&
!        (dsqrt(rfu(1)**2+rfu(2)**2+rfu(3)**2)*dsqrt(dump(1)**2+dump(2)**2+dump(3)**2)))
!  return  
!end

!subroutine calB0Theta5(xr,elev,azim,DATE,ALT,XLT,XLN,b0,theta)
!  implicit none
!  real*8 lat,lon,elev,azim,DATE,ALT,XLT,XLN,b0,theta !pp's lat lon
!  real*8 D,S,IH,IX,IY,IZ,NF
!  real*8 Drad,Srad
!  real*8 pi,ixyz(3)
!  real*8 rfu(3),dump0(3),dump(3),rot(3,3),rea(3),xr(3),geod(3)
!  
!  pi=4.d0*datan(1.d0)
!
!  call IGRF12(DATE,ALT,XLT,XLN,D,S,IH,IX,IY,IZ,NF)
!  b0=dabs(NF) !unit nt
!
!  !ppenu
!  dump0(1)=IY 
!  dump0(2)=IX 
!  dump0(3)=-IZ
!  !dump(1:3)=dump(1:3)/dsqrt(dump(1)**2+dump(2)**2+dump(3)**2)
!
!  !! geocentric cooridinates
!  lat=datan(xr(3)/dsqrt(xr(1)**2+xr(2)**2))
!  lon=datan2(xr(2),xr(1))
!  !call xyzblh(xr(1:3)*1.d3,1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,geod)
!  !lat=geod(1)
!  !lon=geod(2)
!  !call rot_enu2xyz(lat,lon,rot)
!  rot(1,1)=dsin(lon);  rot(1,2)=dsin(lat)*dcos(lon); rot(1,3)=-dcos(lat)*dcos(lon)
!  rot(2,1)=-dcos(lon); rot(2,2)=dsin(lon)*dsin(lat); rot(2,3)=-dsin(lon)*dcos(lat)
!  rot(3,1)=0.d0;       rot(3,2)=-dcos(lat);          rot(3,3)=-dsin(lat)
!  rea(1)=dcos(elev)*dsin(azim)
!  rea(2)=dcos(elev)*dcos(azim)
!  rea(3)=dsin(elev)
!  !rea=-rea 
!  call matmpy(rot,rea,rfu,3,3,1)
!
!  lat=XLT/180.d0*pi
!  lon=XLN/180.d0*pi
!  ixyz(1)=ALT*dsin(lat)*dsin(lon)
!  ixyz(2)=ALT*dsin(lat)*dcos(lon)
!  ixyz(3)=ALT*dcos(lat)
!  call xyzblh(ixyz(1:3)*1.d3,1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,geod)
!  lat=geod(1)
!  lon=geod(2)
!  !call rot_enu2xyz(lat,lon,rot)
!  rot(1,1)=dsin(lon);  rot(1,2)=dsin(lat)*dcos(lon); rot(1,3)=-dcos(lat)*dcos(lon)
!  rot(2,1)=-dcos(lon); rot(2,2)=dsin(lon)*dsin(lat); rot(2,3)=-dsin(lon)*dcos(lat)
!  rot(3,1)=0.d0;       rot(3,2)=-dcos(lat);          rot(3,3)=-dsin(lat)
!  rot=-rot
!  call matmpy(rot,dump0,dump,3,3,1)
!
!  theta=dacos((rfu(1)*dump(1)+rfu(2)*dump(2)+rfu(3)*dump(3))/&
!        (dsqrt(rfu(1)**2+rfu(2)**2+rfu(3)**2)*dsqrt(dump(1)**2+dump(2)**2+dump(3)**2)))
!  return
!end




subroutine calB0Theta3(xr,rf,DATE,ALT,XLT,XLN,b0,theta) 
  implicit none
  real*8 DATE,ALT,XLT,XLN,b0,theta !pp's lat lon
  real*8 D,S,IH,IX,IY,IZ,NF
  real*8 Drad,Srad
  real*8 pi
  real*8 ppenu(3),dump(3)
  real*8 rf(3),rfu(3) ! vector in fixed system
  real*8 xr(3),et(3),nt(3),ut(3)
  pi=4.d0*datan(1.d0)


  call IGRF12(DATE,ALT,XLT,XLN,D,S,IH,IX,IY,IZ,NF)
  b0=dabs(NF) !unit nt

  !ppenu
  ppenu(1)=IY 
  ppenu(2)=IX 
  ppenu(3)=IZ 
  !
  call local_terrestrial_vector(xr,et,nt,ut)
  dump(1)=ppenu(2)*nt(1)+ppenu(1)*et(1)+ppenu(3)*ut(1)
  dump(2)=ppenu(2)*nt(2)+ppenu(1)*et(2)+ppenu(3)*ut(2)
  dump(3)=ppenu(2)*nt(3)+ppenu(1)*et(3)+ppenu(3)*ut(3)
  !
  rfu(1:3)=rf(1:3)/dsqrt(rf(1)**2+rf(2)**2+rf(3)**2)
  theta=dacos((rfu(1)*dump(1)+rfu(2)*dump(2)+rfu(3)*dump(3))/&
       (dsqrt(rfu(1)**2+rfu(2)**2+rfu(3)**2)*dsqrt(dump(1)**2+dump(2)**2+dump(3)**2)))
  return
end
!
subroutine local_terrestrial_vector(xr,et,nt,ut)
  implicit none
  real*8 xr(3),et(3),nt(3),ut(3)
  real*8 xnorm(3),norm,zaxis(3)

  zaxis(1)=0.d0
  zaxis(2)=0.d0
  zaxis(3)=1.d0

  norm=dsqrt(xr(1)**2+xr(2)**2+xr(3)**2)
  xnorm(1)=xr(1)/norm
  xnorm(2)=xr(2)/norm
  xnorm(3)=xr(3)/norm

  ut(1)=xnorm(1)*(-1.d0)
  ut(2)=xnorm(2)*(-1.d0)
  ut(3)=xnorm(3)*(-1.d0)
  
  et(1)=ut(2)*zaxis(3)-ut(3)*zaxis(2)
  et(2)=ut(3)*zaxis(1)-ut(1)*zaxis(3)
  et(3)=ut(1)*zaxis(2)-ut(2)*zaxis(1)
  norm=dsqrt(et(1)**2+et(2)**2+et(3)**2)
  et(1)=et(1)/norm
  et(2)=et(2)/norm
  et(3)=et(3)/norm

  nt(1)=et(2)*ut(3)-et(3)*ut(2)
  nt(2)=et(3)*ut(1)-et(1)*ut(3)
  nt(3)=et(1)*ut(2)-et(2)*ut(1)
  norm=dsqrt(nt(1)**2+nt(2)**2+nt(3)**2)
  nt(1)=nt(1)/norm
  nt(2)=nt(2)/norm
  nt(3)=nt(3)/norm

end
