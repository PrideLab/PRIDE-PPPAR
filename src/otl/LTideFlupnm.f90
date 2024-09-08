!******************************************************************************************
!     Copyright (C) 2024 ZHANG Chuanyin
!
!     This module was written by ZHANG Chuanyin and is used in the PRIDE PPP-AR project.
!     All rights reserved.
!
!     This module is protected by copyright law. Unauthorized copying, modification, publication,
!     transmission, display, performance, sale, distribution of this module's code, in whole or
!     in part, is strictly prohibited without the prior written permission of the original author.
!     You may use this module code in accordance with the license agreement of the original author
!     and PRIDE Lab.
!
! Website: www.zcyphygeodesy.com
! Contact Information: zhangchy@casm.ac.cn
!******************************************************************************************
subroutine LTideFlupnm(rln,maxn,cnm,snm,flv,tdn,GRS,pnm,dpt1,dpt2,gr)
  ! Calculation of deformation field due to surface equivalent water height using spherical harmonic synthesis method
  ! Input associated Legendre functions and their derivatives, normal gravity
  ! tdn(14) - height anomaly, surface gravity, gravity disturbance, ground tilt, deflection of the vertical, 
  ! northeast horizontal, radial, normal height, gravity disturbance gradient (radial mE), 
  ! north-west horizontal gradient (mE)
  implicit none
  integer::maxn,nn,n,m,kk,i,j
  real*8::cnm((maxn+2)**2),snm((maxn+2)**2),pnm((maxn+2)**2),dpt1((maxn+2)**2),dpt2((maxn+2)**2)
  real*8::gm,tdn(3),tn(3),gr,GRS(6),cosml,sinml,kp
  real*8::BLH(3),rln(3),NFD(5),t,u
  real*8::rr,rlat,rlon,ae,pi,RAD,flv(4000,3),dk,dh,dl
!---------------------------------------------------------------------------
  kp=1.d0
  tdn=0.d0;gm=GRS(1);ae=GRS(2);pi=datan(1.d0)*4.d0;RAD=pi/180.d0
  rr=rln(1);rlat=rln(2);rlon=rln(3)
  t=dsin(rlat*RAD);u=dcos(rlat*RAD)
    do n=1,maxn  
      dh=flv(n,1);dl=flv(n,2);dk=flv(n,3);tn=0.d0
      do m=0,n
        j=n*(n+1)/2+m
        cosml=dcos(rlon*RAD*dble(m));sinml=dsin(rlon*RAD*dble(m))
        tn(1)=tn(1)+dble(m)*(cnm(j)*sinml-snm(j)*cosml)*pnm(j+1)*dl
        tn(2)=tn(2)+(cnm(j)*cosml+snm(j)*sinml)*dpt1(j+1)*dl
        tn(3)=tn(3)+(cnm(j)*cosml+snm(j)*sinml)*pnm(j+1)*dh
      enddo
    tn=tn*dexp(dble(n)*dlog(ae/rr))
    tdn(1)=tdn(1)+tn(1)
    tdn(2)=tdn(2)+tn(2)
    tdn(3)=tdn(3)+tn(3)
  enddo
  tdn(1)=-tdn(1)*gm/rr/gr/u
  tdn(2)=-tdn(2)*gm/rr/gr
  tdn(3)=tdn(3)*gm/rr/gr
  return
end
