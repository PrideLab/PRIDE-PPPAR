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
subroutine GNormalfd(BLH,NFD,GRS)
      ! NFD(5) - normal gravity potential, normal gravity, normal gravity gradient,
      ! normal gravity direction, normal gradient direction
      implicit none
      real*8::BLH(3),NFD(5),val(5),rln(3),pn(40),dp1(40),dp2(40)
      integer::i,j,n
      real*8::GRS(6),djn(80),pi,RAD,gm,ae,ww,rr,t,u,vdn(5),atr
!---------------------------------------------------------------------
      pi=datan(1.d0)*4.d0;RAD=pi/180.d0
      call normdjn(GRS,djn)
      gm=GRS(1);ae=GRS(2);ww=GRS(4);vdn=0.d0
      call BLH_RLAT(GRS,BLH,rln)
      rr=rln(1);t=dsin(rln(2)*RAD);u=dcos(rln(2)*RAD)
      call LegPn_dt2(pn,dp1,dp2,40,t)
      do n=1,20
        atr=dexp(dble(2*n)*dlog(ae/rr))*djn(2*n)
        vdn(1)=vdn(1)+atr*pn(2*n)
        vdn(2)=vdn(2)+atr*pn(2*n)*dble(2*n+1)
        vdn(3)=vdn(3)+atr*dp1(2*n)
        vdn(4)=vdn(4)+atr*pn(2*n)*dble((2*n+1)*(n+1))
        vdn(5)=vdn(5)+atr*dp2(2*n)
      enddo
      val(1)=gm/rr*(1.d0-vdn(1))+(ww*rr*u)**2/2.d0
      val(2)=-gm/rr**2*(1.d0-vdn(2))+ww**2*rr*u**2 ! Radial normal gravity
      val(3)=-gm/rr**2*vdn(3)+ww**2*rr*u*t         ! Southward normal gravity
      val(4)=-2.d0*gm/rr**3*(1.d0-vdn(4))+(ww*u)**2
      val(5)=-gm/rr**3*vdn(5)+ww**2*(2.d0*t**2-1.d0)
      NFD(4)=datan(val(3)/val(2))/RAD*36.d2
      NFD(5)=datan(val(5)/val(4))/RAD*36.d2
      val(2)=dsqrt(val(2)**2+val(3)**2)
      val(4)=dsqrt(val(4)**2+val(5)**2)
      NFD(1)=val(1)
      NFD(2)=val(2)
      NFD(3)=val(4)
end
!
!*********************************************************
!
subroutine normdjn(GRS,djn)
! Input GRS - gm, ae, j2, omega, 1/f
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c     djn  array of even zonals of normal ellipsoid
!c     ex. djn(2)=j2=GRS(3)
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit real*8 (a-h,o-z)
      real*8 GRS(6),djn(80)
  ae=GRS(2)
  ff=GRS(5)
  do 101 n2=1,80,1
  djn(n2)=0.d0
  101  continue
  djn(2)=GRS(3)
      esq=(2.d0-ff)*ff!e**2
!c     compute normal even zonals from 4 to 20
      do 200 n2=4,80,2
        n=n2/2
        djn(n2)=(-1.d0)**(n+1)*3.d0*esq**n/(n2+1.d0)/(n2+3.d0) &
               *(1.d0-dble(n)+5.d0*dble(n)*djn(2)/esq)
  200 continue
      return
end
