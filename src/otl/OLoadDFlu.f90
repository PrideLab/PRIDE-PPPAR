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
subroutine OLoadDFlu(mjd,cnm,snm,maxn,fes,nn)
  ! Direct impact of tidal loading on geoid coefficients calculated using tidal loading spherical harmonic model
  !fes(1:nn,1:7)-doodson,n,m,C+(cm),eps+,C-(cm),eps-
  ! cnm, snm return the direct impact of tidal loading on geoid coefficients, maxn is the maximum order of the model
  integer::nn,i,j,n,m,k,maxn,bs,td
  real*8::mjd,fes((maxn+2)**2*40,7),fk,rw,ge,gg,df,du
  real*8::t,pi,RAD,bias,fk1
  real*8::cnm((maxn+2)**2),snm((maxn+2)**2)
  real*8::sigma,th1,th2,astr(6),tt,dod(6),nm,cp,cn,sp,sn
!---------------------------------------------------------------------
  pi=datan(1.d0)*4.d0;RAD=pi/180.d0
  rw=1.025d3;re=5.517d3;rr=6378136.3d0
  ge=9.7803278d0;gg=6.67428d-11
  fk=4.d0*pi*gg*rw/ge;fk1=3.d0*rw/re/rr
  t=(mjd-51544.5d0)/36525.d0 !J2000.0:mjd51544.5d0
  tt=(mjd-51544.d0)-floor(mjd-51544.d0) ! Adjust for 0.5 days
      astr(2)=218.31664563d0+481267.88119575d0*t-0.001466388889d0*t**2 & !s
       -0.000000074112d0*t**3-0.000000153389d0*t**4
  astr(3)=280.4664501606d0+36000.769748805556d0*t+0.000303222222d0*t**2 & !h
       -0.000001905501d0*t**3-0.000000065361d0*t**4
  astr(4)=83.35324312d0+4069.01363525d0*t-0.010321722222d0*t**2 & !p
       -0.000014417168d0*t**3+0.0000000526333d0*t**4
  astr(5)=234.95544499d0+1934.136261972222d0*t-0.002075611111d0*t**2 & !N'
       -0.000000213944d0*t**3+0.000000164972d0*t**4
  astr(6)=282.9373409806d0+1.719457666668d0*t**2+0.000456888889d0*t**2 & !ps
       -0.000001943279d0*t**3-0.0000000033444d0*t**4
  astr(1)=360.d0*tt-astr(2)+astr(3) !tao
      astr=dmod(astr*RAD,2.d0*pi)
  cnm=0.d0;snm=0.d0
  do i=1,nn
    n=nint(fes(i,2));m=nint(fes(i,3))
    if(n<1)goto 9005 ! Three first-order terms
    k=n*(n+1)/2+m
    th1=0.d0;td=nint(fes(i,1));bs=100000
    do j=1,6
      dod(j)=td/bs;td=td-bs*dod(j);bs=bs/10
      if(j>1)dod(j)=dod(j)-5
    enddo
    do j=1,6
      th1=th1+dble(dod(j))*astr(j)
    enddo
    ! Calculate intersection correction parameters
    call CalcTidefu(mjd,fes(i,1),df,du)
    call BiasTide(fes(i,1),bias)   ! Extract phase bias correction
    th1=th1+(du+bias)*RAD
    th1=dmod(th1,2.d0*pi)
    !IERS2010(6.21)
    cp=fes(i,4)*dsin(fes(i,5)*RAD)/dble(2*n+1)
    sp=fes(i,4)*dcos(fes(i,5)*RAD)/dble(2*n+1)
    cn=fes(i,6)*dsin(fes(i,7)*RAD)/dble(2*n+1)
    sn=fes(i,6)*dcos(fes(i,7)*RAD)/dble(2*n+1)
    !IERS2010(6.15)
    cnm(k)=cnm(k)+((cp+cn)*dcos(th1)+(sp+sn)*dsin(th1))*df
    snm(k)=snm(k)-((cp-cn)*dsin(th1)-(sp-sn)*dcos(th1))*df
9005    continue
  enddo
  cnm=cnm*fk*1.d-2; snm=snm*fk*1.d-2
  return
end
