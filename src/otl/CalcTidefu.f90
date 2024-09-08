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
      subroutine BiasTide(doodson,bias)
      real*8::doodson,bias
      integer i
      !IERS2010 Table6.6, Table6.7
!     -------------------------------------
      bias=0.d0 !(n-m)*90
      if(dabs(doodson- 55565)<1.d-3)bias=180   ! omiga1
      if(dabs(doodson-125755)<1.d-3)bias=-90   ! 2Q1
      if(dabs(doodson-127555)<1.d-3)bias=-90   ! sigma1
      if(dabs(doodson-135655)<1.d-3)bias=-90   ! Q1
      if(dabs(doodson-137455)<1.d-3)bias=-90   ! rou1
      if(dabs(doodson-145545)<1.d-3)bias=-90   ! rou1
      if(dabs(doodson-145555)<1.d-3)bias=-90   ! O1
      if(dabs(doodson-155555)<1.d-3)bias=-90   ! M1
      if(dabs(doodson-155655)<1.d-3)bias=90    ! M1
      if(dabs(doodson-157455)<1.d-3)bias=90    ! kama1
      if(dabs(doodson-162556)<1.d-3)bias=-90   ! pi1
      if(dabs(doodson-163555)<1.d-3)bias=-90   ! P1
      if(dabs(doodson-163655)<1.d-3)bias=90    ! theta1
      if(dabs(doodson-164556)<1.d-3)bias=90    ! S1
      if(dabs(doodson-165555)<1.d-3)bias=90    ! K1
      if(dabs(doodson-166554)<1.d-3)bias=90    ! psi1
      if(dabs(doodson-167555)<1.d-3)bias=90    ! phi1
      if(dabs(doodson-175455)<1.d-3)bias=90    ! J1
      if(dabs(doodson-185555)<1.d-3)bias=90    ! OO1
      if(dabs(doodson-253755)<1.d-3)bias=180   ! gama2
      if(dabs(doodson-255545)<1.d-3)bias=180   ! alfa2
      if(dabs(doodson-265455)<1.d-3)bias=180   ! L2
      if(dabs(doodson-263655)<1.d-3)bias=180   ! namta2
      if(dabs(doodson-275554)<1.d-3)bias=180   ! R2
      if(dabs(doodson-355555)<1.d-3)bias=180   ! M3
      end
!
!*********************************************************
!
      subroutine CalcTidefu(mjd,doodson,df,du)
      real*8::mjd,doodson,df,du,pi,RAD
      real*8::shpn(5),td,f(22),u(22),ddsn(22),SINN,SIN2N,COSN,COS2N,OMEGA
      integer i
!     -------------------------------------
      pi=datan(1.d0)*4.d0;RAD=pi/180.d0;df=1.d0;du=0.d0
         CALL ASTRO5(mjd, shpn )
         OMEGA=shpn(4)
         ddsn(1)=135655! Q1
         ddsn(2)=145555! O1
         ddsn(3)=163555! P1
         ddsn(4)=165555! K1
         ddsn(5)=245655! N2
         ddsn(6)=255555! M2
         ddsn(7)=273555! S2
         ddsn(8)=275555! K2
         ddsn(9)=125755! 2Q1
         ddsn(10)=127555! sigma1
         ddsn(11)=137455! rho1
         ddsn(12)=155655! M1
         ddsn(13)=175455! J1
         ddsn(14)=185555! OO1
         ddsn(15)=235755! 2N2
         ddsn(16)=237555! mu2
         ddsn(17)=247455! nu2
         ddsn(18)=265455! L2
         ddsn(19)=164556! S1 (Doodson's phase)
         ddsn(20)=355555! M3
         ddsn(21)=455555! M4
         ddsn(22)=655555! M6
!
!     determine nodal corrections f and u    
!        Note: Update this code next iteration of model  -RDR
!     -----------------------------------
         SINN = DSIN(OMEGA*RAD)
         COSN = DCOS(OMEGA*RAD)
         SIN2N = DSIN(2.d0*OMEGA*RAD)
         COS2N = DCOS(2.d0*OMEGA*RAD)

         F(1) = 1.009 + 0.187*COSN - 0.015*COS2N
         F(2) = F(1)
         F(4) = 1.006 + 0.115*COSN - 0.009*COS2N
         F(3) = F(4)
         F(5) = 1.000 - 0.037*COSN
         F(6) = F(5)
         F(8) = 1.024 + 0.286*COSN + 0.008*COS2N
         F(7) = F(8)
         F(9) = DSQRT((1.0 + 0.189*COSN - 0.0058*COS2N)**2 + &
                     (0.189*SINN - 0.0058*SIN2N)**2)
         F(10) = F(9)
         F(11) = F(9)
         F(12) = DSQRT((1.0 + 0.185*COSN)**2 + (0.185*SINN)**2)
         F(13) = DSQRT((1.0 + 0.198*COSN)**2 + (0.198*SINN)**2)
         F(14) = DSQRT((1.0 + 0.640*COSN + 0.134*COS2N)**2 + &
                     (0.640*SINN + 0.134*SIN2N)**2 )
         F(15) = DSQRT((1.0 - 0.0373*COSN)**2 + (0.0373*SINN)**2)
         F(16) = F(15)
         F(17) = F(15)
         F(18) = F(15)
         F(19) = F(3)
         F(20) = F(6)*dsqrt(F(6))
         F(21) = F(6)*F(6)
         F(22) = F(6)*F(6)*F(6)

         U(1) = 10.8*SINN - 1.3*SIN2N
         U(2) = U(1)
         U(4) = -8.9*SINN + 0.7*SIN2N
         U(3) = U(4)
         U(5) = -2.1*SINN
         U(6) = U(5)
         U(8) = -17.7*SINN + 0.7*SIN2N
         U(7) = U(8)
         U(9) = DATAN2(0.189*SINN - 0.0058*SIN2N, &
                     1.0 + 0.189*COSN - 0.0058*SIN2N)/RAD
         U(10) = U(9)
         U(11) = U(9)
         U(12) = DATAN2( 0.185*SINN, 1.0 + 0.185*COSN)/RAD
         U(13) = DATAN2(-0.198*SINN, 1.0 + 0.198*COSN)/RAD
         U(14) = DATAN2(-0.640*SINN - 0.134*SIN2N, &
                      1.0 + 0.640*COSN + 0.134*COS2N)/RAD
         U(15) = DATAN2(-0.0373*SINN, 1.0 - 0.0373*COSN)/RAD
         U(16) = U(15)
         U(17) = U(15)
         U(18) = U(15)
         U(19) = U(3)
         U(20) = 1.5d0*U(6)
         U(21) = 2.d0*U(6)
         U(22) = 3.d0*U(6)
         do i=1,22
             if(dabs(doodson-ddsn(i))<0.1)then
                 df=f(i);du=u(i);return
             endif
         enddo
      return
      end
!
!*********************************************************************************
!

      SUBROUTINE ASTRO5( TIME, SHPNP )
!
!---------------------------------------------------------------------
!  Computes the 5 basic astronomical mean longitudes  s, h, p, N, p'.
!
!  Note N is not N', i.e. N is decreasing with time.
!
!  TIME is UTC in decimal Modified Julian Day (MJD).
!  All longitudes returned in degrees.
!
!  R. D. Ray, NASA/GSFC   August 2003
!
!  Most of the formulae for mean longitudes are extracted from 
!  Jean Meeus, Astronomical Algorithms, 2nd ed., 1998.  
!  Page numbers below refer to this book.
!
!  Note: This routine uses TIME in UT and does not distinguish between
!    the subtle differences of UTC, UT1, etc.  This is more than adequate
!    for the calculation of these arguments, especially in tidal studies.
!---------------------------------------------------------------------
!
      IMPLICIT NONE
      DOUBLE PRECISION TIME, SHPNP(5)

      DOUBLE PRECISION TJD,T,CIRCLE,DEL,TJLAST,D
      PARAMETER       (CIRCLE=360.0D0)

!     Convert to Julian Day and to Ephemeris Time
!     -------------------------------------------
      TJD = TIME + 2400000.5D0

!     Compute time argument in centuries relative to J2000 
!     ----------------------------------------------------
      T = ( TJD - 2451545.d0 )/36525.d0
!
!
!     mean longitude of moon (p.338)
!     ------------------------------
      SHPNP(1) = (((-1.53388d-8*T + 1.855835d-6)*T - 1.5786d-3)*T + &
                 481267.88123421d0)*T + 218.3164477d0
!
!     mean elongation of moon (p.338)
!     -------------------------------
      D = (((-8.8445d-9*T + 1.83195d-6)*T - 1.8819d-3)*T + &
               445267.1114034d0)*T + 297.8501921d0
!
!     mean longitude of sun
!     ---------------------
      SHPNP(2) = SHPNP(1) - D
!
!     mean longitude of lunar perigee (p.343)
!     ---------------------------------------
      SHPNP(3) = ((-1.249172d-5*T - 1.032d-2)*T + 4069.0137287d0)*T + &
                 83.3532465d0
!
!     mean longitude of ascending lunar node (p.144)
!     ----------------------------------------------
      SHPNP(4) = ((2.22222d-6*T + 2.0708d-3)*T - 1934.136261d0)*T + &
                 125.04452d0
!
!     mean longitude of solar perigee (Simon et al., 1994)
!     ----------------------------------------------------
      SHPNP(5) = 282.94d0 + 1.7192d0 * T

      SHPNP(1) = DMOD(SHPNP(1),CIRCLE)
      SHPNP(2) = DMOD(SHPNP(2),CIRCLE)
      SHPNP(3) = DMOD(SHPNP(3),CIRCLE)
      SHPNP(4) = DMOD(SHPNP(4),CIRCLE)

      IF (SHPNP(1).LT.0.D0) SHPNP(1) = SHPNP(1) + CIRCLE
      IF (SHPNP(2).LT.0.D0) SHPNP(2) = SHPNP(2) + CIRCLE
      IF (SHPNP(3).LT.0.D0) SHPNP(3) = SHPNP(3) + CIRCLE
      IF (SHPNP(4).LT.0.D0) SHPNP(4) = SHPNP(4) + CIRCLE
      RETURN
      END
