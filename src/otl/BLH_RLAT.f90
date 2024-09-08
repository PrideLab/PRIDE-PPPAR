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
SUBROUTINE BLH_RLAT(GRS,BLH,RLAT)
      ! March 27, 2005, Zhang Chuanyin
       IMPLICIT REAL*8(A-H,O-Z)
       real*8 GRS(6),a,f,BLH(3),RLAT(3)
       real*8 N,K,flat,RLATB
!----------------------------------------------------------
       a=GRS(2); f=GRS(5) !1/f
       AE=a
       E2=(2.d0-f)*f
       pi=datan(1.d0)*4.d0
       RAD=datan(1.d0)/45.d0
       FLAT=BLH(1)
       FLON=BLH(2)
       HT=BLH(3)

       FLATR=FLAT*RAD
       FLONR=FLON*RAD
       T1=DSIN(FLATR)**2
       N=AE/DSQRT(1.d0-E2*T1)
       T2=(N+HT)*DCOS(FLATR)
       X=T2*DCOS(FLONR)
       Y=T2*DSIN(FLONR)
       Z=(N*(1.d0-E2)+HT)*DSIN(FLATR)
       N=AE/DSQRT(1.d0-E2*T1)
! COMPUTE THE GEOCENTRIC RADIUS
       RE=DSQRT(X**2+Y**2+Z**2)
! COMPUTE THE GEOCENTRIC LATITUDE
       RLATB=DATAN(Z/DSQRT(X**2+Y**2))
! COMPUTE NORMAL GRAVITY:UNITS ARE M/SEC**2
!      GR=GEQT*(1.+K*T1)/DSQRT(1.-E2*T1)
       RLAT(1)=RE
       RLAT(2)=RLATB/RAD
       RLAT(3)=FLON
       RETURN
END
