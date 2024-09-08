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
subroutine PlmBar_d(p,dp,lmax,rlat)
  implicit none
  integer::lmax,m,n,kk
  real*8 ::p((lmax+2)*(lmax+3)/2+5),dp((lmax+2)*(lmax+3)/2+5)
  real*8 ::rleg(lmax+5),dleg(lmax+5),rlat
!---------------------------------------------------------------------------
! Calculate all Legendre functions pnm; the indices of array pnm are 3 greater than those of cnm and snm
  do n=1,lmax+1
     call Legendre(lmax,n-1,rlat,rleg,dleg)
     do m=n,lmax+1
        kk=m*(m-1)/2+n
        p(kk)=rleg(m)
        dp(kk)=dleg(m)
     enddo
  enddo
end
!
!******************************************************************************
!
      subroutine Legendre(maxn,m,rlat,RLEG,DLEG)
!        
!        THIS SUBROUTINE COMPUTES  ALL NORMALIZED LEGENDRE FUNCTION
!        IN "RLEG" AND THEIR DERIVATIVES IN "DLEG". ORDER IS ALWAYS
!        M , AND COLATITUDE IS ALWAYS THETA  (RADIANS). MAXIMUM DEG
!        IS  NMX  . ALL CALCULATIONS IN DOUBLE PRECISION.
!        IR  MUST BE SET TO ZERO BEFORE THE FIRST CALL TO THIS SUB.
!        THE DIMENSIONS OF ARRAYS  RLEG, DLEG, AND RLNN  MUST BE
!        AT LEAST EQUAL TO  NMX+1  .
!
!        THIS PROGRAM DOES NOT COMPUTE DERIVATIVES AT THE POLES .
!
!        IF    IFLAG = 1  , ONLY THE LEGENDRE FUNCTIONS ARE
!        COMPUTED.
!
!    ORIGINAL PROGRAMMER :OSCAR L. COLOMBO, DEPT. OF GEODETIC SCIENCE
!    THE OHIO STATE UNIVERSITY, AUGUST 1980 . ******************
!
!     LEGFDN(M,RLAT,RLEG,DLEG,NMX,IR,RLNN,IFLAG)
      IMPLICIT REAL*8 (A-H,O-Z)
      integer::iflag,maxn,m,nmx
      real*8::RLEG(maxn+5),DLEG(maxn+5),RLNN(2*maxn+5), DRTS(2*maxn+5),DIRT(2*maxn+5)
      real*8::rlat,pi,dtr

      iflag=2
      pi=datan(1.d0)*4.d0
      dtr=pi/180.d0
      NMX=maxn
      THETA=pi/2.d0-(rlat*dtr)
      IR=0

      NMX1 = NMX+1
      NMX2P = 2*NMX+1
      M1 = M+1
      M2 = M+2
      M3 = M+3
      IF(IR.EQ.1) GO TO 10
      IR = 1
      DO 5     N = 1,NMX2P
      DRTS(N) = DSQRT(N*1.D0)
    5 DIRT(N) = 1.D0/DRTS(N)
   10 COTHET = DCOS(THETA)
      SITHET = DSIN(THETA)
      IF(IFLAG.NE.1.AND.THETA.NE.0.D0)SITHI = 1.D0/SITHET
!
!            COMPUTE THE LEGENDRE FUNCTIONS .
!
      RLNN(1) = 1.D0
      RLNN(2) = SITHET*DRTS(3)
      DO 15    N1 = 3,M1
      N = N1-1
      N2 = 2*N
   15 RLNN(N1) = DRTS(N2+1)*DIRT(N2)*SITHET*RLNN(N1-1)
      IF(M.GT.1) GO TO 20
      IF(M.EQ.0) GO TO 16
      RLEG(2) = RLNN(2)
      RLEG(3) = DRTS(5)*COTHET*RLEG(2)
      GO TO 20
   16 RLEG(1) = 1.D0
      RLEG(2) = COTHET*DRTS(3)
   20 CONTINUE
      RLEG(M1) = RLNN(M1)
        IF(M2.GT.NMX1) GO TO 35
      RLEG(M2) = DRTS(M1*2+1)*COTHET*RLEG(M1)
        IF(M3.GT.NMX1) GO TO 35
      DO 30     N1 = M3,NMX1
      N = N1-1
      IF(M.EQ.0.AND.N.LT.2.OR.M.EQ.1.AND.N.LT.3) GO TO 30
      N2 = 2*N
      RLEG(N1) = DRTS(N2+1)*DIRT(N+M)*DIRT(N-M)*(DRTS(N2-1)*COTHET*RLEG(N1-1)-DRTS(N+M-1)*DRTS(N-M-1)*DIRT(N2-3)*RLEG(N1-2))
      GO TO 30
   30 CONTINUE
   35 IF(IFLAG.EQ.1) RETURN
!      IF(SITHET.EQ.0.D0) WRITE(6,99)
!   99 FORMAT(//' *** LEGFDN  DOES NOT COMPUTE DERIVATIVES AT THE POLES
!     2 *****************'//)
      IF(SITHET.EQ.0.D0) RETURN
!
!            COMPUTE ALL THE DERIVATIVES OF THE LEGENDRE FUNCTIONS.
!
      RLNN(1) = 0.D0
      RLN = RLNN(2)
      RLNN(2) = DRTS(3)*COTHET
      DO 40      N1 = 3, M1
      N = N1-1
      N2 = 2*N
      RLN1 = RLNN(N1)
      RLNN(N1) = DRTS(N2+1)*DIRT(N2)*(SITHET*RLNN(N)+COTHET*RLN)
      RLN = RLN1
   40 CONTINUE
      DLEG(M1) = RLNN(M1)
        IF(M2.GT.NMX1) RETURN
      DO 60      N1 = M2,NMX1
      N = N1-1
      N2 = N*2
      DLEG(N1) = SITHI*( N   *RLEG(N1)*COTHET-DRTS(N-M)*DRTS(N+M)*DRTS(N2+1)*DIRT(N2-1)*RLEG(N))
   60 CONTINUE
      RETURN
      END
