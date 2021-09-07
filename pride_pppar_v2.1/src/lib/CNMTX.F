      SUBROUTINE CNMTX ( DMJD,H )
*+
*  - - - - - - - - - -
*   C N M T X 
*  - - - - - - - - - -
*
*  This routine is part of the International Earth Rotation and
*  Reference Systems Service (IERS) Conventions software collection.
*
*  The purpose of the subroutine is to compute the time dependent part
*  of second degree diurnal and semidiurnal tidal potential from
*  the dominant spectral lines in the Cartwright-Tayler-Edden harmonic
*  decomposition.
*   
*  In general, Class 1, 2, and 3 models represent physical effects that
*  act on geodetic parameters while canonical models provide lower-level
*  representations or basic computations that are used by Class 1, 2, or
*  3 models.
*
*  Status: Canonical model
* 
*     Class 1 models are those recommended to be used a priori in the
*     reduction of raw space geodetic data in order to determine
*     geodetic parameter estimates.
*     Class 2 models are those that eliminate an observational
*     singularity and are purely conventional in nature.
*     Class 3 models are those that are not required as either Class
*     1 or 2.
*     Canonical models are accepted as is and cannot be classified as a
*     Class 1, 2, or 3 model.
*
*  Given:
*     dmjd           d     Modified Julian Date 
*
*  Returned:
*     h              d     vector of length 12 with partials of the 
*                          tidal variation with respect to the
*                          orthoweights (Note 1) 
*
*  Notes:
*
*  1) The diurnal and semidiurnal orthoweights fit to the 8 constituents
*     are listed in Reference Ray et al.          
*
*  Test case:
*     given input: dmjd = 54964.0D0
*
*     expected output: h(1) = 15.35873641938967360D0
*                      h(2) = 9.784941251812741214D0
*                      h(3) = -5.520740128266865554D0
*                      h(4) = 3.575314211234633888D0
*                      h(5) = -13.93717453496387648D0
*                      h(6) = -9.167400321705855504D0
*                      h(7) = 5.532815475865292321D0
*                      h(8) = 9.558741883500834646D0
*                      h(9) = -10.22541212627272600D0
*                      h(10)= 0.8367570529461261231D0
*                      h(11)= 1.946355176475630611D0
*                      h(12)= -13.55702062247304696D0
*
*  References:
*
*     Ray,R. D., Steinberg, D. J., Chao, B. F., and Cartwright, D. E.,
*     "Diurnal and Semidiurnal Variations in the Earth's Rotation
*     Rate Induced by Ocean Tides", 1994, Science, 264, pp. 830-832
*
*     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
*     IERS Technical Note No. 36, BKG (2010)
*
*  Revisions:
*  2008 November 07 B.E. Stetzler    Added header and copyright
*  2008 November 21 B.E. Stetzler    Redefined variables as double
*                                    precision and changed twopi
*  2009 May      12 B.E. Stetzler    Added D0 to DATA structure values
*  2009 May      19 B.E. Stetzler    Provided validated test case
*  2009 June     08 B.E. Stetzler    Redefined pinm and alpha, used
*                                    double precision exclusively, and
*                                    updated validated test case values
*  2009 September 03 B.E.Stetzler    Capitalized all variables for Fortran
*                                    77 compatibility
*  2010 March     17 B.E.Stetzler    Aligned table array values in DATA
*                                    structures
*-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER I,J,K,M,N,NLINES,NMAX
      PARAMETER (NLINES=71)
      CHARACTER*7 NUMARG(NLINES)
      DOUBLE PRECISION H(12),DMJD, DT60, D1960, DT, TWOPI
      DOUBLE PRECISION NJ(NLINES),MJ(NLINES)
      DOUBLE PRECISION HS(NLINES),PHASE(NLINES),FREQ(NLINES)
      DOUBLE PRECISION ANM(2:3,0:3,-1:1),BNM(2:3,0:3,-1:1)
      DOUBLE PRECISION AP, AM, BP, BM, PINM, ALPHA
      DOUBLE PRECISION P(0:2,2),Q(0:2,2),SP(6,2)

* Define the orthotide weight factors
      DATA ((SP(I,M),I=1,6),M=1,2) /
     . 0.0298D0,0.1408D0,+0.0805D0, 0.6002D0,+0.3025D0, 0.1517D0,
     . 0.0200D0,0.0905D0,+0.0638D0, 0.3476D0,+0.1645D0, 0.0923D0/

      PARAMETER (TWOPI = 6.283185307179586476925287D0)
      DATA DT /2D0/
      DATA NMAX /2/

* tidal potential model for 71 diurnal and semidiurnal lines

      DATA D1960/37076.5D0/
      DATA (NJ(J),MJ(J),HS(J),PHASE(J),FREQ(J),NUMARG(J),J=1,15)
     ./2D0,1D0,  -1.94D0,  9.0899831D0,  5.18688050D0,'117.655',
     . 2D0,1D0,  -1.25D0,  8.8234208D0,  5.38346657D0,'125.745',
     . 2D0,1D0,  -6.64D0, 12.1189598D0,  5.38439079D0,'125.755',
     . 2D0,1D0,  -1.51D0,  1.4425700D0,  5.41398343D0,'127.545',
     . 2D0,1D0,  -8.02D0,  4.7381090D0,  5.41490765D0,'127.555',
     . 2D0,1D0,  -9.47D0,  4.4715466D0,  5.61149372D0,'135.645',
     . 2D0,1D0, -50.20D0,  7.7670857D0,  5.61241794D0,'135.655',
     . 2D0,1D0,  -1.80D0, -2.9093042D0,  5.64201057D0,'137.445',
     . 2D0,1D0,  -9.54D0,  0.3862349D0,  5.64293479D0,'137.455',
     . 2D0,1D0,   1.52D0, -3.1758666D0,  5.83859664D0,'145.535',
     . 2D0,1D0, -49.45D0,  0.1196725D0,  5.83952086D0,'145.545',
     . 2D0,1D0,-262.21D0,  3.4152116D0,  5.84044508D0,'145.555',
     . 2D0,1D0,   1.70D0, 12.8946194D0,  5.84433381D0,'145.755',
     . 2D0,1D0,   3.43D0,  5.5137686D0,  5.87485066D0,'147.555',
     . 2D0,1D0,   1.94D0,  6.4441883D0,  6.03795537D0,'153.655'/
      DATA (NJ(J),MJ(J),HS(J),PHASE(J),FREQ(J),NUMARG(J),J=16,30)
     ./2D0,1D0,   1.37D0, -4.2322016D0,  6.06754801D0,'155.445',
     . 2D0,1D0,   7.41D0, -0.9366625D0,  6.06847223D0,'155.455',
     . 2D0,1D0,  20.62D0,  8.5427453D0,  6.07236095D0,'155.655',
     . 2D0,1D0,   4.14D0, 11.8382843D0,  6.07328517D0,'155.665',
     . 2D0,1D0,   3.94D0,  1.1618945D0,  6.10287781D0,'157.455',
     . 2D0,1D0,  -7.14D0,  5.9693878D0,  6.24878055D0,'162.556',
     . 2D0,1D0,   1.37D0, -1.2032249D0,  6.26505830D0,'163.545',
     . 2D0,1D0,-122.03D0,  2.0923141D0,  6.26598252D0,'163.555',
     . 2D0,1D0,   1.02D0, -1.7847596D0,  6.28318449D0,'164.554',
     . 2D0,1D0,   2.89D0,  8.0679449D0,  6.28318613D0,'164.556',
     . 2D0,1D0,  -7.30D0,  0.8953321D0,  6.29946388D0,'165.545',
     . 2D0,1D0, 368.78D0,  4.1908712D0,  6.30038810D0,'165.555',
     . 2D0,1D0,  50.01D0,  7.4864102D0,  6.30131232D0,'165.565',
     . 2D0,1D0,  -1.08D0, 10.7819493D0,  6.30223654D0,'165.575',
     . 2D0,1D0,   2.93D0,  0.3137975D0,  6.31759007D0,'166.554'/
      DATA (NJ(J),MJ(J),HS(J),PHASE(J),FREQ(J),NUMARG(J),J=31,45)
     ./2D0,1D0,   5.25D0,  6.2894282D0,  6.33479368D0,'167.555',
     . 2D0,1D0,   3.95D0,  7.2198478D0,  6.49789839D0,'173.655',
     . 2D0,1D0,  20.62D0, -0.1610030D0,  6.52841524D0,'175.455',
     . 2D0,1D0,   4.09D0,  3.1345361D0,  6.52933946D0,'175.465',
     . 2D0,1D0,   3.42D0,  2.8679737D0,  6.72592553D0,'183.555',
     . 2D0,1D0,   1.69D0, -4.5128771D0,  6.75644239D0,'185.355',
     . 2D0,1D0,  11.29D0,  4.9665307D0,  6.76033111D0,'185.555',
     . 2D0,1D0,   7.23D0,  8.2620698D0,  6.76125533D0,'185.565',
     . 2D0,1D0,   1.51D0, 11.5576089D0,  6.76217955D0,'185.575',
     . 2D0,1D0,   2.16D0,  0.6146566D0,  6.98835826D0,'195.455',
     . 2D0,1D0,   1.38D0,  3.9101957D0,  6.98928248D0,'195.465',
     . 2D0,2D0,   1.80D0, 20.6617051D0, 11.45675174D0,'225.855',
     . 2D0,2D0,   4.67D0, 13.2808543D0, 11.48726860D0,'227.655',
     . 2D0,2D0,  16.01D0, 16.3098310D0, 11.68477889D0,'235.755',
     . 2D0,2D0,  19.32D0,  8.9289802D0, 11.71529575D0,'237.555'/
      DATA (NJ(J),MJ(J),HS(J),PHASE(J),FREQ(J),NUMARG(J),J=46,60)
     ./2D0,2D0,   1.30D0,  5.0519065D0, 11.73249771D0,'238.554',
     . 2D0,2D0,  -1.02D0, 15.8350306D0, 11.89560406D0,'244.656',
     . 2D0,2D0,  -4.51D0,  8.6624178D0, 11.91188181D0,'245.645',
     . 2D0,2D0, 120.99D0, 11.9579569D0, 11.91280603D0,'245.655',
     . 2D0,2D0,   1.13D0,  8.0808832D0, 11.93000800D0,'246.654',
     . 2D0,2D0,  22.98D0,  4.5771061D0, 11.94332289D0,'247.455',
     . 2D0,2D0,   1.06D0,  0.7000324D0, 11.96052486D0,'248.454',
     . 2D0,2D0,  -1.90D0, 14.9869335D0, 12.11031632D0,'253.755',
     . 2D0,2D0,  -2.18D0, 11.4831564D0, 12.12363121D0,'254.556',
     . 2D0,2D0, -23.58D0,  4.3105437D0, 12.13990896D0,'255.545',
     . 2D0,2D0, 631.92D0,  7.6060827D0, 12.14083318D0,'255.555',
     . 2D0,2D0,   1.92D0,  3.7290090D0, 12.15803515D0,'256.554',
     . 2D0,2D0,  -4.66D0, 10.6350594D0, 12.33834347D0,'263.655',
     . 2D0,2D0, -17.86D0,  3.2542086D0, 12.36886033D0,'265.455',
     . 2D0,2D0,   4.47D0, 12.7336164D0, 12.37274905D0,'265.655'/
      DATA (NJ(J),MJ(J),HS(J),PHASE(J),FREQ(J),NUMARG(J),J=61,71)
     ./2D0,2D0,   1.97D0, 16.0291555D0, 12.37367327D0,'265.665',
     . 2D0,2D0,  17.20D0, 10.1602590D0, 12.54916865D0,'272.556',
     . 2D0,2D0, 294.00D0,  6.2831853D0, 12.56637061D0,'273.555',
     . 2D0,2D0,  -2.46D0,  2.4061116D0, 12.58357258D0,'274.554',
     . 2D0,2D0,  -1.02D0,  5.0862033D0, 12.59985198D0,'275.545',
     . 2D0,2D0,  79.96D0,  8.3817423D0, 12.60077620D0,'275.555',
     . 2D0,2D0,  23.83D0, 11.6772814D0, 12.60170041D0,'275.565',
     . 2D0,2D0,   2.59D0, 14.9728205D0, 12.60262463D0,'275.575',
     . 2D0,2D0,   4.47D0,  4.0298682D0, 12.82880334D0,'285.455',
     . 2D0,2D0,   1.95D0,  7.3254073D0, 12.82972756D0,'285.465',
     . 2D0,2D0,   1.17D0,  9.1574019D0, 13.06071921D0,'295.555'/

* Compute the time dependent potential matrix

      DO 50 K=-1,1
         DT60 = (DMJD - K*DT) - D1960
         ANM(2,1:2,K) = 0D0
         BNM(2,1:2,K) = 0D0
         DO 100 J=1,NLINES
            N = NJ(J)
            M = MJ(J)
            PINM = DFLOAT(MOD(N+M,2))*TWOPI/4.0D0
            !alpha = phase(j) + freq(j)*dt60 - pinm
            ALPHA = DMOD(PHASE(J) - PINM,TWOPI) 
     .              + DMOD(FREQ(J)*DT60,TWOPI)
            ANM(N,M,K) = ANM(N,M,K) + HS(J)*DCOS(ALPHA)
            BNM(N,M,K) = BNM(N,M,K) - HS(J)*DSIN(ALPHA)
100      CONTINUE
50    CONTINUE

* orthogonalize the response terms

      DO 150 M = 1,2
        AP = ANM(2,M,1) + ANM(2,M,-1)
        AM = ANM(2,M,1) - ANM(2,M,-1)
        BP = BNM(2,M,1) + BNM(2,M,-1)
        BM = BNM(2,M,1) - BNM(2,M,-1)
        P(0,M) = SP(1,M)*ANM(2,M,0)
        P(1,M) = SP(2,M)*ANM(2,M,0) - SP(3,M)*AP
        P(2,M) = SP(4,M)*ANM(2,M,0) - SP(5,M)*AP + SP(6,M)*BM
        Q(0,M) = SP(1,M)*BNM(2,M,0)
        Q(1,M) = SP(2,M)*BNM(2,M,0) - SP(3,M)*BP
        Q(2,M) = SP(4,M)*BNM(2,M,0) - SP(5,M)*BP - SP(6,M)*AM
        ANM(2,M,-1) = P(0,M)
        ANM(2,M,0) = P(1,M)
        ANM(2,M,1) = P(2,M)
        BNM(2,M,-1) = Q(0,M)
        BNM(2,M,0) = Q(1,M)
        BNM(2,M,1) = Q(2,M)
150   CONTINUE

* fill partials vector
      J = 1
      DO 200 N=2,NMAX
         DO 250 M = 1,N
            DO 300 K = -1,1
               H(J)  = ANM(N,M,K)
               H(J+1)= BNM(N,M,K)
               J = J + 2
300         CONTINUE
250      CONTINUE 
200   CONTINUE

      RETURN

* Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2008
*  IERS Conventions Center
*
*  ==================================
*  IERS Conventions Software License
*  ==================================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is provided by the IERS Conventions Center ("the
*     Center").
*
*  2. Permission is granted to anyone to use the Software for any
*     purpose, including commercial applications, free of charge,
*     subject to the conditions and restrictions listed below.
*
*  3. You (the user) may adapt the Software and its algorithms for your
*     own purposes and you may distribute the resulting "derived work"
*     to others,provided that the derived work complies with the
*     following requirements:
*
*     a) Your work shall be clearly identified so that it cannot be
*        mistaken for IERS Conventions software and that it has been
*        neither distributed by nor endorsed by the Center.
*
*     b) Your work (including source code) must contain descriptions of
*        how the derived work is based upon and/or differs from the
*        original Software.
*
*     c) The name(s) of all modified routine(s) that you distribute
*        shall be changed.
* 
*     d) The origin of the IERS Conventions components of your derived
*        work must not be misrepresented; you must not claim that you
*        wrote the original Software.
*
*     e) The source code must be included for all routine(s) that you
*        distribute.  This notice must be reproduced intact in any
*        source distribution. 
*
*  4. In any published work produced by the user and which includes
*     results achieved by using the Software, you shall acknowledge
*     that the Software was used in obtaining those results.
*
*  5. The Software is provided to the user "as is" and the Center makes
*     no warranty as to its use or performance.   The Center does not
*     and cannot warrant the performance or results which the user may
*     obtain by using the Software.  The Center makes no warranties,
*     express or implied, as to non-infringement of third party rights,
*     merchantability, or fitness for any particular purpose.  In no
*     event will the Center be liable to the user for any consequential,
*     incidental, or special damages, including any lost profits or lost
*     savings, even if a Center representative has been advised of such
*     damages, or for any claim by any third party.
*
*  Correspondence concerning IERS Conventions software should be
*  addressed as follows:
*
*                     Gerard Petit
*     Internet email: gpetit[at]bipm.org
*     Postal address: IERS Conventions Center
*                     Time, frequency and gravimetry section, BIPM
*                     Pavillon de Breteuil
*                     92312 Sevres  FRANCE
*
*     or
*
*                     Brian Luzum
*     Internet email: brian.luzum[at]usno.navy.mil
*     Postal address: IERS Conventions Center
*                     Earth Orientation Department
*                     3450 Massachusetts Ave, NW
*                     Washington, DC 20392
*
*
*-----------------------------------------------------------------------
      END
