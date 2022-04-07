      SUBROUTINE UTLIBR (RMJD, DUT1, DLOD)
*+
*  - - - - - - - - 
*   U T L I B R
*  - - - - - - - -
*
*  This routine is part of the International Earth Rotation and
*  Reference Systems Service (IERS) Conventions software collection.
*
*  This subroutine evaluates the model of subdiurnal libration
*  in the axial component of rotation, expressed by UT1 and LOD.
*  This effect is due to the influence of tidal gravitation on the
*  departures of the Earth's mass distribution from the rotational
*  symmetry, expressed by the non-zonal components of geopotential.
*  The amplitudes have been computed for an elastic Earth with liquid
*  core. The adopted truncation level is 0.033 microseconds in UT1
*  corresponding to the angular displacement of 0.5 microarcseconds
*  or to 0.015 mm at the planet surface. With this truncation level
*  the model contains 11 semidiurnal terms. The coefficients of
*  the model are given in Table 5.1b of the IERS Conventions (2010).
*
*  In general, Class 1, 2, and 3 models represent physical effects that
*  act on geodetic parameters while canonical models provide lower-level
*  representations or basic computations that are used by Class 1, 2, or
*  3 models.
* 
*  Status:  Class 3 model
*
*     Class 1 models are those recommended to be used a priori in the
*     reduction of raw space geodetic data in order to determine
*     geodetic parameter estimates.
*     Class 2 models are those that eliminate an observational
*     singularity and are purely conventional in nature.
*     Class 3 models are those that are not required as either Class
*     1 or 2.
*     Canonical models are accepted as is and cannot be classified as
*     a Class 1, 2, or 3 model.
*
*  Given:
*     rmjd        d      Time expressed as modified Julian date
*
*  Returned:
*     dUT1        d      Incremental UT1 in microseconds
*     dLOD        d      Incremental LOD in microseconds per day
*
*  Notes:
*  1) The procedure FUNDARG.F is the same as used by the program PMSDNUT2.F
*     which implements the corresponding model of the lunisolar libration in
*     polar motion.
*
*  Called:
*     FUNDARG             Compute the angular fundamental arguments
*
*  Test cases:
*     given input:  rmjd_a = 44239.1 ( January 1, 1980 2:24.00 )
*                   rmjd_b = 55227.4 ( January 31, 2010 9:35.59 )
*
*     expected output: dUT1_a =   2.441143834386761746D0 mus;
*                      dLOD_a = -14.78971247349449492D0 mus / day
*                      dUT1_b = - 2.655705844335680244D0 mus;
*                      dLOD_b =  27.39445826599846967D0 mus / day
*
*  References:
*
*     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
*     IERS Technical Note No. 36, BKG (2010)
*
*  Revisions:
*  2010 May       A.Brzezinski   Original code
*  2010 June  1   B.E.Stetzler   Initial changes to code
*  2010 June  2   B.E.Stetzler   Provided test case
*  2010 June  2   B.E.Stetzler   Capitalized all variables for FORTRAN
*                                77 compatibility
*  2010 June  2   B.E.Stetzler   Replaced call to PMARGS to FUNDARG
*                                for universal fundamental argument
*                                subroutine
*  2010 June  2   B.E.Stetzler   Validated test case using internally
*                                computed GMST and call to FUNDARG
*                                matched previous external call to
*                                PMARGS for all six parameters
*  2010 June  23  B.E.Stetzler   Modified coefficients of semi-diurnal
*                                variations in UT1 and LOD due to
*                                libration for a non-rigid Earth to
*                                coincide with Table 5.1b
*-----------------------------------------------------------------------

      IMPLICIT NONE
      DOUBLE PRECISION RMJD, DUT1, DLOD

*         ----------------------------
*           D E F I N I T I O N S
*         ----------------------------
*  iarg   - array defining for each of the 11 trigonometric terms a set
*           of 6 integer multipliers of the fundamental angular arguments
*  arg    - vector of the following 6 fundamental arguments used to
*           compute the angular argument of the trigonometric functions
*           arg(1:6) = [ GMST+pi, el, elp, f, d, om ]; this vector is
*           evaluated by the subroutine FUNDARG which is called as an 
*           external subroutine.  Originally evaluated by the subroutine
*           PMARGS. 
*  period - array of periods of the trigonometric terms of expansion, in
*           mean solar days; only for a check - not used in computations
*  dUT1s, dUT1c - sine and cosine coefficients of dUT1, in microseconds
*  dLODs, dLODc - sine and cosine coefficients of dLOD, in microseconds
*                 per day
*  angle  - angular argument of the trigonometric functions
*           angle = Sum(i=1:6) iarg(i,j)*arg(i), for j=1,11

      INTEGER I, J
      INTEGER IARG(6,11)
      DOUBLE PRECISION T, GMST, L, LP, F, D, OM
      DOUBLE PRECISION ARG(6)
      DOUBLE PRECISION PER(11), DUT1S(11), DUT1C(11), DLODS(11),
     .                 DLODC(11) 
      DOUBLE PRECISION ANGLE

* Set constants

*  Arcseconds to radians
      DOUBLE PRECISION DAS2R
      PARAMETER ( DAS2R = 4.848136811095359935899141D-6 )

*  Arcseconds in a full circle
      DOUBLE PRECISION TURNAS
      PARAMETER ( TURNAS = 1296000D0 )

*  rmjd0   - modified Julian date of J2000
*  twopi   - 2*pi

      DOUBLE PRECISION RMJD0, PI, TWOPI
      PARAMETER ( RMJD0   = 51544.5D0                )
      PARAMETER ( PI      = 3.141592653589793238462643D0 )
      PARAMETER ( TWOPI   = 6.283185307179586476925287D0 )

*  Radians to seconds
      DOUBLE PRECISION RAD2SEC
      PARAMETER ( RAD2SEC = 86400D0/TWOPI            )

* Coefficients of the quasi semidiurnal terms in dUT1, dLOD 
* Source: IERS Conventions (2010), Table 5.1b

      DATA 
     .((IARG(I,J),I=1,6), PER(J), DUT1S(J),DUT1C(J), DLODS(J), DLODC(J),
     .                                                           J=1,11)
     ./2,-2, 0,-2, 0,-2,0.5377239d0, 0.05d0,-0.03d0, -0.3d0, -0.6d0,
     . 2, 0, 0,-2,-2,-2,0.5363232d0, 0.06d0,-0.03d0, -0.4d0, -0.7d0,
     . 2,-1, 0,-2, 0,-2,0.5274312d0, 0.35d0,-0.20d0, -2.4d0, -4.1d0,
     . 2, 1, 0,-2,-2,-2,0.5260835d0, 0.07d0,-0.04d0, -0.5d0, -0.8d0,
     . 2, 0, 0,-2, 0,-1,0.5175645d0,-0.07d0, 0.04d0,  0.5d0,  0.8d0,
     . 2, 0, 0,-2, 0,-2,0.5175251d0, 1.75d0,-1.01d0,-12.2d0,-21.3d0,
     . 2, 1, 0,-2, 0,-2,0.5079842d0,-0.05d0, 0.03d0,  0.3d0,  0.6d0,
     . 2, 0,-1,-2, 2,-2,0.5006854d0, 0.04d0,-0.03d0, -0.3d0, -0.6d0,
     . 2, 0, 0,-2, 2,-2,0.5000000d0, 0.76d0,-0.44d0, -5.5d0, -9.6d0,
     . 2, 0, 0, 0, 0, 0,0.4986348d0, 0.21d0,-0.12d0, -1.5d0, -2.6d0,
     . 2, 0, 0, 0, 0,-1,0.4985982d0, 0.06d0,-0.04d0, -0.4d0, -0.8d0/

* Compute the harmonic model of dUT1 and dLOD 
* dUT1 and dLOD are set to zero first 
      DUT1 = 0D0
      DLOD = 0D0

* Evaluate the vector of the fundamental arguments
* arg(1:6) = [ GMST+pi, el, elp, f, d, om ] at t = rmjd

*  Convert the input epoch to Julian centuries of TDB since J2000
      T = (RMJD-RMJD0)/36525D0

*  Compute GMST + pi
      GMST = MOD (   67310.54841D0 +
     .               T*( (8640184.812866D0 + 3155760000D0) +
     .               T*( 0.093104D0 +
     .               T*( -0.0000062d0 ))), 86400D0 )

      CALL FUNDARG ( T, L, LP, F, D, OM )

      ARG(1) = GMST / RAD2SEC + PI
      ARG(1) = DMOD( ARG(1), TWOPI )
      ARG(2) = L
      ARG(3) = LP
      ARG(4) = F
      ARG(5) = D
      ARG(6) = OM 

      DO 20 J=1,11

* For the j-th term of the trigonometric expansion, compute the angular
* argument angle of sine and cosine functions as a linear integer
* combination of the 6 fundamental arguments
        ANGLE = 0D0
        DO 10 I=1,6
          ANGLE = ANGLE + IARG(I,J) * ARG(I)
   10   CONTINUE
        ANGLE = DMOD( ANGLE, TWOPI )

* Compute contribution from the j-th term of expansion to dUT1 and dLOD 
        DUT1 = DUT1 + DUT1S(J)*DSIN(ANGLE) + DUT1C(J)*DCOS(ANGLE)
        DLOD = DLOD + DLODS(J)*DSIN(ANGLE) + DLODC(J)*DCOS(ANGLE)
   20 CONTINUE
      RETURN

*  Finished.

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
*     to others, provided that the derived work complies with the
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
