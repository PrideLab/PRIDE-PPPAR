      SUBROUTINE PMSDNUT2 (RMJD, PM)
*+
*  - - - - - - - - - - -
*   P M S D N U T 2
*  - - - - - - - - - - -
*
*  This routine is part of the International Earth Rotation and
*  Reference Systems Service (IERS) Conventions software collection.
*
*  This subroutine evaluates the model of polar motion for
*  a nonrigid Earth due to tidal gravitation. This polar motion
*  is equivalent to the so-called "subdiurnal nutation." The model
*  is a sum of a first order polynomial and 25 trigonometric terms
*  (15 long periodic and 10 quasi diurnal) with coefficients given
*  in Table 5.1a of the IERS Conventions (2010).
*
*     :------------------------------------------:
*     :                                          :
*     :                 IMPORTANT                :
*     :                                          :
*     : In the present version this subroutine   :
*     : neglects the linear trend and the long   :
*     : periodic terms of the expansion, for the :
*     : reasons explained in Section 5.5.1.1 of  :
*     : the IERS Conventions (2010).             :
*     : If the full expansion is needed, set the :
*     : parameter iband to 0 instead of 1,       :
*     : that is, replace the statement           :
*     :     PARAMETER ( iband = 1 )              :
*     : to  PARAMETER ( iband = 0 )              :            
*     :                                          :
*     :__________________________________________:
*
*  In general, Class 1, 2, and 3 models represent physical effects that
*  act on geodetic parameters while canonical models provide lower-level
*  representations or basic computations that are used by Class 1, 2, or
*  3 models.
* 
*  Status:  Class 1 model
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
*     pm          d(2)      Vector of length 2 (Note 1)
*
*  Notes:
*
*  1) The polar motion coordinates (dx, dy) are expressed in
*     microarcseconds.
*
*  Called:
*     FUNDARG             Compute the angular fundamental arguments
*
*  Test case:
*     given input: rmjd = 54335D0 ( August 23, 2007 ) 
*
*     expected output: (dx) pm(1)  = 24.83144238273364834D0 microarcseconds
*                      (dy) pm(2) = -14.09240692041837661D0 microarcseconds
*
*  References:
*
*     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
*     IERS Technical Note No. 36, BKG (2010)
*
*  Revisions:
*  2005 March       A.Brzezinski   Original code
*  2008 November 26 B.E.Stetzler   Initial changes to code
*  2008 December 01 B.E.Stetzler   Provided test case
*  2009 August   18 B.E.Stetzler   Capitalized all variables for FORTRAN
*                                  77 compatibility
*  2010 May      14 B.E.Stetzler   Replaced call to PMARGS to FUNDARG
*                                  for universal fundamental argument
*                                  subroutine
*  2010 May      17 B.E.Stetzler   Validated test case using internally
*                                  computed GMST and call to FUNDARG
*                                  matched previous external call to
*                                  PMARGS for all six parameters
*  2010 June     23 B.E.Stetzler   Modified coefficients of the long
*                                  and short period terms in polar 
*                                  motion and secular polar motion
*                                  rate to coincide with Table 5.1a
*  2011 October  13 B.E. Stetzler  Modification above yields difference
*                                  in output between the original PMsdnut
*                                  and PMSDNUT2.F at the 10^-2 muas level
*                                  Test case header corrected to match
*                                  output from PMSDNUT2.F
*-----------------------------------------------------------------------

      IMPLICIT NONE
      DOUBLE PRECISION RMJD, PM(2)

*         ----------------------------
*           D E F I N I T I O N S
*         ----------------------------
*  iband  - parameter defining the range of periods for the terms which
*           are included in computations; if equal to 1 only the quasi 
*           diurnal terms are computed, otherwise the full model
*  iarg   - array defining for each of the 25 trigonometric terms a set
*           of 6 integer multipliers of the fundamental angular arguments
*  arg    - vector of the following 6 fundamental arguments used to
*           compute the angular argument of the trigonometric functions
*           arg(1:6) = [ GMST+pi, el, elp, f, d, om ]; this vector is
*           evaluated by the subroutine FUNDARG which is called as an 
*           external subroutine.  Originally evaluated by the subroutine
*           PMARGS. 
*  period - array of periods of the trigonometric terms of expansion, in
*           mean solar days; only for a check - not used in computations
*  xs, xc - sine and cosine coefficients of the x coordinate of the pole,
*           in microarcseconds
*  ys, yc - sine and cosine coefficients of the y coordinate of the pole,
*           in microarcseconds
*  angle  - angular argument of the trigonometric functions
*           angle = Sum(i=1:6) iarg(i,j)*arg(i), for j=1,25

      INTEGER IBAND, I, J, JSTART
      PARAMETER ( IBAND = 1 )
      INTEGER IARG(6,25)
      DOUBLE PRECISION T, GMST, L, LP, F, D, OM
      DOUBLE PRECISION ARG(6)
      DOUBLE PRECISION PER(25), XS(25), XC(25), YS(25), YC(25)
      DOUBLE PRECISION ANGLE, XRATE, YRATE

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

* Coefficients of the long periodic terms in polar motion
* Source: IERS Conventions (2010), Table 5.1a

      DATA 
     . (  (IARG(I,J),I=1,6),    PER(J),   XS(J),  XC(J),  YS(J),  YC(J),
     .                                                           J=1,15)
     ./ 0,  0, 0,  0,  0, -1, 6798.3837d0,  0.0d0, 0.6d0,-0.1d0,-0.1d0,
     .  0, -1, 0,  1,  0,  2, 6159.1355d0,  1.5d0, 0.0d0,-0.2d0, 0.1d0,
     .  0, -1, 0,  1,  0,  1, 3231.4956d0,-28.5d0,-0.2d0, 3.4d0,-3.9d0,
     .  0, -1, 0,  1,  0,  0, 2190.3501d0, -4.7d0,-0.1d0, 0.6d0,-0.9d0,
     .  0,  1, 1, -1,  0,  0, 438.35990d0, -0.7d0, 0.2d0,-0.2d0,-0.7d0,
     .  0,  1, 1, -1,  0, -1, 411.80661d0,  1.0d0, 0.3d0,-0.3d0, 1.0d0,
     .  0,  0, 0,  1, -1,  1, 365.24219d0,  1.2d0, 0.2d0,-0.2d0, 1.4d0,
     .  0,  1, 0,  1, -2,  1, 193.55971d0,  1.3d0, 0.4d0,-0.2d0, 2.9d0,
     .  0,  0, 0,  1,  0,  2, 27.431826d0, -0.1d0,-0.2d0, 0.0d0,-1.7d0,
     .  0,  0, 0,  1,  0,  1, 27.321582d0,  0.9d0, 4.0d0,-0.1d0,32.4d0,
     .  0,  0, 0,  1,  0,  0, 27.212221d0,  0.1d0, 0.6d0, 0.0d0, 5.1d0,
     .  0, -1, 0,  1,  2,  1, 14.698136d0,  0.0d0, 0.1d0, 0.0d0, 0.6d0,
     .  0,  1, 0,  1,  0,  1, 13.718786d0, -0.1d0, 0.3d0, 0.0d0, 2.7d0,
     .  0,  0, 0,  3,  0,  3, 9.1071941d0, -0.1d0, 0.1d0, 0.0d0, 0.9d0,
     .  0,  0, 0,  3,  0,  2, 9.0950103d0, -0.1d0, 0.1d0, 0.0d0, 0.6d0/

* Coefficients of the quasi diurnal terms in polar motion
* Source: IERS Conventions (2010), Table 5.1a

      DATA 
     .(  (IARG(I,J),i=1,6),     PER(J),   XS(J),  XC(J),  YS(J),  YC(J),
     .                                                          J=16,25)
     ./ 1, -1, 0, -2,  0, -1,1.1196992d0, -0.4d0, 0.3d0,-0.3d0, -0.4d0,
     .  1, -1, 0, -2,  0, -2,1.1195149d0, -2.3d0, 1.3d0,-1.3d0, -2.3d0,
     .  1,  1, 0, -2, -2, -2,1.1134606d0, -0.4d0, 0.3d0,-0.3d0, -0.4d0,
     .  1,  0, 0, -2,  0, -1,1.0759762d0, -2.1d0, 1.2d0,-1.2d0, -2.1d0,
     .  1,  0, 0, -2,  0, -2,1.0758059d0,-11.4d0, 6.5d0,-6.5d0,-11.4d0,
     .  1, -1, 0,  0,  0,  0,1.0347187d0,  0.8d0,-0.5d0, 0.5d0,  0.8d0,
     .  1,  0, 0, -2,  2, -2,1.0027454d0, -4.8d0, 2.7d0,-2.7d0, -4.8d0,
     .  1,  0, 0,  0,  0,  0,0.9972696d0, 14.3d0,-8.2d0, 8.2d0, 14.3d0,
     .  1,  0, 0,  0,  0, -1,0.9971233d0,  1.9d0,-1.1d0, 1.1d0,  1.9d0,
     .  1,  1, 0,  0,  0,  0,0.9624365d0,  0.8d0,-0.4d0, 0.4d0,  0.8d0/

* Rate of secular polar motion, in microarcseconds per year
* Source: IERS Conventions (2010), Table 5.1a

      DATA XRATE, YRATE / -3.8d0, -4.3d0/

* Compute the periodical part of the model
* Coordinates of the pole are set to zero first
      PM(1) = 0D0
      PM(2) = 0D0

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

      IF (IBAND.EQ.1) THEN
        JSTART = 16
      ELSE
        JSTART = 1
      ENDIF
      DO 20 J=JSTART,25

* For the j-th term of the trigonometric expansion, compute the angular
* argument angle of sine and cosine functions as a linear integer
* combination of the 6 fundamental arguments
        ANGLE = 0D0
        DO 10 I=1,6
          ANGLE = ANGLE + IARG(I,J) * ARG(I)
   10   CONTINUE
        ANGLE = DMOD( ANGLE, TWOPI )


* Compute contribution from the j-th term to the polar motion coordinates
        PM(1) = PM(1) + XS(J)*DSIN(ANGLE) + XC(J)*DCOS(ANGLE)
        PM(2) = PM(2) + YS(J)*DSIN(ANGLE) + YC(J)*DCOS(ANGLE)
   20 CONTINUE
      IF (IBAND.EQ.1) RETURN

* Add the secular term of the model
      PM(1) = PM(1) + XRATE * (RMJD-RMJD0) / 365.25D0
      PM(2) = PM(2) + YRATE * (RMJD-RMJD0) / 365.25D0

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
