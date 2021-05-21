      SUBROUTINE SPLINE (NN,X,U,S,A)
*+
*  - - - - - - - - -
*   S P L I N E
*  - - - - - - - - -
*
*  This routine is part of the International Earth Rotation and
*  Reference Systems Service (IERS) Conventions software collection.
*
*  The purpose of the subroutine is to find an array s for the spline
*  interpolator function EVAL.
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
*  Given: This is a support routine of the main program HARDISP.F.
*     nn             i      number of data points supplied, which may be
*                           negative (Note 1)
*     x              d      array containing x-coordinates where function
*                           is sampled (Note 2)
*     u              d      array containing sample values that are to be
*                           interpolated 
*     a              d      working space array of dimension at least nn
*
*  Returned:
*     s              d      output array of 2nd derivative at sample points 
*
*  Notes:
*
*  1) If the user wishes to force the derivatives at the ends of the series
*     to assume specified values, he or she should put du(1)/dx and du(n)/dx
*     in the variables s1 and s2 and call the subroutine with nn = -(number
*     of terms in the series).  Normally a parabola is fitted through the 
*     1st and last 3 points to find the slopes.  If less than 4 points are
*     supplied, straight lines are fitted.
* 
*  2) The sequence xx(1), xx(2), ... xx(nn) must be strictly increasing.
*
*  Called:
*     None
*
*  Test case:
*     Not provided for this subroutine.
*
*  References:
*
*     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
*     IERS Technical Note No. 36, BKG (2010)
*
*  Revisions:
*  2009 June   08 B.E.Stetzler    Added header and copyright
*  2009 August 19 B.E.Stetzler    Capitalized all variables for FORTRAN
*                                 77 compatibility
*  2009 August 26 B.E.Stetzler    Used IMPLICIT NONE and defined all variables
*-----------------------------------------------------------------------

      IMPLICIT NONE
      INTEGER I,J,N,N1,NN,NMAX
      PARAMETER (NMAX = 20)
      REAL A,C,Q,Q1,QN,S,U,X,U1,U2,X1,X2
      DIMENSION X(NMAX),U(NMAX),S(NMAX),A(NMAX)

      Q(U1,X1,U2,X2)=(U1/X1**2-U2/X2**2)/(1.0/X1-1.0/X2)
      

      N = IABS(NN)

      IF (N.LE.3) THEN

*  series too short for cubic spline - use straight lines.
         DO I=1,N
            S(I)=0.0
         ENDDO
         RETURN
      ENDIF

      Q1=Q(U(2)-U(1),X(2)-X(1),U(3)-U(1),X(3)-X(1))
      QN=Q(U(N-1)-U(N),X(N-1)-X(N),U(N-2)-U(N),X(N-2)-X(N))

      IF (NN.LE.0) THEN
         Q1=S(1)
         QN=S(2)
      ENDIF

      S(1)=6.0*((U(2)-U(1))/(X(2)-X(1)) - Q1)
      N1= N - 1

      DO I=2,N1
         S(I)= (U(I-1)/(X(I)-X(I-1)) - U(I)*(1.0/(X(I)-X(I-1))
     .   + 1.0/(X(I+1)-X(I))) + U(I+1)/(X(I+1)-X(I)))*6.0
      ENDDO

      S(N)=6.0*(QN + (U(N1)-U(N))/(X(N)-X(N1)))
      A(1)=2.0*(X(2)-X(1))
      A(2)=1.5*(X(2)-X(1)) + 2.0*(X(3)-X(2))
      S(2)=S(2) - 0.5*S(1)

      DO I=3,N1
         C=(X(I)-X(I-1))/A(I-1)
         A(I)=2.0*(X(I+1)-X(I-1)) - C*(X(I)-X(I-1))
         S(I)=S(I) - C*S(I-1)
      ENDDO

      C=(X(N)-X(N1))/A(N1)
      A(N)=(2.0-C)*(X(N)-X(N1))
      S(N)=S(N) - C*S(N1)

*  Back substitute
      S(N)= S(N)/A(N)

      DO J=1,N1
         I=N-J
         S(I) =(S(I) - (X(I+1)-X(I))*S(I+1))/A(I)
      ENDDO
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
