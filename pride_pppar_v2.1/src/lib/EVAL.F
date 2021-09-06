      REAL FUNCTION EVAL (Y,NN,X,U,S)
*+
*  - - - - - - - - - - -
*   E V A L 
*  - - - - - - - - - - -
*
*  This routine is part of the International Earth Rotation and
*  Reference Systems Service (IERS) Conventions software collection.
*
*  This function performs cubic spline interpolation of a given function
*  sampled at unequally spaced intervals.  The subroutine SPLINE needs
*  to be called beforehand to set up the array s.
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
*     y            d     the coordinate at which a function value is
*                        desired (Note 1)
*     nn           i     number of samples of the original function
*     x            d     array containing sample coordinates x(1),x(2),...
*                        x(nn) (Note 2)
*     s            d     array containing the 2nd derivatives at the sample
*                        points (Note 3)
*
*  Returned:
*     u            d     array containing samples of a function at the
*                        coordinates x(1),x(2),...x(nn)
*
*  Notes:
*
*  1) If y falls outside the range (x(1),x(nn)), the value at the nearest
*     endpoint of the series is used.
*
*  2) The sequence x(1),x(2),...x(nn) must be strictly increasing.
*
*  3) This array is found by the subroutine SPLINE, which must be called
*     once before beginning this interpolation.
*
*  Called:
*     None
*
*  Test case:
*     
*  Not provided for this function.  This is a support routine of the main
*  program HARDISP.F.
*
*  References:
*
*     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
*     IERS Technical Note No. 36, BKG (2010)
*
*  Revisions:
*  2009 June   03 B.E.Stetzler Initial standardization of function
*  2009 August 19 B.E.Stetzler Capitalized all variables for FORTRAN 77
*                               compatibility
*-----------------------------------------------------------------------
      
      IMPLICIT NONE
      INTEGER K,K1,K2,NN
      REAL DELI,DK,DY,DY1,F1,F2,F3,FF1,FF2,S,U,X,Y
      DIMENSION X(*),U(*),S(*)

      NN = IABS(NN)

*     If y is out of range, substitute endpoint values

      IF (Y.LE.X(1)) THEN
         EVAL=U(1)
         RETURN
      ENDIF

      IF (Y.GE.X(NN)) THEN
         EVAL=U(NN)
         RETURN
      ENDIF

*    Locate interval (x(k1),x(k2)) which contains y
      DO 100 K=2,NN
         IF(X(K-1).LT.Y.AND.X(K).GE.Y) THEN
           K1=K-1
           K2=K
         ENDIF
100   CONTINUE

*    Evaluate and then interpolate
      DY=X(K2)-Y
      DY1=Y-X(K1)
      DK=X(K2)-X(K1)
      DELI=1.0D0/(6.0D0*DK)
      FF1=S(K1)*DY*DY*DY
      FF2=S(K2)*DY1*DY1*DY1
      F1=(FF1+FF2)*DELI
      F2=DY1*((U(K2)/DK)-(S(K2)*DK)/6.0D0)
      F3= DY*((U(K1)/DK)-(S(K1)*DK)/6.0D0)
      EVAL=F1+F2+F3
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
