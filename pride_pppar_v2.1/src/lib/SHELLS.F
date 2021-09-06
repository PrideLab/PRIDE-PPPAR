      SUBROUTINE SHELLS (X,K,N)
*+
*  - - - - - - - - -
*   S H E L L S
*  - - - - - - - - -
*
*  This routine is part of the International Earth Rotation and
*  Reference Systems Service (IERS) Conventions software collection.
*
*  The subroutine sorts an array x, of length n, sorting upward,
*  and returns an array k which may be used to key another array
*  to the sorted pattern (i.e., if we had an array f to which x
*  corresponded before sorting, then after calling SHELLS,
*  f(k(1)) will be the element of f corresponding to the
*  smallest x, f(k(2)) the next smallest, and so on).
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
*     x              d      array to be sorted (Note 1)
*     n              i      length of the input array x
*
*  Returned:
*     k              i      sorted array that may be used to key another 
*                           array
*  Notes:
*
*  1) See the subroutine ADMINT.F header comments for detailed information.
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
*  1982 December 29              Revised so that array k is sorted in turn
*                                
*  2009 June   05 B.E. Stetzler    Added header and copyright
*  2009 August 19 B.E. Stetzler    Capitalized all variables for FORTRAN
*                                  77 compatibility
*-----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER I,IGAP,IEX,IK,IMAX,IPL,J,K,L,N
      REAL SV,X
      DIMENSION X(N),K(N)

      IGAP = N

      DO 1 I = 1,N
 1    K(I) = I
 5    IF(IGAP.LE.1) GO TO 25

      IGAP = IGAP/2
      IMAX = N - IGAP
 10   IEX = 0
      DO 20 I = 1,IMAX
      IPL = I + IGAP
      IF(X(I).LE.X(IPL)) GO TO 20
      SV = X(I)
      IK = K(I)
      X(I) = X(IPL)
      K(I) = K(IPL)
      X(IPL) = SV
      K(IPL) = IK
      IEX = IEX + 1
 20   CONTINUE

      IF(IEX.GT.0) GO TO 10
      GO TO 5

*  Now sort k's (for identical values of x, if any)

 25   J = 1
 30   IF(J.GE.N) RETURN
      IF(X(J).EQ.X(J+1)) GO TO 33
      J = J + 1
      GO TO 30
*  Have at least two x's with the same value. See how long this is true
 33   L = J
 35   IF(X(L).NE.X(L+1)) GO TO 38
      L = L + 1
      IF(L.LT.N) GO TO 35
*  j and l are the indices within which x(i) does not change - sort k
 38   IGAP = L - J + 1
 40   IF(IGAP.LE.1) J = L + 1
      IF(IGAP.LE.1) GO TO 30

      IGAP = IGAP/2
      IMAX = L-J+1 - IGAP
 45   IEX = 0

      DO 50 I=1,IMAX
      IPL = I + IGAP + J - 1
      IF(K(I+J-1).LE.K(IPL)) GO TO 50
      IK = K(I+J-1)
      K(I+J-1) = K(IPL)
      K(IPL) = IK
      IEX = IEX + 1
 50   CONTINUE
      IF(IEX.GT.0) GO TO 45
      GO TO 40

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
