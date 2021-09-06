      SUBROUTINE RECURS(X,N,HC,NF,OM,SCR)
*+
*  - - - - - - - - -
*   R E C U R S
*  - - - - - - - - -
*
*  This routine is part of the International Earth Rotation and
*  Reference Systems Service (IERS) Conventions software collection.
*
*  The purpose of the subroutine is to perform sine and cosine recursion
*  to fill in data x, of length n, for nf sines and cosines with frequencies
*  om. 
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
*     x              d      data provided from a file given as standard
*                           input from the MAIN program HARDISP.F (Note 1)
*     n              i      length of the data file x
*     hc             d      array containing alternating cosine and sine
*                           coefficients
*     nf             i      number of sine and cosine terms
*     om             d      sine and cosine frequencies (Note 2)  
*
*  Returned:
*     scr            d      scratch array of length 3 times nf which is
*                           returned as the recursion cr
*  Notes:
*
*  1) See the MAIN program HARDISP.F header comments for detailed information.
* 
*  2) The frequencies are normalized so that the Nyquist frequency is pi.
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
*  2009 June 05 B.E. Stetzler    Added header and copyright, used DCOS
*                                and DSIN exclusively, and replaced END 
*                                DO statements with CONTINUE statements
*  2009 August 19 B.E. Stetzler  Capitalized all variables for FORTRAN
*                                77 compatibility
*-----------------------------------------------------------------------

      IMPLICIT NONE
      INTEGER I,J,N,NF
      REAL X,HC,OM
      DOUBLE PRECISION SC,SCR
      DIMENSION X(*),HC(*),SCR(*),OM(*)

*  Set up for start of recursion by computing harmonic values
*  at starting point and just before it

      DO I = 1,NF
         SCR(3*I-2) = HC(2*I-1)
         SCR(3*I-1) = HC(2*I-1)*COS(OM(I)) -HC(2*I)*SIN(OM(I))
         SCR(3*I) = 2.*DCOS(DBLE(OM(I)))
      ENDDO

*  Do recursion over data
      DO I = 1,N
         X(I) = 0.
*  Then do recursive computation for each harmonic
         DO J  = 1,NF
            X(I) = X(I) + SCR(3*J-2)
            SC = SCR(3*J-2)
            SCR(3*J-2) = SCR(3*J)*SC-SCR(3*J-1)
            SCR(3*J-1) = SC
         ENDDO
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
