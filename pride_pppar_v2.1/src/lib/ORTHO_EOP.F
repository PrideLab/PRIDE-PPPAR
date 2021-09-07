      SUBROUTINE ORTHO_EOP ( TIME, EOP )
*+
*  - - - - - - - - - -
*   O R T H O _ E O P 
*  - - - - - - - - - -
*
*  This routine is part of the International Earth Rotation and
*  Reference Systems Service (IERS) Conventions software collection.
*
*  The purpose of the subroutine is to compute the diurnal and semi-
*  diurnal variations in Earth Orientation Parameters (x,y, UT1) from
*  ocean tides.
* 
*  In general, Class 1, 2, and 3 models represent physical effects that
*  act on geodetic parameters while canonical models provide lower-level
*  representations or basic computations that are used by Class 1, 2, or
*  3 models.
* 
*  Status: Class 1 model
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
*     TIME           d     Modified Julian Date 
*
*  Returned:
*     EOP            d     delta_x, in microarcseconds
*                          delta_y, in microarcseconds
*                          delta_UT1, in microseconds 
*
*  Notes:
*
*  1) The diurnal and semidiurnal orthoweights fit to the 8 constituents
*     are listed in Reference 1.
*     
*  Called:
*     CNMTX                Compute time dependent part of second degree
*                          diurnal and semidiurnal tidal potential from
*                          the dominant spectral lines in the Cartwright-
*                          Tayler-Edden harmonic decomposition                  *            
*  Test case:
*     given input: MJD = 47100D0
*
*     expected output: delta_x = -162.8386373279636530D0 microarcseconds
*                      delta_y = 117.7907525842668974D0 microarcseconds
*                      delta_UT1 = -23.39092370609808214D0 microseconds
*
*  References:
*
*     Ray, R. D., Steinberg, D. J., Chao, B. F., and Cartwright, D. E.,
*     "Diurnal and Semidiurnal Variations in the Earth's Rotation
*     Rate Induced by Ocean Tides", 1994, Science, 264, pp. 830-832
*
*     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
*     IERS Technical Note No. 36, BKG (2010)
*
*  Revisions:
*  1997 February    R. Eanes         Original code
*  2008 November 07 B.E. Stetzler    Added header and copyright
*  2008 November 07 B.E. Stetzler    Provided test case
*  2009 May      12 B.E. Stetzler    Replaced ENDDO statements with
*                                    CONTINUE statements
*  2009 June     09 B.E. Stetzler    Updated validated test case values
*                                    based on changes to CNMTX.F
*  2010 March    19 B.E. Stetzler    Capitalized variables for FORTRAN
*                                    77 backwards compatibility
*-----------------------------------------------------------------------

      IMPLICIT NONE
      INTEGER K, J
      DOUBLE PRECISION EOP(3), ORTHOW(12,3), TIME, H(12)
      DATA ORTHOW /
     .-6.77832D0,-14.86323D0,0.47884D0,-1.45303D0,0.16406D0,  0.42030D0,
     . 0.09398D0,25.73054D0,-4.77974D0, 0.28080D0,1.94539D0, -0.73089D0,
     .14.86283D0,-6.77846D0, 1.45234D0, 0.47888D0,-0.42056D0, 0.16469D0,
     .15.30276D0,-4.30615D0, 0.07564D0, 2.28321D0,-0.45717D0,-1.62010D0,
     .-1.76335D0,1.03364D0, -0.27553D0, 0.34569D0,-0.12343D0,-0.10146D0,
     .-0.47119D0,1.28997D0, -0.19336D0, 0.02724D0, 0.08955D0, 0.04726D0/

* Compute the partials of the tidal variations to the orthoweights
      CALL CNMTX (TIME, H)

* Compute eop changes
      DO 20 K=1,3
         EOP(K) = 0D0
         DO 40 J=1,12
            EOP(K) = EOP(K) + H(J)*ORTHOW(J,K)
40       CONTINUE         

20    CONTINUE
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
