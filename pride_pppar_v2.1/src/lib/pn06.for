      SUBROUTINE iau_PN06 ( DATE1, DATE2, DPSI, DEPS,
C     :                      EPSA, RB, RP, RBP, RN, RBPN )
     :                      RBPN )
*+
*  - - - - - - - - -
*   i a u _ P N 0 6
*  - - - - - - - - -
*
*  Precession-nutation, IAU 2006 model:  a multi-purpose routine,
*  supporting classical (equinox-based) use directly and CIO-based use
*  indirectly.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  support routine.
*
*  Given:
*     DATE1,DATE2   d       TT as a 2-part Julian Date (Note 1)
*     DPSI,DEPS     d       nutation (Note 2)
*
*  Returned:
*     EPSA          d       mean obliquity (Note 3)
*     RB          d(3,3)    frame bias matrix (Note 4)
*     RP          d(3,3)    precession matrix (Note 5)
*     RBP         d(3,3)    bias-precession matrix (Note 6)
*     RN          d(3,3)    nutation matrix (Note 7)
*     RBPN        d(3,3)    GCRS-to-true matrix (Note 8)
*
*  Notes:
*
*  1) The TT date DATE1+DATE2 is a Julian Date, apportioned in any
*     convenient way between the two arguments.  For example,
*     JD(TT)=2450123.7 could be expressed in any of these ways,
*     among others:
*
*            DATE1          DATE2
*
*         2450123.7D0        0D0        (JD method)
*          2451545D0      -1421.3D0     (J2000 method)
*         2400000.5D0     50123.2D0     (MJD method)
*         2450123.5D0       0.2D0       (date & time method)
*
*     The JD method is the most natural and convenient to use in
*     cases where the loss of several decimal digits of resolution
*     is acceptable.  The J2000 method is best matched to the way
*     the argument is handled internally and will deliver the
*     optimum resolution.  The MJD method and the date & time methods
*     are both good compromises between resolution and convenience.
*
*  2) The caller is responsible for providing the nutation components;
*     they are in longitude and obliquity, in radians and are with
*     respect to the equinox and ecliptic of date.  For high-accuracy
*     applications, free core nutation should be included as well as
*     any other relevant corrections to the position of the CIP.
*
*  3) The returned mean obliquity is consistent with the IAU 2006
*     precession.
*
*  4) The matrix RB transforms vectors from GCRS to mean J2000.0 by
*     applying frame bias.
*
*  5) The matrix RP transforms vectors from mean J2000.0 to mean of date
*     by applying precession.
*
*  6) The matrix RBP transforms vectors from GCRS to mean of date by
*     applying frame bias then precession.  It is the product RP x RB.
*
*  7) The matrix RN transforms vectors from mean of date to true of date
*     by applying the nutation (luni-solar + planetary).
*
*  8) The matrix RBPN transforms vectors from GCRS to true of date
*     (CIP/equinox).  It is the product RN x RBP, applying frame bias,
*     precession and nutation in that order.
*
*  9) The X,Y,Z coordinates of the Celestial Intermediate Pole are
*     elements (3,1-3) of the matrix RBPN.
*
*  Called:
*     iau_PFW06    bias-precession F-W angles, IAU 2006
*     iau_FW2M     F-W angles to r-matrix
*     iau_TR       transpose r-matrix
*     iau_RXR      product of two r-matrices
*
*  References:
*
*     Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855
*
*     Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981
*
*  This revision:  2013 November 14
*
*  SOFA release 2016-05-03
*
*  Copyright (C) 2016 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION DATE1, DATE2, DPSI, DEPS,
     :                 EPSA, RB(3,3), RP(3,3), RBP(3,3),
     :                 RN(3,3), RBPN(3,3)

*  JD for MJD 0
      DOUBLE PRECISION DJM0
      PARAMETER (DJM0 = 2400000.5D0 )

*  Reference epoch (J2000.0), MJD
      DOUBLE PRECISION DJM00
      PARAMETER ( DJM00 = 51544.5D0 )

      DOUBLE PRECISION GAMB, PHIB, PSIB, EPS, RT(3,3)

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Bias-precession Fukushima-Williams angles of J2000.0 = frame bias.
      CALL iau_PFW06 ( DJM0, DJM00, GAMB, PHIB, PSIB, EPS )

*  B matrix.
      CALL iau_FW2M ( GAMB, PHIB, PSIB, EPS, RB )

*  Bias-precession Fukushima-Williams angles of date.
      CALL iau_PFW06 ( DATE1, DATE2, GAMB, PHIB, PSIB, EPS )

*  Bias-precession matrix.
      CALL iau_FW2M ( GAMB, PHIB, PSIB, EPS, RBP )

*  Solve for precession matrix.
      CALL iau_TR ( RB, RT )
      CALL iau_RXR ( RBP, RT, RP )

*  Equinox-based bias-precession-nutation matrix.
      CALL iau_FW2M ( GAMB, PHIB, PSIB+DPSI, EPS+DEPS, RBPN )

*  Solve for nutation matrix.
      CALL iau_TR ( RBP, RT )
      CALL iau_RXR ( RBPN, RT, RN )

*  Obliquity, mean of date.
      EPSA = EPS

*  Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2016
*  Standards Of Fundamental Astronomy Board
*  of the International Astronomical Union.
*
*  =====================
*  SOFA Software License
*  =====================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING SIX TERMS AND
*  CONDITIONS WHICH APPLY TO ITS USE.
*
*  1. The Software is owned by the IAU SOFA Board ("SOFA").
*
*  2. Permission is granted to anyone to use the SOFA software for any
*     purpose, including commercial applications, free of charge and
*     without payment of royalties, subject to the conditions and
*     restrictions listed below.
*
*  3. You (the user) may copy and distribute SOFA source code to others,
*     and use and adapt its code and algorithms in your own software,
*     on a world-wide, royalty-free basis.  That portion of your
*     distribution that does not consist of intact and unchanged copies
*     of SOFA source code files is a "derived work" that must comply
*     with the following requirements:
*
*     a) Your work shall be marked or carry a statement that it
*        (i) uses routines and computations derived by you from
*        software provided by SOFA under license to you; and
*        (ii) does not itself constitute software provided by and/or
*        endorsed by SOFA.
*
*     b) The source code of your derived work must contain descriptions
*        of how the derived work is based upon, contains and/or differs
*        from the original SOFA software.
*
*     c) The names of all routines in your derived work shall not
*        include the prefix "iau" or "sofa" or trivial modifications
*        thereof such as changes of case.
*
*     d) The origin of the SOFA components of your derived work must
*        not be misrepresented;  you must not claim that you wrote the
*        original software, nor file a patent application for SOFA
*        software or algorithms embedded in the SOFA software.
*
*     e) These requirements must be reproduced intact in any source
*        distribution and shall apply to anyone to whom you have
*        granted a further right to modify the source code of your
*        derived work.
*
*     Note that, as originally distributed, the SOFA software is
*     intended to be a definitive implementation of the IAU standards,
*     and consequently third-party modifications are discouraged.  All
*     variations, no matter how minor, must be explicitly marked as
*     such, as explained above.
*
*  4. You shall not cause the SOFA software to be brought into
*     disrepute, either by misuse, or use for inappropriate tasks, or
*     by inappropriate modification.
*
*  5. The SOFA software is provided "as is" and SOFA makes no warranty
*     as to its use or performance.   SOFA does not and cannot warrant
*     the performance or results which the user may obtain by using the
*     SOFA software.  SOFA makes no warranties, express or implied, as
*     to non-infringement of third party rights, merchantability, or
*     fitness for any particular purpose.  In no event will SOFA be
*     liable to the user for any consequential, incidental, or special
*     damages, including any lost profits or lost savings, even if a
*     SOFA representative has been advised of such damages, or for any
*     claim by any third party.
*
*  6. The provision of any version of the SOFA software under the terms
*     and conditions specified herein does not imply that future
*     versions will also be made available under the same terms and
*     conditions.
*
*  In any published work or commercial product which uses the SOFA
*  software directly, acknowledgement (see www.iausofa.org) is
*  appreciated.
*
*  Correspondence concerning SOFA software should be addressed as
*  follows:
*
*      By email:  sofa@ukho.gov.uk
*      By post:   IAU SOFA Center
*                 HM Nautical Almanac Office
*                 UK Hydrographic Office
*                 Admiralty Way, Taunton
*                 Somerset, TA1 2DN
*                 United Kingdom
*
*-----------------------------------------------------------------------

      END
