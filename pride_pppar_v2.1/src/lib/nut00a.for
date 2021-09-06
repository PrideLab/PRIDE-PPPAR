      SUBROUTINE iau_NUT00A ( DATE1, DATE2, DPSI, DEPS )
*+
*  - - - - - - - - - - -
*   i a u _ N U T 0 0 A
*  - - - - - - - - - - -
*
*  Nutation, IAU 2000A model (MHB2000 luni-solar and planetary nutation
*  with free core nutation omitted).
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  canonical model.
*
*  Given:
*     DATE1,DATE2   d    TT as a 2-part Julian Date (Note 1)
*
*  Returned:
*     DPSI,DEPS     d    nutation, luni-solar + planetary (Note 2)
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
*  2) The nutation components in longitude and obliquity are in radians
*     and with respect to the equinox and ecliptic of date.  The
*     obliquity at J2000.0 is assumed to be the Lieske et al. (1977)
*     value of 84381.448 arcsec.
*
*     Both the luni-solar and planetary nutations are included.  The
*     latter are due to direct planetary nutations and the perturbations
*     of the lunar and terrestrial orbits.
*
*  3) The routine computes the MHB2000 nutation series with the
*     associated corrections for planetary nutations.  It is an
*     implementation of the nutation part of the IAU 2000A precession-
*     nutation model, formally adopted by the IAU General Assembly in
*     2000, namely MHB2000 (Mathews et al. 2002), but with the free core
*     nutation (FCN - see Note 4) omitted.
*
*  4) The full MHB2000 model also contains contributions to the
*     nutations in longitude and obliquity due to the free-excitation of
*     the free-core-nutation during the period 1979-2000.  These FCN
*     terms, which are time-dependent and unpredictable, are NOT
*     included in the present routine and, if required, must be
*     independently computed.  With the FCN corrections included, the
*     present routine delivers a pole which is at current epochs
*     accurate to a few hundred microarcseconds.  The omission of FCN
*     introduces further errors of about that size.
*
*  5) The present routine provides classical nutation.  The MHB2000
*     algorithm, from which it is adapted, deals also with (i) the
*     offsets between the GCRS and mean poles and (ii) the adjustments
*     in longitude and obliquity due to the changed precession rates.
*     These additional functions, namely frame bias and precession
*     adjustments, are supported by the SOFA routines iau_BI00 and
*     iau_PR00.
*
*  6) The MHB2000 algorithm also provides "total" nutations, comprising
*     the arithmetic sum of the frame bias, precession adjustments,
*     luni-solar nutation and planetary nutation.  These total nutations
*     can be used in combination with an existing IAU 1976 precession
*     implementation, such as iau_PMAT76, to deliver GCRS-to-true
*     predictions of sub-mas accuracy at current epochs.  However, there
*     are three shortcomings in the MHB2000 model that must be taken
*     into account if more accurate or definitive results are required
*     (see Wallace 2002):
*
*       (i) The MHB2000 total nutations are simply arithmetic sums,
*           yet in reality the various components are successive Euler
*           rotations.  This slight lack of rigor leads to cross terms
*           that exceed 1 mas after a century.  The rigorous procedure
*           is to form the GCRS-to-true rotation matrix by applying the
*           bias, precession and nutation in that order.
*
*      (ii) Although the precession adjustments are stated to be with
*           respect to Lieske et al. (1977), the MHB2000 model does
*           not specify which set of Euler angles are to be used and
*           how the adjustments are to be applied.  The most literal and
*           straightforward procedure is to adopt the 4-rotation
*           epsilon_0, psi_A, omega_A, xi_A option, and to add DPSIPR to
*           psi_A and DEPSPR to both omega_A and eps_A.
*
*     (iii) The MHB2000 model predates the determination by Chapront
*           et al. (2002) of a 14.6 mas displacement between the J2000.0
*           mean equinox and the origin of the ICRS frame.  It should,
*           however, be noted that neglecting this displacement when
*           calculating star coordinates does not lead to a 14.6 mas
*           change in right ascension, only a small second-order
*           distortion in the pattern of the precession-nutation effect.
*
*     For these reasons, the SOFA routines do not generate the "total
*     nutations" directly, though they can of course easily be generated
*     by calling iau_BI00, iau_PR00 and the present routine and adding
*     the results.
*
*  7) The MHB2000 model contains 41 instances where the same frequency
*     appears multiple times, of which 38 are duplicates and three are
*     triplicates.  To keep the present code close to the original MHB
*     algorithm, this small inefficiency has not been corrected.
*
*  Called:
*     iau_FAL03    mean anomaly of the Moon
*     iau_FAF03    mean argument of the latitude of the Moon
*     iau_FAOM03   mean longitude of the Moon's ascending node
*     iau_FAME03   mean longitude of Mercury
*     iau_FAVE03   mean longitude of Venus
*     iau_FAE03    mean longitude of Earth
*     iau_FAMA03   mean longitude of Mars
*     iau_FAJU03   mean longitude of Jupiter
*     iau_FASA03   mean longitude of Saturn
*     iau_FAUR03   mean longitude of Uranus
*     iau_FAPA03   general accumulated precession in longitude
*
*  References:
*
*     Chapront, J., Chapront-Touze, M. & Francou, G. 2002,
*     Astron.Astrophys. 387, 700
*
*     Lieske, J.H., Lederle, T., Fricke, W. & Morando, B. 1977,
*     Astron.Astrophys. 58, 1-16
*
*     Mathews, P.M., Herring, T.A., Buffet, B.A. 2002, J.Geophys.Res.
*     107, B4.  The MHB_2000 code itself was obtained on 9th September
*     2002 from ftp//maia.usno.navy.mil/conv2000/chapter5/IAU2000A.
*
*     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
*     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
*
*     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
*     Astron.Astrophys.Supp.Ser. 135, 111
*
*     Wallace, P.T., "Software for Implementing the IAU 2000
*     Resolutions", in IERS Workshop 5.1 (2002)
*
*  This revision:  2009 December 15
*
*  SOFA release 2016-05-03
*
*  Copyright (C) 2016 IAU SOFA Board.  See notes at end.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION DATE1, DATE2, DPSI, DEPS

*  Arcseconds to radians
      DOUBLE PRECISION DAS2R
      PARAMETER ( DAS2R = 4.848136811095359935899141D-6 )

*  Arcseconds in a full circle
      DOUBLE PRECISION TURNAS
      PARAMETER ( TURNAS = 1296000D0 )

*  2Pi
      DOUBLE PRECISION D2PI
      PARAMETER ( D2PI = 6.283185307179586476925287D0 )

*  Units of 0.1 microarcsecond to radians
      DOUBLE PRECISION U2R
      PARAMETER ( U2R = DAS2R/1D7 )

*  Reference epoch (J2000.0), JD
      DOUBLE PRECISION DJ00
      PARAMETER ( DJ00 = 2451545D0 )

*  Days per Julian century
      DOUBLE PRECISION DJC
      PARAMETER ( DJC = 36525D0 )

*  Miscellaneous
      INTEGER I, J
      DOUBLE PRECISION T, EL, ELP, F, D, OM, ARG, DP, DE, SARG, CARG,
     :                 AL, ALSU, AF, AD, AOM, ALME, ALVE, ALEA, ALMA,
     :                 ALJU, ALSA, ALUR, ALNE, APA, DPSILS, DEPSLS,
     :                 DPSIPL, DEPSPL
      DOUBLE PRECISION iau_FAL03, iau_FAF03, iau_FAOM03, iau_FAME03,
     :                 iau_FAVE03, iau_FAE03, iau_FAMA03, iau_FAJU03,
     :                 iau_FASA03, iau_FAUR03, iau_FAPA03

*  -------------------------
*  Luni-Solar nutation model
*  -------------------------

*  Number of terms in the luni-solar nutation model
      INTEGER NLS
      PARAMETER ( NLS = 678 )

*  Coefficients for fundamental arguments
      INTEGER NALS(5,NLS)

*  Longitude and obliquity coefficients
      DOUBLE PRECISION CLS(6,NLS)

*  ---------------
*  Planetary terms
*  ---------------

*  Number of terms in the planetary nutation model
      INTEGER NPL
      PARAMETER ( NPL = 687 )

*  Coefficients for fundamental arguments
      INTEGER NAPL(14,NPL)

*  Longitude and obliquity coefficients
      INTEGER ICPL(4,NPL)

*  ----------------------------------------
*  Tables of argument and term coefficients
*  ----------------------------------------

*
*  Luni-Solar argument multipliers
*               L     L'    F     D     Om

      DATA ( ( NALS(I,J), I=1,5 ), J=  1, 10 ) /
     :          0,    0,    0,    0,    1,
     :          0,    0,    2,   -2,    2,
     :          0,    0,    2,    0,    2,
     :          0,    0,    0,    0,    2,
     :          0,    1,    0,    0,    0,
     :          0,    1,    2,   -2,    2,
     :          1,    0,    0,    0,    0,
     :          0,    0,    2,    0,    1,
     :          1,    0,    2,    0,    2,
     :          0,   -1,    2,   -2,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J= 11, 20 ) /
     :          0,    0,    2,   -2,    1,
     :         -1,    0,    2,    0,    2,
     :         -1,    0,    0,    2,    0,
     :          1,    0,    0,    0,    1,
     :         -1,    0,    0,    0,    1,
     :         -1,    0,    2,    2,    2,
     :          1,    0,    2,    0,    1,
     :         -2,    0,    2,    0,    1,
     :          0,    0,    0,    2,    0,
     :          0,    0,    2,    2,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J= 21, 30 ) /
     :          0,   -2,    2,   -2,    2,
     :         -2,    0,    0,    2,    0,
     :          2,    0,    2,    0,    2,
     :          1,    0,    2,   -2,    2,
     :         -1,    0,    2,    0,    1,
     :          2,    0,    0,    0,    0,
     :          0,    0,    2,    0,    0,
     :          0,    1,    0,    0,    1,
     :         -1,    0,    0,    2,    1,
     :          0,    2,    2,   -2,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J= 31, 40 ) /
     :          0,    0,   -2,    2,    0,
     :          1,    0,    0,   -2,    1,
     :          0,   -1,    0,    0,    1,
     :         -1,    0,    2,    2,    1,
     :          0,    2,    0,    0,    0,
     :          1,    0,    2,    2,    2,
     :         -2,    0,    2,    0,    0,
     :          0,    1,    2,    0,    2,
     :          0,    0,    2,    2,    1,
     :          0,   -1,    2,    0,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J= 41, 50 ) /
     :          0,    0,    0,    2,    1,
     :          1,    0,    2,   -2,    1,
     :          2,    0,    2,   -2,    2,
     :         -2,    0,    0,    2,    1,
     :          2,    0,    2,    0,    1,
     :          0,   -1,    2,   -2,    1,
     :          0,    0,    0,   -2,    1,
     :         -1,   -1,    0,    2,    0,
     :          2,    0,    0,   -2,    1,
     :          1,    0,    0,    2,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J= 51, 60 ) /
     :          0,    1,    2,   -2,    1,
     :          1,   -1,    0,    0,    0,
     :         -2,    0,    2,    0,    2,
     :          3,    0,    2,    0,    2,
     :          0,   -1,    0,    2,    0,
     :          1,   -1,    2,    0,    2,
     :          0,    0,    0,    1,    0,
     :         -1,   -1,    2,    2,    2,
     :         -1,    0,    2,    0,    0,
     :          0,   -1,    2,    2,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J= 61, 70 ) /
     :         -2,    0,    0,    0,    1,
     :          1,    1,    2,    0,    2,
     :          2,    0,    0,    0,    1,
     :         -1,    1,    0,    1,    0,
     :          1,    1,    0,    0,    0,
     :          1,    0,    2,    0,    0,
     :         -1,    0,    2,   -2,    1,
     :          1,    0,    0,    0,    2,
     :         -1,    0,    0,    1,    0,
     :          0,    0,    2,    1,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J= 71, 80 ) /
     :         -1,    0,    2,    4,    2,
     :         -1,    1,    0,    1,    1,
     :          0,   -2,    2,   -2,    1,
     :          1,    0,    2,    2,    1,
     :         -2,    0,    2,    2,    2,
     :         -1,    0,    0,    0,    2,
     :          1,    1,    2,   -2,    2,
     :         -2,    0,    2,    4,    2,
     :         -1,    0,    4,    0,    2,
     :          2,    0,    2,   -2,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J= 81, 90 ) /
     :          2,    0,    2,    2,    2,
     :          1,    0,    0,    2,    1,
     :          3,    0,    0,    0,    0,
     :          3,    0,    2,   -2,    2,
     :          0,    0,    4,   -2,    2,
     :          0,    1,    2,    0,    1,
     :          0,    0,   -2,    2,    1,
     :          0,    0,    2,   -2,    3,
     :         -1,    0,    0,    4,    0,
     :          2,    0,   -2,    0,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J= 91,100 ) /
     :         -2,    0,    0,    4,    0,
     :         -1,   -1,    0,    2,    1,
     :         -1,    0,    0,    1,    1,
     :          0,    1,    0,    0,    2,
     :          0,    0,   -2,    0,    1,
     :          0,   -1,    2,    0,    1,
     :          0,    0,    2,   -1,    2,
     :          0,    0,    2,    4,    2,
     :         -2,   -1,    0,    2,    0,
     :          1,    1,    0,   -2,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=101,110 ) /
     :         -1,    1,    0,    2,    0,
     :         -1,    1,    0,    1,    2,
     :          1,   -1,    0,    0,    1,
     :          1,   -1,    2,    2,    2,
     :         -1,    1,    2,    2,    2,
     :          3,    0,    2,    0,    1,
     :          0,    1,   -2,    2,    0,
     :         -1,    0,    0,   -2,    1,
     :          0,    1,    2,    2,    2,
     :         -1,   -1,    2,    2,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=111,120 ) /
     :          0,   -1,    0,    0,    2,
     :          1,    0,    2,   -4,    1,
     :         -1,    0,   -2,    2,    0,
     :          0,   -1,    2,    2,    1,
     :          2,   -1,    2,    0,    2,
     :          0,    0,    0,    2,    2,
     :          1,   -1,    2,    0,    1,
     :         -1,    1,    2,    0,    2,
     :          0,    1,    0,    2,    0,
     :          0,   -1,   -2,    2,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=121,130 ) /
     :          0,    3,    2,   -2,    2,
     :          0,    0,    0,    1,    1,
     :         -1,    0,    2,    2,    0,
     :          2,    1,    2,    0,    2,
     :          1,    1,    0,    0,    1,
     :          1,    1,    2,    0,    1,
     :          2,    0,    0,    2,    0,
     :          1,    0,   -2,    2,    0,
     :         -1,    0,    0,    2,    2,
     :          0,    1,    0,    1,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=131,140 ) /
     :          0,    1,    0,   -2,    1,
     :         -1,    0,    2,   -2,    2,
     :          0,    0,    0,   -1,    1,
     :         -1,    1,    0,    0,    1,
     :          1,    0,    2,   -1,    2,
     :          1,   -1,    0,    2,    0,
     :          0,    0,    0,    4,    0,
     :          1,    0,    2,    1,    2,
     :          0,    0,    2,    1,    1,
     :          1,    0,    0,   -2,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=141,150 ) /
     :         -1,    0,    2,    4,    1,
     :          1,    0,   -2,    0,    1,
     :          1,    1,    2,   -2,    1,
     :          0,    0,    2,    2,    0,
     :         -1,    0,    2,   -1,    1,
     :         -2,    0,    2,    2,    1,
     :          4,    0,    2,    0,    2,
     :          2,   -1,    0,    0,    0,
     :          2,    1,    2,   -2,    2,
     :          0,    1,    2,    1,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=151,160 ) /
     :          1,    0,    4,   -2,    2,
     :         -1,   -1,    0,    0,    1,
     :          0,    1,    0,    2,    1,
     :         -2,    0,    2,    4,    1,
     :          2,    0,    2,    0,    0,
     :          1,    0,    0,    1,    0,
     :         -1,    0,    0,    4,    1,
     :         -1,    0,    4,    0,    1,
     :          2,    0,    2,    2,    1,
     :          0,    0,    2,   -3,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=161,170 ) /
     :         -1,   -2,    0,    2,    0,
     :          2,    1,    0,    0,    0,
     :          0,    0,    4,    0,    2,
     :          0,    0,    0,    0,    3,
     :          0,    3,    0,    0,    0,
     :          0,    0,    2,   -4,    1,
     :          0,   -1,    0,    2,    1,
     :          0,    0,    0,    4,    1,
     :         -1,   -1,    2,    4,    2,
     :          1,    0,    2,    4,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=171,180 ) /
     :         -2,    2,    0,    2,    0,
     :         -2,   -1,    2,    0,    1,
     :         -2,    0,    0,    2,    2,
     :         -1,   -1,    2,    0,    2,
     :          0,    0,    4,   -2,    1,
     :          3,    0,    2,   -2,    1,
     :         -2,   -1,    0,    2,    1,
     :          1,    0,    0,   -1,    1,
     :          0,   -2,    0,    2,    0,
     :         -2,    0,    0,    4,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=181,190 ) /
     :         -3,    0,    0,    0,    1,
     :          1,    1,    2,    2,    2,
     :          0,    0,    2,    4,    1,
     :          3,    0,    2,    2,    2,
     :         -1,    1,    2,   -2,    1,
     :          2,    0,    0,   -4,    1,
     :          0,    0,    0,   -2,    2,
     :          2,    0,    2,   -4,    1,
     :         -1,    1,    0,    2,    1,
     :          0,    0,    2,   -1,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=191,200 ) /
     :          0,   -2,    2,    2,    2,
     :          2,    0,    0,    2,    1,
     :          4,    0,    2,   -2,    2,
     :          2,    0,    0,   -2,    2,
     :          0,    2,    0,    0,    1,
     :          1,    0,    0,   -4,    1,
     :          0,    2,    2,   -2,    1,
     :         -3,    0,    0,    4,    0,
     :         -1,    1,    2,    0,    1,
     :         -1,   -1,    0,    4,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=201,210 ) /
     :         -1,   -2,    2,    2,    2,
     :         -2,   -1,    2,    4,    2,
     :          1,   -1,    2,    2,    1,
     :         -2,    1,    0,    2,    0,
     :         -2,    1,    2,    0,    1,
     :          2,    1,    0,   -2,    1,
     :         -3,    0,    2,    0,    1,
     :         -2,    0,    2,   -2,    1,
     :         -1,    1,    0,    2,    2,
     :          0,   -1,    2,   -1,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=211,220 ) /
     :         -1,    0,    4,   -2,    2,
     :          0,   -2,    2,    0,    2,
     :         -1,    0,    2,    1,    2,
     :          2,    0,    0,    0,    2,
     :          0,    0,    2,    0,    3,
     :         -2,    0,    4,    0,    2,
     :         -1,    0,   -2,    0,    1,
     :         -1,    1,    2,    2,    1,
     :          3,    0,    0,    0,    1,
     :         -1,    0,    2,    3,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=221,230 ) /
     :          2,   -1,    2,    0,    1,
     :          0,    1,    2,    2,    1,
     :          0,   -1,    2,    4,    2,
     :          2,   -1,    2,    2,    2,
     :          0,    2,   -2,    2,    0,
     :         -1,   -1,    2,   -1,    1,
     :          0,   -2,    0,    0,    1,
     :          1,    0,    2,   -4,    2,
     :          1,   -1,    0,   -2,    1,
     :         -1,   -1,    2,    0,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=231,240 ) /
     :          1,   -1,    2,   -2,    2,
     :         -2,   -1,    0,    4,    0,
     :         -1,    0,    0,    3,    0,
     :         -2,   -1,    2,    2,    2,
     :          0,    2,    2,    0,    2,
     :          1,    1,    0,    2,    0,
     :          2,    0,    2,   -1,    2,
     :          1,    0,    2,    1,    1,
     :          4,    0,    0,    0,    0,
     :          2,    1,    2,    0,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=241,250 ) /
     :          3,   -1,    2,    0,    2,
     :         -2,    2,    0,    2,    1,
     :          1,    0,    2,   -3,    1,
     :          1,    1,    2,   -4,    1,
     :         -1,   -1,    2,   -2,    1,
     :          0,   -1,    0,   -1,    1,
     :          0,   -1,    0,   -2,    1,
     :         -2,    0,    0,    0,    2,
     :         -2,    0,   -2,    2,    0,
     :         -1,    0,   -2,    4,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=251,260 ) /
     :          1,   -2,    0,    0,    0,
     :          0,    1,    0,    1,    1,
     :         -1,    2,    0,    2,    0,
     :          1,   -1,    2,   -2,    1,
     :          1,    2,    2,   -2,    2,
     :          2,   -1,    2,   -2,    2,
     :          1,    0,    2,   -1,    1,
     :          2,    1,    2,   -2,    1,
     :         -2,    0,    0,   -2,    1,
     :          1,   -2,    2,    0,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=261,270 ) /
     :          0,    1,    2,    1,    1,
     :          1,    0,    4,   -2,    1,
     :         -2,    0,    4,    2,    2,
     :          1,    1,    2,    1,    2,
     :          1,    0,    0,    4,    0,
     :          1,    0,    2,    2,    0,
     :          2,    0,    2,    1,    2,
     :          3,    1,    2,    0,    2,
     :          4,    0,    2,    0,    1,
     :         -2,   -1,    2,    0,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=271,280 ) /
     :          0,    1,   -2,    2,    1,
     :          1,    0,   -2,    1,    0,
     :          0,   -1,   -2,    2,    1,
     :          2,   -1,    0,   -2,    1,
     :         -1,    0,    2,   -1,    2,
     :          1,    0,    2,   -3,    2,
     :          0,    1,    2,   -2,    3,
     :          0,    0,    2,   -3,    1,
     :         -1,    0,   -2,    2,    1,
     :          0,    0,    2,   -4,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=281,290 ) /
     :         -2,    1,    0,    0,    1,
     :         -1,    0,    0,   -1,    1,
     :          2,    0,    2,   -4,    2,
     :          0,    0,    4,   -4,    4,
     :          0,    0,    4,   -4,    2,
     :         -1,   -2,    0,    2,    1,
     :         -2,    0,    0,    3,    0,
     :          1,    0,   -2,    2,    1,
     :         -3,    0,    2,    2,    2,
     :         -3,    0,    2,    2,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=291,300 ) /
     :         -2,    0,    2,    2,    0,
     :          2,   -1,    0,    0,    1,
     :         -2,    1,    2,    2,    2,
     :          1,    1,    0,    1,    0,
     :          0,    1,    4,   -2,    2,
     :         -1,    1,    0,   -2,    1,
     :          0,    0,    0,   -4,    1,
     :          1,   -1,    0,    2,    1,
     :          1,    1,    0,    2,    1,
     :         -1,    2,    2,    2,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=301,310 ) /
     :          3,    1,    2,   -2,    2,
     :          0,   -1,    0,    4,    0,
     :          2,   -1,    0,    2,    0,
     :          0,    0,    4,    0,    1,
     :          2,    0,    4,   -2,    2,
     :         -1,   -1,    2,    4,    1,
     :          1,    0,    0,    4,    1,
     :          1,   -2,    2,    2,    2,
     :          0,    0,    2,    3,    2,
     :         -1,    1,    2,    4,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=311,320 ) /
     :          3,    0,    0,    2,    0,
     :         -1,    0,    4,    2,    2,
     :          1,    1,    2,    2,    1,
     :         -2,    0,    2,    6,    2,
     :          2,    1,    2,    2,    2,
     :         -1,    0,    2,    6,    2,
     :          1,    0,    2,    4,    1,
     :          2,    0,    2,    4,    2,
     :          1,    1,   -2,    1,    0,
     :         -3,    1,    2,    1,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=321,330 ) /
     :          2,    0,   -2,    0,    2,
     :         -1,    0,    0,    1,    2,
     :         -4,    0,    2,    2,    1,
     :         -1,   -1,    0,    1,    0,
     :          0,    0,   -2,    2,    2,
     :          1,    0,    0,   -1,    2,
     :          0,   -1,    2,   -2,    3,
     :         -2,    1,    2,    0,    0,
     :          0,    0,    2,   -2,    4,
     :         -2,   -2,    0,    2,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=331,340 ) /
     :         -2,    0,   -2,    4,    0,
     :          0,   -2,   -2,    2,    0,
     :          1,    2,    0,   -2,    1,
     :          3,    0,    0,   -4,    1,
     :         -1,    1,    2,   -2,    2,
     :          1,   -1,    2,   -4,    1,
     :          1,    1,    0,   -2,    2,
     :         -3,    0,    2,    0,    0,
     :         -3,    0,    2,    0,    2,
     :         -2,    0,    0,    1,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=341,350 ) /
     :          0,    0,   -2,    1,    0,
     :         -3,    0,    0,    2,    1,
     :         -1,   -1,   -2,    2,    0,
     :          0,    1,    2,   -4,    1,
     :          2,    1,    0,   -4,    1,
     :          0,    2,    0,   -2,    1,
     :          1,    0,    0,   -3,    1,
     :         -2,    0,    2,   -2,    2,
     :         -2,   -1,    0,    0,    1,
     :         -4,    0,    0,    2,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=351,360 ) /
     :          1,    1,    0,   -4,    1,
     :         -1,    0,    2,   -4,    1,
     :          0,    0,    4,   -4,    1,
     :          0,    3,    2,   -2,    2,
     :         -3,   -1,    0,    4,    0,
     :         -3,    0,    0,    4,    1,
     :          1,   -1,   -2,    2,    0,
     :         -1,   -1,    0,    2,    2,
     :          1,   -2,    0,    0,    1,
     :          1,   -1,    0,    0,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=361,370 ) /
     :          0,    0,    0,    1,    2,
     :         -1,   -1,    2,    0,    0,
     :          1,   -2,    2,   -2,    2,
     :          0,   -1,    2,   -1,    1,
     :         -1,    0,    2,    0,    3,
     :          1,    1,    0,    0,    2,
     :         -1,    1,    2,    0,    0,
     :          1,    2,    0,    0,    0,
     :         -1,    2,    2,    0,    2,
     :         -1,    0,    4,   -2,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=371,380 ) /
     :          3,    0,    2,   -4,    2,
     :          1,    2,    2,   -2,    1,
     :          1,    0,    4,   -4,    2,
     :         -2,   -1,    0,    4,    1,
     :          0,   -1,    0,    2,    2,
     :         -2,    1,    0,    4,    0,
     :         -2,   -1,    2,    2,    1,
     :          2,    0,   -2,    2,    0,
     :          1,    0,    0,    1,    1,
     :          0,    1,    0,    2,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=381,390 ) /
     :          1,   -1,    2,   -1,    2,
     :         -2,    0,    4,    0,    1,
     :          2,    1,    0,    0,    1,
     :          0,    1,    2,    0,    0,
     :          0,   -1,    4,   -2,    2,
     :          0,    0,    4,   -2,    4,
     :          0,    2,    2,    0,    1,
     :         -3,    0,    0,    6,    0,
     :         -1,   -1,    0,    4,    1,
     :          1,   -2,    0,    2,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=391,400 ) /
     :         -1,    0,    0,    4,    2,
     :         -1,   -2,    2,    2,    1,
     :         -1,    0,    0,   -2,    2,
     :          1,    0,   -2,   -2,    1,
     :          0,    0,   -2,   -2,    1,
     :         -2,    0,   -2,    0,    1,
     :          0,    0,    0,    3,    1,
     :          0,    0,    0,    3,    0,
     :         -1,    1,    0,    4,    0,
     :         -1,   -1,    2,    2,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=401,410 ) /
     :         -2,    0,    2,    3,    2,
     :          1,    0,    0,    2,    2,
     :          0,   -1,    2,    1,    2,
     :          3,   -1,    0,    0,    0,
     :          2,    0,    0,    1,    0,
     :          1,   -1,    2,    0,    0,
     :          0,    0,    2,    1,    0,
     :          1,    0,    2,    0,    3,
     :          3,    1,    0,    0,    0,
     :          3,   -1,    2,   -2,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=411,420 ) /
     :          2,    0,    2,   -1,    1,
     :          1,    1,    2,    0,    0,
     :          0,    0,    4,   -1,    2,
     :          1,    2,    2,    0,    2,
     :         -2,    0,    0,    6,    0,
     :          0,   -1,    0,    4,    1,
     :         -2,   -1,    2,    4,    1,
     :          0,   -2,    2,    2,    1,
     :          0,   -1,    2,    2,    0,
     :         -1,    0,    2,    3,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=421,430 ) /
     :         -2,    1,    2,    4,    2,
     :          2,    0,    0,    2,    2,
     :          2,   -2,    2,    0,    2,
     :         -1,    1,    2,    3,    2,
     :          3,    0,    2,   -1,    2,
     :          4,    0,    2,   -2,    1,
     :         -1,    0,    0,    6,    0,
     :         -1,   -2,    2,    4,    2,
     :         -3,    0,    2,    6,    2,
     :         -1,    0,    2,    4,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=431,440 ) /
     :          3,    0,    0,    2,    1,
     :          3,   -1,    2,    0,    1,
     :          3,    0,    2,    0,    0,
     :          1,    0,    4,    0,    2,
     :          5,    0,    2,   -2,    2,
     :          0,   -1,    2,    4,    1,
     :          2,   -1,    2,    2,    1,
     :          0,    1,    2,    4,    2,
     :          1,   -1,    2,    4,    2,
     :          3,   -1,    2,    2,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=441,450 ) /
     :          3,    0,    2,    2,    1,
     :          5,    0,    2,    0,    2,
     :          0,    0,    2,    6,    2,
     :          4,    0,    2,    2,    2,
     :          0,   -1,    1,   -1,    1,
     :         -1,    0,    1,    0,    3,
     :          0,   -2,    2,   -2,    3,
     :          1,    0,   -1,    0,    1,
     :          2,   -2,    0,   -2,    1,
     :         -1,    0,    1,    0,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=451,460 ) /
     :         -1,    0,    1,    0,    1,
     :         -1,   -1,    2,   -1,    2,
     :         -2,    2,    0,    2,    2,
     :         -1,    0,    1,    0,    0,
     :         -4,    1,    2,    2,    2,
     :         -3,    0,    2,    1,    1,
     :         -2,   -1,    2,    0,    2,
     :          1,    0,   -2,    1,    1,
     :          2,   -1,   -2,    0,    1,
     :         -4,    0,    2,    2,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=461,470 ) /
     :         -3,    1,    0,    3,    0,
     :         -1,    0,   -1,    2,    0,
     :          0,   -2,    0,    0,    2,
     :          0,   -2,    0,    0,    2,
     :         -3,    0,    0,    3,    0,
     :         -2,   -1,    0,    2,    2,
     :         -1,    0,   -2,    3,    0,
     :         -4,    0,    0,    4,    0,
     :          2,    1,   -2,    0,    1,
     :          2,   -1,    0,   -2,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=471,480 ) /
     :          0,    0,    1,   -1,    0,
     :         -1,    2,    0,    1,    0,
     :         -2,    1,    2,    0,    2,
     :          1,    1,    0,   -1,    1,
     :          1,    0,    1,   -2,    1,
     :          0,    2,    0,    0,    2,
     :          1,   -1,    2,   -3,    1,
     :         -1,    1,    2,   -1,    1,
     :         -2,    0,    4,   -2,    2,
     :         -2,    0,    4,   -2,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=481,490 ) /
     :         -2,   -2,    0,    2,    1,
     :         -2,    0,   -2,    4,    0,
     :          1,    2,    2,   -4,    1,
     :          1,    1,    2,   -4,    2,
     :         -1,    2,    2,   -2,    1,
     :          2,    0,    0,   -3,    1,
     :         -1,    2,    0,    0,    1,
     :          0,    0,    0,   -2,    0,
     :         -1,   -1,    2,   -2,    2,
     :         -1,    1,    0,    0,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=491,500 ) /
     :          0,    0,    0,   -1,    2,
     :         -2,    1,    0,    1,    0,
     :          1,   -2,    0,   -2,    1,
     :          1,    0,   -2,    0,    2,
     :         -3,    1,    0,    2,    0,
     :         -1,    1,   -2,    2,    0,
     :         -1,   -1,    0,    0,    2,
     :         -3,    0,    0,    2,    0,
     :         -3,   -1,    0,    2,    0,
     :          2,    0,    2,   -6,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=501,510 ) /
     :          0,    1,    2,   -4,    2,
     :          2,    0,    0,   -4,    2,
     :         -2,    1,    2,   -2,    1,
     :          0,   -1,    2,   -4,    1,
     :          0,    1,    0,   -2,    2,
     :         -1,    0,    0,   -2,    0,
     :          2,    0,   -2,   -2,    1,
     :         -4,    0,    2,    0,    1,
     :         -1,   -1,    0,   -1,    1,
     :          0,    0,   -2,    0,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=511,520 ) /
     :         -3,    0,    0,    1,    0,
     :         -1,    0,   -2,    1,    0,
     :         -2,    0,   -2,    2,    1,
     :          0,    0,   -4,    2,    0,
     :         -2,   -1,   -2,    2,    0,
     :          1,    0,    2,   -6,    1,
     :         -1,    0,    2,   -4,    2,
     :          1,    0,    0,   -4,    2,
     :          2,    1,    2,   -4,    2,
     :          2,    1,    2,   -4,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=521,530 ) /
     :          0,    1,    4,   -4,    4,
     :          0,    1,    4,   -4,    2,
     :         -1,   -1,   -2,    4,    0,
     :         -1,   -3,    0,    2,    0,
     :         -1,    0,   -2,    4,    1,
     :         -2,   -1,    0,    3,    0,
     :          0,    0,   -2,    3,    0,
     :         -2,    0,    0,    3,    1,
     :          0,   -1,    0,    1,    0,
     :         -3,    0,    2,    2,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=531,540 ) /
     :          1,    1,   -2,    2,    0,
     :         -1,    1,    0,    2,    2,
     :          1,   -2,    2,   -2,    1,
     :          0,    0,    1,    0,    2,
     :          0,    0,    1,    0,    1,
     :          0,    0,    1,    0,    0,
     :         -1,    2,    0,    2,    1,
     :          0,    0,    2,    0,    2,
     :         -2,    0,    2,    0,    2,
     :          2,    0,    0,   -1,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=541,550 ) /
     :          3,    0,    0,   -2,    1,
     :          1,    0,    2,   -2,    3,
     :          1,    2,    0,    0,    1,
     :          2,    0,    2,   -3,    2,
     :         -1,    1,    4,   -2,    2,
     :         -2,   -2,    0,    4,    0,
     :          0,   -3,    0,    2,    0,
     :          0,    0,   -2,    4,    0,
     :         -1,   -1,    0,    3,    0,
     :         -2,    0,    0,    4,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=551,560 ) /
     :         -1,    0,    0,    3,    1,
     :          2,   -2,    0,    0,    0,
     :          1,   -1,    0,    1,    0,
     :         -1,    0,    0,    2,    0,
     :          0,   -2,    2,    0,    1,
     :         -1,    0,    1,    2,    1,
     :         -1,    1,    0,    3,    0,
     :         -1,   -1,    2,    1,    2,
     :          0,   -1,    2,    0,    0,
     :         -2,    1,    2,    2,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=561,570 ) /
     :          2,   -2,    2,   -2,    2,
     :          1,    1,    0,    1,    1,
     :          1,    0,    1,    0,    1,
     :          1,    0,    1,    0,    0,
     :          0,    2,    0,    2,    0,
     :          2,   -1,    2,   -2,    1,
     :          0,   -1,    4,   -2,    1,
     :          0,    0,    4,   -2,    3,
     :          0,    1,    4,   -2,    1,
     :          4,    0,    2,   -4,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=571,580 ) /
     :          2,    2,    2,   -2,    2,
     :          2,    0,    4,   -4,    2,
     :         -1,   -2,    0,    4,    0,
     :         -1,   -3,    2,    2,    2,
     :         -3,    0,    2,    4,    2,
     :         -3,    0,    2,   -2,    1,
     :         -1,   -1,    0,   -2,    1,
     :         -3,    0,    0,    0,    2,
     :         -3,    0,   -2,    2,    0,
     :          0,    1,    0,   -4,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=581,590 ) /
     :         -2,    1,    0,   -2,    1,
     :         -4,    0,    0,    0,    1,
     :         -1,    0,    0,   -4,    1,
     :         -3,    0,    0,   -2,    1,
     :          0,    0,    0,    3,    2,
     :         -1,    1,    0,    4,    1,
     :          1,   -2,    2,    0,    1,
     :          0,    1,    0,    3,    0,
     :         -1,    0,    2,    2,    3,
     :          0,    0,    2,    2,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=591,600 ) /
     :         -2,    0,    2,    2,    2,
     :         -1,    1,    2,    2,    0,
     :          3,    0,    0,    0,    2,
     :          2,    1,    0,    1,    0,
     :          2,   -1,    2,   -1,    2,
     :          0,    0,    2,    0,    1,
     :          0,    0,    3,    0,    3,
     :          0,    0,    3,    0,    2,
     :         -1,    2,    2,    2,    1,
     :         -1,    0,    4,    0,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=601,610 ) /
     :          1,    2,    2,    0,    1,
     :          3,    1,    2,   -2,    1,
     :          1,    1,    4,   -2,    2,
     :         -2,   -1,    0,    6,    0,
     :          0,   -2,    0,    4,    0,
     :         -2,    0,    0,    6,    1,
     :         -2,   -2,    2,    4,    2,
     :          0,   -3,    2,    2,    2,
     :          0,    0,    0,    4,    2,
     :         -1,   -1,    2,    3,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=611,620 ) /
     :         -2,    0,    2,    4,    0,
     :          2,   -1,    0,    2,    1,
     :          1,    0,    0,    3,    0,
     :          0,    1,    0,    4,    1,
     :          0,    1,    0,    4,    0,
     :          1,   -1,    2,    1,    2,
     :          0,    0,    2,    2,    3,
     :          1,    0,    2,    2,    2,
     :         -1,    0,    2,    2,    2,
     :         -2,    0,    4,    2,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=621,630 ) /
     :          2,    1,    0,    2,    1,
     :          2,    1,    0,    2,    0,
     :          2,   -1,    2,    0,    0,
     :          1,    0,    2,    1,    0,
     :          0,    1,    2,    2,    0,
     :          2,    0,    2,    0,    3,
     :          3,    0,    2,    0,    2,
     :          1,    0,    2,    0,    2,
     :          1,    0,    3,    0,    3,
     :          1,    1,    2,    1,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=631,640 ) /
     :          0,    2,    2,    2,    2,
     :          2,    1,    2,    0,    0,
     :          2,    0,    4,   -2,    1,
     :          4,    1,    2,   -2,    2,
     :         -1,   -1,    0,    6,    0,
     :         -3,   -1,    2,    6,    2,
     :         -1,    0,    0,    6,    1,
     :         -3,    0,    2,    6,    1,
     :          1,   -1,    0,    4,    1,
     :          1,   -1,    0,    4,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=641,650 ) /
     :         -2,    0,    2,    5,    2,
     :          1,   -2,    2,    2,    1,
     :          3,   -1,    0,    2,    0,
     :          1,   -1,    2,    2,    0,
     :          0,    0,    2,    3,    1,
     :         -1,    1,    2,    4,    1,
     :          0,    1,    2,    3,    2,
     :         -1,    0,    4,    2,    1,
     :          2,    0,    2,    1,    1,
     :          5,    0,    0,    0,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=651,660 ) /
     :          2,    1,    2,    1,    2,
     :          1,    0,    4,    0,    1,
     :          3,    1,    2,    0,    1,
     :          3,    0,    4,   -2,    2,
     :         -2,   -1,    2,    6,    2,
     :          0,    0,    0,    6,    0,
     :          0,   -2,    2,    4,    2,
     :         -2,    0,    2,    6,    1,
     :          2,    0,    0,    4,    1,
     :          2,    0,    0,    4,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=661,670 ) /
     :          2,   -2,    2,    2,    2,
     :          0,    0,    2,    4,    0,
     :          1,    0,    2,    3,    2,
     :          4,    0,    0,    2,    0,
     :          2,    0,    2,    2,    0,
     :          0,    0,    4,    2,    2,
     :          4,   -1,    2,    0,    2,
     :          3,    0,    2,    1,    2,
     :          2,    1,    2,    2,    1,
     :          4,    1,    2,    0,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=671,678 ) /
     :         -1,   -1,    2,    6,    2,
     :         -1,    0,    2,    6,    1,
     :          1,   -1,    2,    4,    1,
     :          1,    1,    2,    4,    2,
     :          3,    1,    2,    2,    2,
     :          5,    0,    2,    0,    1,
     :          2,   -1,    2,    4,    2,
     :          2,    0,    2,    4,    1 /

*
*  Luni-Solar nutation coefficients, unit 1e-7 arcsec
*  longitude (sin, t*sin, cos), obliquity (cos, t*cos, sin)
*

      DATA ( ( CLS(I,J), I=1,6 ), J=  1, 10 ) /
     : -172064161D0, -174666D0,  33386D0, 92052331D0,  9086D0, 15377D0,
     :  -13170906D0,   -1675D0, -13696D0,  5730336D0, -3015D0, -4587D0,
     :   -2276413D0,    -234D0,   2796D0,   978459D0,  -485D0,  1374D0,
     :    2074554D0,     207D0,   -698D0,  -897492D0,   470D0,  -291D0,
     :    1475877D0,   -3633D0,  11817D0,    73871D0,  -184D0, -1924D0,
     :    -516821D0,    1226D0,   -524D0,   224386D0,  -677D0,  -174D0,
     :     711159D0,      73D0,   -872D0,    -6750D0,     0D0,   358D0,
     :    -387298D0,    -367D0,    380D0,   200728D0,    18D0,   318D0,
     :    -301461D0,     -36D0,    816D0,   129025D0,   -63D0,   367D0,
     :     215829D0,    -494D0,    111D0,   -95929D0,   299D0,   132D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J= 11, 20 ) /
     :     128227D0,     137D0,    181D0,   -68982D0,    -9D0,    39D0,
     :     123457D0,      11D0,     19D0,   -53311D0,    32D0,    -4D0,
     :     156994D0,      10D0,   -168D0,    -1235D0,     0D0,    82D0,
     :      63110D0,      63D0,     27D0,   -33228D0,     0D0,    -9D0,
     :     -57976D0,     -63D0,   -189D0,    31429D0,     0D0,   -75D0,
     :     -59641D0,     -11D0,    149D0,    25543D0,   -11D0,    66D0,
     :     -51613D0,     -42D0,    129D0,    26366D0,     0D0,    78D0,
     :      45893D0,      50D0,     31D0,   -24236D0,   -10D0,    20D0,
     :      63384D0,      11D0,   -150D0,    -1220D0,     0D0,    29D0,
     :     -38571D0,      -1D0,    158D0,    16452D0,   -11D0,    68D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J= 21, 30 ) /
     :      32481D0,       0D0,      0D0,   -13870D0,     0D0,     0D0,
     :     -47722D0,       0D0,    -18D0,      477D0,     0D0,   -25D0,
     :     -31046D0,      -1D0,    131D0,    13238D0,   -11D0,    59D0,
     :      28593D0,       0D0,     -1D0,   -12338D0,    10D0,    -3D0,
     :      20441D0,      21D0,     10D0,   -10758D0,     0D0,    -3D0,
     :      29243D0,       0D0,    -74D0,     -609D0,     0D0,    13D0,
     :      25887D0,       0D0,    -66D0,     -550D0,     0D0,    11D0,
     :     -14053D0,     -25D0,     79D0,     8551D0,    -2D0,   -45D0,
     :      15164D0,      10D0,     11D0,    -8001D0,     0D0,    -1D0,
     :     -15794D0,      72D0,    -16D0,     6850D0,   -42D0,    -5D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J= 31, 40 ) /
     :      21783D0,       0D0,     13D0,     -167D0,     0D0,    13D0,
     :     -12873D0,     -10D0,    -37D0,     6953D0,     0D0,   -14D0,
     :     -12654D0,      11D0,     63D0,     6415D0,     0D0,    26D0,
     :     -10204D0,       0D0,     25D0,     5222D0,     0D0,    15D0,
     :      16707D0,     -85D0,    -10D0,      168D0,    -1D0,    10D0,
     :      -7691D0,       0D0,     44D0,     3268D0,     0D0,    19D0,
     :     -11024D0,       0D0,    -14D0,      104D0,     0D0,     2D0,
     :       7566D0,     -21D0,    -11D0,    -3250D0,     0D0,    -5D0,
     :      -6637D0,     -11D0,     25D0,     3353D0,     0D0,    14D0,
     :      -7141D0,      21D0,      8D0,     3070D0,     0D0,     4D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J= 41, 50 ) /
     :      -6302D0,     -11D0,      2D0,     3272D0,     0D0,     4D0,
     :       5800D0,      10D0,      2D0,    -3045D0,     0D0,    -1D0,
     :       6443D0,       0D0,     -7D0,    -2768D0,     0D0,    -4D0,
     :      -5774D0,     -11D0,    -15D0,     3041D0,     0D0,    -5D0,
     :      -5350D0,       0D0,     21D0,     2695D0,     0D0,    12D0,
     :      -4752D0,     -11D0,     -3D0,     2719D0,     0D0,    -3D0,
     :      -4940D0,     -11D0,    -21D0,     2720D0,     0D0,    -9D0,
     :       7350D0,       0D0,     -8D0,      -51D0,     0D0,     4D0,
     :       4065D0,       0D0,      6D0,    -2206D0,     0D0,     1D0,
     :       6579D0,       0D0,    -24D0,     -199D0,     0D0,     2D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J= 51, 60 ) /
     :       3579D0,       0D0,      5D0,    -1900D0,     0D0,     1D0,
     :       4725D0,       0D0,     -6D0,      -41D0,     0D0,     3D0,
     :      -3075D0,       0D0,     -2D0,     1313D0,     0D0,    -1D0,
     :      -2904D0,       0D0,     15D0,     1233D0,     0D0,     7D0,
     :       4348D0,       0D0,    -10D0,      -81D0,     0D0,     2D0,
     :      -2878D0,       0D0,      8D0,     1232D0,     0D0,     4D0,
     :      -4230D0,       0D0,      5D0,      -20D0,     0D0,    -2D0,
     :      -2819D0,       0D0,      7D0,     1207D0,     0D0,     3D0,
     :      -4056D0,       0D0,      5D0,       40D0,     0D0,    -2D0,
     :      -2647D0,       0D0,     11D0,     1129D0,     0D0,     5D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J= 61, 70 ) /
     :      -2294D0,       0D0,    -10D0,     1266D0,     0D0,    -4D0,
     :       2481D0,       0D0,     -7D0,    -1062D0,     0D0,    -3D0,
     :       2179D0,       0D0,     -2D0,    -1129D0,     0D0,    -2D0,
     :       3276D0,       0D0,      1D0,       -9D0,     0D0,     0D0,
     :      -3389D0,       0D0,      5D0,       35D0,     0D0,    -2D0,
     :       3339D0,       0D0,    -13D0,     -107D0,     0D0,     1D0,
     :      -1987D0,       0D0,     -6D0,     1073D0,     0D0,    -2D0,
     :      -1981D0,       0D0,      0D0,      854D0,     0D0,     0D0,
     :       4026D0,       0D0,   -353D0,     -553D0,     0D0,  -139D0,
     :       1660D0,       0D0,     -5D0,     -710D0,     0D0,    -2D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J= 71, 80 ) /
     :      -1521D0,       0D0,      9D0,      647D0,     0D0,     4D0,
     :       1314D0,       0D0,      0D0,     -700D0,     0D0,     0D0,
     :      -1283D0,       0D0,      0D0,      672D0,     0D0,     0D0,
     :      -1331D0,       0D0,      8D0,      663D0,     0D0,     4D0,
     :       1383D0,       0D0,     -2D0,     -594D0,     0D0,    -2D0,
     :       1405D0,       0D0,      4D0,     -610D0,     0D0,     2D0,
     :       1290D0,       0D0,      0D0,     -556D0,     0D0,     0D0,
     :      -1214D0,       0D0,      5D0,      518D0,     0D0,     2D0,
     :       1146D0,       0D0,     -3D0,     -490D0,     0D0,    -1D0,
     :       1019D0,       0D0,     -1D0,     -527D0,     0D0,    -1D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J= 81, 90 ) /
     :      -1100D0,       0D0,      9D0,      465D0,     0D0,     4D0,
     :       -970D0,       0D0,      2D0,      496D0,     0D0,     1D0,
     :       1575D0,       0D0,     -6D0,      -50D0,     0D0,     0D0,
     :        934D0,       0D0,     -3D0,     -399D0,     0D0,    -1D0,
     :        922D0,       0D0,     -1D0,     -395D0,     0D0,    -1D0,
     :        815D0,       0D0,     -1D0,     -422D0,     0D0,    -1D0,
     :        834D0,       0D0,      2D0,     -440D0,     0D0,     1D0,
     :       1248D0,       0D0,      0D0,     -170D0,     0D0,     1D0,
     :       1338D0,       0D0,     -5D0,      -39D0,     0D0,     0D0,
     :        716D0,       0D0,     -2D0,     -389D0,     0D0,    -1D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J= 91,100 ) /
     :       1282D0,       0D0,     -3D0,      -23D0,     0D0,     1D0,
     :        742D0,       0D0,      1D0,     -391D0,     0D0,     0D0,
     :       1020D0,       0D0,    -25D0,     -495D0,     0D0,   -10D0,
     :        715D0,       0D0,     -4D0,     -326D0,     0D0,     2D0,
     :       -666D0,       0D0,     -3D0,      369D0,     0D0,    -1D0,
     :       -667D0,       0D0,      1D0,      346D0,     0D0,     1D0,
     :       -704D0,       0D0,      0D0,      304D0,     0D0,     0D0,
     :       -694D0,       0D0,      5D0,      294D0,     0D0,     2D0,
     :      -1014D0,       0D0,     -1D0,        4D0,     0D0,    -1D0,
     :       -585D0,       0D0,     -2D0,      316D0,     0D0,    -1D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=101,110 ) /
     :       -949D0,       0D0,      1D0,        8D0,     0D0,    -1D0,
     :       -595D0,       0D0,      0D0,      258D0,     0D0,     0D0,
     :        528D0,       0D0,      0D0,     -279D0,     0D0,     0D0,
     :       -590D0,       0D0,      4D0,      252D0,     0D0,     2D0,
     :        570D0,       0D0,     -2D0,     -244D0,     0D0,    -1D0,
     :       -502D0,       0D0,      3D0,      250D0,     0D0,     2D0,
     :       -875D0,       0D0,      1D0,       29D0,     0D0,     0D0,
     :       -492D0,       0D0,     -3D0,      275D0,     0D0,    -1D0,
     :        535D0,       0D0,     -2D0,     -228D0,     0D0,    -1D0,
     :       -467D0,       0D0,      1D0,      240D0,     0D0,     1D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=111,120 ) /
     :        591D0,       0D0,      0D0,     -253D0,     0D0,     0D0,
     :       -453D0,       0D0,     -1D0,      244D0,     0D0,    -1D0,
     :        766D0,       0D0,      1D0,        9D0,     0D0,     0D0,
     :       -446D0,       0D0,      2D0,      225D0,     0D0,     1D0,
     :       -488D0,       0D0,      2D0,      207D0,     0D0,     1D0,
     :       -468D0,       0D0,      0D0,      201D0,     0D0,     0D0,
     :       -421D0,       0D0,      1D0,      216D0,     0D0,     1D0,
     :        463D0,       0D0,      0D0,     -200D0,     0D0,     0D0,
     :       -673D0,       0D0,      2D0,       14D0,     0D0,     0D0,
     :        658D0,       0D0,      0D0,       -2D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=121,130 ) /
     :       -438D0,       0D0,      0D0,      188D0,     0D0,     0D0,
     :       -390D0,       0D0,      0D0,      205D0,     0D0,     0D0,
     :        639D0,     -11D0,     -2D0,      -19D0,     0D0,     0D0,
     :        412D0,       0D0,     -2D0,     -176D0,     0D0,    -1D0,
     :       -361D0,       0D0,      0D0,      189D0,     0D0,     0D0,
     :        360D0,       0D0,     -1D0,     -185D0,     0D0,    -1D0,
     :        588D0,       0D0,     -3D0,      -24D0,     0D0,     0D0,
     :       -578D0,       0D0,      1D0,        5D0,     0D0,     0D0,
     :       -396D0,       0D0,      0D0,      171D0,     0D0,     0D0,
     :        565D0,       0D0,     -1D0,       -6D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=131,140 ) /
     :       -335D0,       0D0,     -1D0,      184D0,     0D0,    -1D0,
     :        357D0,       0D0,      1D0,     -154D0,     0D0,     0D0,
     :        321D0,       0D0,      1D0,     -174D0,     0D0,     0D0,
     :       -301D0,       0D0,     -1D0,      162D0,     0D0,     0D0,
     :       -334D0,       0D0,      0D0,      144D0,     0D0,     0D0,
     :        493D0,       0D0,     -2D0,      -15D0,     0D0,     0D0,
     :        494D0,       0D0,     -2D0,      -19D0,     0D0,     0D0,
     :        337D0,       0D0,     -1D0,     -143D0,     0D0,    -1D0,
     :        280D0,       0D0,     -1D0,     -144D0,     0D0,     0D0,
     :        309D0,       0D0,      1D0,     -134D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=141,150 ) /
     :       -263D0,       0D0,      2D0,      131D0,     0D0,     1D0,
     :        253D0,       0D0,      1D0,     -138D0,     0D0,     0D0,
     :        245D0,       0D0,      0D0,     -128D0,     0D0,     0D0,
     :        416D0,       0D0,     -2D0,      -17D0,     0D0,     0D0,
     :       -229D0,       0D0,      0D0,      128D0,     0D0,     0D0,
     :        231D0,       0D0,      0D0,     -120D0,     0D0,     0D0,
     :       -259D0,       0D0,      2D0,      109D0,     0D0,     1D0,
     :        375D0,       0D0,     -1D0,       -8D0,     0D0,     0D0,
     :        252D0,       0D0,      0D0,     -108D0,     0D0,     0D0,
     :       -245D0,       0D0,      1D0,      104D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=151,160 ) /
     :        243D0,       0D0,     -1D0,     -104D0,     0D0,     0D0,
     :        208D0,       0D0,      1D0,     -112D0,     0D0,     0D0,
     :        199D0,       0D0,      0D0,     -102D0,     0D0,     0D0,
     :       -208D0,       0D0,      1D0,      105D0,     0D0,     0D0,
     :        335D0,       0D0,     -2D0,      -14D0,     0D0,     0D0,
     :       -325D0,       0D0,      1D0,        7D0,     0D0,     0D0,
     :       -187D0,       0D0,      0D0,       96D0,     0D0,     0D0,
     :        197D0,       0D0,     -1D0,     -100D0,     0D0,     0D0,
     :       -192D0,       0D0,      2D0,       94D0,     0D0,     1D0,
     :       -188D0,       0D0,      0D0,       83D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=161,170 ) /
     :        276D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :       -286D0,       0D0,      1D0,        6D0,     0D0,     0D0,
     :        186D0,       0D0,     -1D0,      -79D0,     0D0,     0D0,
     :       -219D0,       0D0,      0D0,       43D0,     0D0,     0D0,
     :        276D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :       -153D0,       0D0,     -1D0,       84D0,     0D0,     0D0,
     :       -156D0,       0D0,      0D0,       81D0,     0D0,     0D0,
     :       -154D0,       0D0,      1D0,       78D0,     0D0,     0D0,
     :       -174D0,       0D0,      1D0,       75D0,     0D0,     0D0,
     :       -163D0,       0D0,      2D0,       69D0,     0D0,     1D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=171,180 ) /
     :       -228D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :         91D0,       0D0,     -4D0,      -54D0,     0D0,    -2D0,
     :        175D0,       0D0,      0D0,      -75D0,     0D0,     0D0,
     :       -159D0,       0D0,      0D0,       69D0,     0D0,     0D0,
     :        141D0,       0D0,      0D0,      -72D0,     0D0,     0D0,
     :        147D0,       0D0,      0D0,      -75D0,     0D0,     0D0,
     :       -132D0,       0D0,      0D0,       69D0,     0D0,     0D0,
     :        159D0,       0D0,    -28D0,      -54D0,     0D0,    11D0,
     :        213D0,       0D0,      0D0,       -4D0,     0D0,     0D0,
     :        123D0,       0D0,      0D0,      -64D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=181,190 ) /
     :       -118D0,       0D0,     -1D0,       66D0,     0D0,     0D0,
     :        144D0,       0D0,     -1D0,      -61D0,     0D0,     0D0,
     :       -121D0,       0D0,      1D0,       60D0,     0D0,     0D0,
     :       -134D0,       0D0,      1D0,       56D0,     0D0,     1D0,
     :       -105D0,       0D0,      0D0,       57D0,     0D0,     0D0,
     :       -102D0,       0D0,      0D0,       56D0,     0D0,     0D0,
     :        120D0,       0D0,      0D0,      -52D0,     0D0,     0D0,
     :        101D0,       0D0,      0D0,      -54D0,     0D0,     0D0,
     :       -113D0,       0D0,      0D0,       59D0,     0D0,     0D0,
     :       -106D0,       0D0,      0D0,       61D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=191,200 ) /
     :       -129D0,       0D0,      1D0,       55D0,     0D0,     0D0,
     :       -114D0,       0D0,      0D0,       57D0,     0D0,     0D0,
     :        113D0,       0D0,     -1D0,      -49D0,     0D0,     0D0,
     :       -102D0,       0D0,      0D0,       44D0,     0D0,     0D0,
     :        -94D0,       0D0,      0D0,       51D0,     0D0,     0D0,
     :       -100D0,       0D0,     -1D0,       56D0,     0D0,     0D0,
     :         87D0,       0D0,      0D0,      -47D0,     0D0,     0D0,
     :        161D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :         96D0,       0D0,      0D0,      -50D0,     0D0,     0D0,
     :        151D0,       0D0,     -1D0,       -5D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=201,210 ) /
     :       -104D0,       0D0,      0D0,       44D0,     0D0,     0D0,
     :       -110D0,       0D0,      0D0,       48D0,     0D0,     0D0,
     :       -100D0,       0D0,      1D0,       50D0,     0D0,     0D0,
     :         92D0,       0D0,     -5D0,       12D0,     0D0,    -2D0,
     :         82D0,       0D0,      0D0,      -45D0,     0D0,     0D0,
     :         82D0,       0D0,      0D0,      -45D0,     0D0,     0D0,
     :        -78D0,       0D0,      0D0,       41D0,     0D0,     0D0,
     :        -77D0,       0D0,      0D0,       43D0,     0D0,     0D0,
     :          2D0,       0D0,      0D0,       54D0,     0D0,     0D0,
     :         94D0,       0D0,      0D0,      -40D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=211,220 ) /
     :        -93D0,       0D0,      0D0,       40D0,     0D0,     0D0,
     :        -83D0,       0D0,     10D0,       40D0,     0D0,    -2D0,
     :         83D0,       0D0,      0D0,      -36D0,     0D0,     0D0,
     :        -91D0,       0D0,      0D0,       39D0,     0D0,     0D0,
     :        128D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :        -79D0,       0D0,      0D0,       34D0,     0D0,     0D0,
     :        -83D0,       0D0,      0D0,       47D0,     0D0,     0D0,
     :         84D0,       0D0,      0D0,      -44D0,     0D0,     0D0,
     :         83D0,       0D0,      0D0,      -43D0,     0D0,     0D0,
     :         91D0,       0D0,      0D0,      -39D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=221,230 ) /
     :        -77D0,       0D0,      0D0,       39D0,     0D0,     0D0,
     :         84D0,       0D0,      0D0,      -43D0,     0D0,     0D0,
     :        -92D0,       0D0,      1D0,       39D0,     0D0,     0D0,
     :        -92D0,       0D0,      1D0,       39D0,     0D0,     0D0,
     :        -94D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         68D0,       0D0,      0D0,      -36D0,     0D0,     0D0,
     :        -61D0,       0D0,      0D0,       32D0,     0D0,     0D0,
     :         71D0,       0D0,      0D0,      -31D0,     0D0,     0D0,
     :         62D0,       0D0,      0D0,      -34D0,     0D0,     0D0,
     :        -63D0,       0D0,      0D0,       33D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=231,240 ) /
     :        -73D0,       0D0,      0D0,       32D0,     0D0,     0D0,
     :        115D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :       -103D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         63D0,       0D0,      0D0,      -28D0,     0D0,     0D0,
     :         74D0,       0D0,      0D0,      -32D0,     0D0,     0D0,
     :       -103D0,       0D0,     -3D0,        3D0,     0D0,    -1D0,
     :        -69D0,       0D0,      0D0,       30D0,     0D0,     0D0,
     :         57D0,       0D0,      0D0,      -29D0,     0D0,     0D0,
     :         94D0,       0D0,      0D0,       -4D0,     0D0,     0D0,
     :         64D0,       0D0,      0D0,      -33D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=241,250 ) /
     :        -63D0,       0D0,      0D0,       26D0,     0D0,     0D0,
     :        -38D0,       0D0,      0D0,       20D0,     0D0,     0D0,
     :        -43D0,       0D0,      0D0,       24D0,     0D0,     0D0,
     :        -45D0,       0D0,      0D0,       23D0,     0D0,     0D0,
     :         47D0,       0D0,      0D0,      -24D0,     0D0,     0D0,
     :        -48D0,       0D0,      0D0,       25D0,     0D0,     0D0,
     :         45D0,       0D0,      0D0,      -26D0,     0D0,     0D0,
     :         56D0,       0D0,      0D0,      -25D0,     0D0,     0D0,
     :         88D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :        -75D0,       0D0,      0D0,        0D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=251,260 ) /
     :         85D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         49D0,       0D0,      0D0,      -26D0,     0D0,     0D0,
     :        -74D0,       0D0,     -3D0,       -1D0,     0D0,    -1D0,
     :        -39D0,       0D0,      0D0,       21D0,     0D0,     0D0,
     :         45D0,       0D0,      0D0,      -20D0,     0D0,     0D0,
     :         51D0,       0D0,      0D0,      -22D0,     0D0,     0D0,
     :        -40D0,       0D0,      0D0,       21D0,     0D0,     0D0,
     :         41D0,       0D0,      0D0,      -21D0,     0D0,     0D0,
     :        -42D0,       0D0,      0D0,       24D0,     0D0,     0D0,
     :        -51D0,       0D0,      0D0,       22D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=261,270 ) /
     :        -42D0,       0D0,      0D0,       22D0,     0D0,     0D0,
     :         39D0,       0D0,      0D0,      -21D0,     0D0,     0D0,
     :         46D0,       0D0,      0D0,      -18D0,     0D0,     0D0,
     :        -53D0,       0D0,      0D0,       22D0,     0D0,     0D0,
     :         82D0,       0D0,      0D0,       -4D0,     0D0,     0D0,
     :         81D0,       0D0,     -1D0,       -4D0,     0D0,     0D0,
     :         47D0,       0D0,      0D0,      -19D0,     0D0,     0D0,
     :         53D0,       0D0,      0D0,      -23D0,     0D0,     0D0,
     :        -45D0,       0D0,      0D0,       22D0,     0D0,     0D0,
     :        -44D0,       0D0,      0D0,       -2D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=271,280 ) /
     :        -33D0,       0D0,      0D0,       16D0,     0D0,     0D0,
     :        -61D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :         28D0,       0D0,      0D0,      -15D0,     0D0,     0D0,
     :        -38D0,       0D0,      0D0,       19D0,     0D0,     0D0,
     :        -33D0,       0D0,      0D0,       21D0,     0D0,     0D0,
     :        -60D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         48D0,       0D0,      0D0,      -10D0,     0D0,     0D0,
     :         27D0,       0D0,      0D0,      -14D0,     0D0,     0D0,
     :         38D0,       0D0,      0D0,      -20D0,     0D0,     0D0,
     :         31D0,       0D0,      0D0,      -13D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=281,290 ) /
     :        -29D0,       0D0,      0D0,       15D0,     0D0,     0D0,
     :         28D0,       0D0,      0D0,      -15D0,     0D0,     0D0,
     :        -32D0,       0D0,      0D0,       15D0,     0D0,     0D0,
     :         45D0,       0D0,      0D0,       -8D0,     0D0,     0D0,
     :        -44D0,       0D0,      0D0,       19D0,     0D0,     0D0,
     :         28D0,       0D0,      0D0,      -15D0,     0D0,     0D0,
     :        -51D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -36D0,       0D0,      0D0,       20D0,     0D0,     0D0,
     :         44D0,       0D0,      0D0,      -19D0,     0D0,     0D0,
     :         26D0,       0D0,      0D0,      -14D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=291,300 ) /
     :        -60D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         35D0,       0D0,      0D0,      -18D0,     0D0,     0D0,
     :        -27D0,       0D0,      0D0,       11D0,     0D0,     0D0,
     :         47D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :         36D0,       0D0,      0D0,      -15D0,     0D0,     0D0,
     :        -36D0,       0D0,      0D0,       20D0,     0D0,     0D0,
     :        -35D0,       0D0,      0D0,       19D0,     0D0,     0D0,
     :        -37D0,       0D0,      0D0,       19D0,     0D0,     0D0,
     :         32D0,       0D0,      0D0,      -16D0,     0D0,     0D0,
     :         35D0,       0D0,      0D0,      -14D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=301,310 ) /
     :         32D0,       0D0,      0D0,      -13D0,     0D0,     0D0,
     :         65D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :         47D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :         32D0,       0D0,      0D0,      -16D0,     0D0,     0D0,
     :         37D0,       0D0,      0D0,      -16D0,     0D0,     0D0,
     :        -30D0,       0D0,      0D0,       15D0,     0D0,     0D0,
     :        -32D0,       0D0,      0D0,       16D0,     0D0,     0D0,
     :        -31D0,       0D0,      0D0,       13D0,     0D0,     0D0,
     :         37D0,       0D0,      0D0,      -16D0,     0D0,     0D0,
     :         31D0,       0D0,      0D0,      -13D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=311,320 ) /
     :         49D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :         32D0,       0D0,      0D0,      -13D0,     0D0,     0D0,
     :         23D0,       0D0,      0D0,      -12D0,     0D0,     0D0,
     :        -43D0,       0D0,      0D0,       18D0,     0D0,     0D0,
     :         26D0,       0D0,      0D0,      -11D0,     0D0,     0D0,
     :        -32D0,       0D0,      0D0,       14D0,     0D0,     0D0,
     :        -29D0,       0D0,      0D0,       14D0,     0D0,     0D0,
     :        -27D0,       0D0,      0D0,       12D0,     0D0,     0D0,
     :         30D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -11D0,       0D0,      0D0,        5D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=321,330 ) /
     :        -21D0,       0D0,      0D0,       10D0,     0D0,     0D0,
     :        -34D0,       0D0,      0D0,       15D0,     0D0,     0D0,
     :        -10D0,       0D0,      0D0,        6D0,     0D0,     0D0,
     :        -36D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -9D0,       0D0,      0D0,        4D0,     0D0,     0D0,
     :        -12D0,       0D0,      0D0,        5D0,     0D0,     0D0,
     :        -21D0,       0D0,      0D0,        5D0,     0D0,     0D0,
     :        -29D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :        -15D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :        -20D0,       0D0,      0D0,        0D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=331,340 ) /
     :         28D0,       0D0,      0D0,        0D0,     0D0,    -2D0,
     :         17D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -22D0,       0D0,      0D0,       12D0,     0D0,     0D0,
     :        -14D0,       0D0,      0D0,        7D0,     0D0,     0D0,
     :         24D0,       0D0,      0D0,      -11D0,     0D0,     0D0,
     :         11D0,       0D0,      0D0,       -6D0,     0D0,     0D0,
     :         14D0,       0D0,      0D0,       -6D0,     0D0,     0D0,
     :         24D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         18D0,       0D0,      0D0,       -8D0,     0D0,     0D0,
     :        -38D0,       0D0,      0D0,        0D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=341,350 ) /
     :        -31D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -16D0,       0D0,      0D0,        8D0,     0D0,     0D0,
     :         29D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -18D0,       0D0,      0D0,       10D0,     0D0,     0D0,
     :        -10D0,       0D0,      0D0,        5D0,     0D0,     0D0,
     :        -17D0,       0D0,      0D0,       10D0,     0D0,     0D0,
     :          9D0,       0D0,      0D0,       -4D0,     0D0,     0D0,
     :         16D0,       0D0,      0D0,       -6D0,     0D0,     0D0,
     :         22D0,       0D0,      0D0,      -12D0,     0D0,     0D0,
     :         20D0,       0D0,      0D0,        0D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=351,360 ) /
     :        -13D0,       0D0,      0D0,        6D0,     0D0,     0D0,
     :        -17D0,       0D0,      0D0,        9D0,     0D0,     0D0,
     :        -14D0,       0D0,      0D0,        8D0,     0D0,     0D0,
     :          0D0,       0D0,      0D0,       -7D0,     0D0,     0D0,
     :         14D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         19D0,       0D0,      0D0,      -10D0,     0D0,     0D0,
     :        -34D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -20D0,       0D0,      0D0,        8D0,     0D0,     0D0,
     :          9D0,       0D0,      0D0,       -5D0,     0D0,     0D0,
     :        -18D0,       0D0,      0D0,        7D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=361,370 ) /
     :         13D0,       0D0,      0D0,       -6D0,     0D0,     0D0,
     :         17D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -12D0,       0D0,      0D0,        5D0,     0D0,     0D0,
     :         15D0,       0D0,      0D0,       -8D0,     0D0,     0D0,
     :        -11D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :         13D0,       0D0,      0D0,       -5D0,     0D0,     0D0,
     :        -18D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -35D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          9D0,       0D0,      0D0,       -4D0,     0D0,     0D0,
     :        -19D0,       0D0,      0D0,       10D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=371,380 ) /
     :        -26D0,       0D0,      0D0,       11D0,     0D0,     0D0,
     :          8D0,       0D0,      0D0,       -4D0,     0D0,     0D0,
     :        -10D0,       0D0,      0D0,        4D0,     0D0,     0D0,
     :         10D0,       0D0,      0D0,       -6D0,     0D0,     0D0,
     :        -21D0,       0D0,      0D0,        9D0,     0D0,     0D0,
     :        -15D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          9D0,       0D0,      0D0,       -5D0,     0D0,     0D0,
     :        -29D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -19D0,       0D0,      0D0,       10D0,     0D0,     0D0,
     :         12D0,       0D0,      0D0,       -5D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=381,390 ) /
     :         22D0,       0D0,      0D0,       -9D0,     0D0,     0D0,
     :        -10D0,       0D0,      0D0,        5D0,     0D0,     0D0,
     :        -20D0,       0D0,      0D0,       11D0,     0D0,     0D0,
     :        -20D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -17D0,       0D0,      0D0,        7D0,     0D0,     0D0,
     :         15D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :          8D0,       0D0,      0D0,       -4D0,     0D0,     0D0,
     :         14D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -12D0,       0D0,      0D0,        6D0,     0D0,     0D0,
     :         25D0,       0D0,      0D0,        0D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=391,400 ) /
     :        -13D0,       0D0,      0D0,        6D0,     0D0,     0D0,
     :        -14D0,       0D0,      0D0,        8D0,     0D0,     0D0,
     :         13D0,       0D0,      0D0,       -5D0,     0D0,     0D0,
     :        -17D0,       0D0,      0D0,        9D0,     0D0,     0D0,
     :        -12D0,       0D0,      0D0,        6D0,     0D0,     0D0,
     :        -10D0,       0D0,      0D0,        5D0,     0D0,     0D0,
     :         10D0,       0D0,      0D0,       -6D0,     0D0,     0D0,
     :        -15D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -22D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         28D0,       0D0,      0D0,       -1D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=401,410 ) /
     :         15D0,       0D0,      0D0,       -7D0,     0D0,     0D0,
     :         23D0,       0D0,      0D0,      -10D0,     0D0,     0D0,
     :         12D0,       0D0,      0D0,       -5D0,     0D0,     0D0,
     :         29D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :        -25D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :         22D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -18D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         15D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :        -23D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         12D0,       0D0,      0D0,       -5D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=411,420 ) /
     :         -8D0,       0D0,      0D0,        4D0,     0D0,     0D0,
     :        -19D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -10D0,       0D0,      0D0,        4D0,     0D0,     0D0,
     :         21D0,       0D0,      0D0,       -9D0,     0D0,     0D0,
     :         23D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :        -16D0,       0D0,      0D0,        8D0,     0D0,     0D0,
     :        -19D0,       0D0,      0D0,        9D0,     0D0,     0D0,
     :        -22D0,       0D0,      0D0,       10D0,     0D0,     0D0,
     :         27D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :         16D0,       0D0,      0D0,       -8D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=421,430 ) /
     :         19D0,       0D0,      0D0,       -8D0,     0D0,     0D0,
     :          9D0,       0D0,      0D0,       -4D0,     0D0,     0D0,
     :         -9D0,       0D0,      0D0,        4D0,     0D0,     0D0,
     :         -9D0,       0D0,      0D0,        4D0,     0D0,     0D0,
     :         -8D0,       0D0,      0D0,        4D0,     0D0,     0D0,
     :         18D0,       0D0,      0D0,       -9D0,     0D0,     0D0,
     :         16D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :        -10D0,       0D0,      0D0,        4D0,     0D0,     0D0,
     :        -23D0,       0D0,      0D0,        9D0,     0D0,     0D0,
     :         16D0,       0D0,      0D0,       -1D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=431,440 ) /
     :        -12D0,       0D0,      0D0,        6D0,     0D0,     0D0,
     :         -8D0,       0D0,      0D0,        4D0,     0D0,     0D0,
     :         30D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :         24D0,       0D0,      0D0,      -10D0,     0D0,     0D0,
     :         10D0,       0D0,      0D0,       -4D0,     0D0,     0D0,
     :        -16D0,       0D0,      0D0,        7D0,     0D0,     0D0,
     :        -16D0,       0D0,      0D0,        7D0,     0D0,     0D0,
     :         17D0,       0D0,      0D0,       -7D0,     0D0,     0D0,
     :        -24D0,       0D0,      0D0,       10D0,     0D0,     0D0,
     :        -12D0,       0D0,      0D0,        5D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=441,450 ) /
     :        -24D0,       0D0,      0D0,       11D0,     0D0,     0D0,
     :        -23D0,       0D0,      0D0,        9D0,     0D0,     0D0,
     :        -13D0,       0D0,      0D0,        5D0,     0D0,     0D0,
     :        -15D0,       0D0,      0D0,        7D0,     0D0,     0D0,
     :          0D0,       0D0,  -1988D0,        0D0,     0D0, -1679D0,
     :          0D0,       0D0,    -63D0,        0D0,     0D0,   -27D0,
     :         -4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          0D0,       0D0,      5D0,        0D0,     0D0,     4D0,
     :          5D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :          0D0,       0D0,    364D0,        0D0,     0D0,   176D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=451,460 ) /
     :          0D0,       0D0,  -1044D0,        0D0,     0D0,  -891D0,
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          0D0,       0D0,    330D0,        0D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :         -5D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,        0D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=461,470 ) /
     :          3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          0D0,       0D0,      5D0,        0D0,     0D0,     0D0,
     :          0D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          6D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :         -7D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -12D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=471,480 ) /
     :         -5D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -7D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :          7D0,       0D0,      0D0,       -4D0,     0D0,     0D0,
     :          0D0,       0D0,    -12D0,        0D0,     0D0,   -10D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -7D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=481,490 ) /
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :          0D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :          7D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :         -5D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -5D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -2D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=491,500 ) /
     :         -8D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :          9D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          6D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :         -5D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -7D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        2D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=501,510 ) /
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :         -5D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          9D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          9D0,       0D0,      0D0,       -3D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=511,520 ) /
     :         -4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          8D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :          6D0,       0D0,      0D0,       -3D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=521,530 ) /
     :          3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :         -7D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          9D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -5D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :        -13D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -7D0,       0D0,      0D0,        0D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=531,540 ) /
     :         10D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :         10D0,       0D0,     13D0,        6D0,     0D0,    -5D0,
     :          0D0,       0D0,     30D0,        0D0,     0D0,    14D0,
     :          0D0,       0D0,   -162D0,        0D0,     0D0,  -138D0,
     :          0D0,       0D0,     75D0,        0D0,     0D0,     0D0,
     :         -7D0,       0D0,      0D0,        4D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -2D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=541,550 ) /
     :          5D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -5D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          6D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          9D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -7D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=551,560 ) /
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          7D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -6D0,       0D0,     -3D0,        3D0,     0D0,     1D0,
     :          0D0,       0D0,     -3D0,        0D0,     0D0,    -2D0,
     :         11D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :         11D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        2D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=561,570 ) /
     :         -1D0,       0D0,      3D0,        3D0,     0D0,    -1D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          0D0,       0D0,    -13D0,        0D0,     0D0,   -11D0,
     :          3D0,       0D0,      6D0,        0D0,     0D0,     0D0,
     :         -7D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :         -7D0,       0D0,      0D0,        3D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=571,580 ) /
     :          8D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         11D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          8D0,       0D0,      0D0,       -4D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :         11D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -6D0,       0D0,      0D0,        3D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=581,590 ) /
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -8D0,       0D0,      0D0,        4D0,     0D0,     0D0,
     :         -7D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :          6D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :         -6D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :          6D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          6D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -2D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=591,600 ) /
     :         -5D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          6D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          0D0,       0D0,    -26D0,        0D0,     0D0,   -11D0,
     :          0D0,       0D0,    -10D0,        0D0,     0D0,    -5D0,
     :          5D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :        -13D0,       0D0,      0D0,        0D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=601,610 ) /
     :          3D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          7D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -6D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -5D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -7D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -2D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=611,620 ) /
     :         13D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :        -11D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          6D0,       0D0,      0D0,       -3D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=621,630 ) /
     :          3D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :        -12D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :          0D0,       0D0,     -5D0,        0D0,     0D0,    -2D0,
     :         -7D0,       0D0,      0D0,        4D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=631,640 ) /
     :          6D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :         -5D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         12D0,       0D0,      0D0,        0D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=641,650 ) /
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          6D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :         -6D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          6D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :          6D0,       0D0,      0D0,        0D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=651,660 ) /
     :         -6D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          7D0,       0D0,      0D0,       -4D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :         -5D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -6D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :         -6D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         10D0,       0D0,      0D0,        0D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=661,670 ) /
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          7D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          7D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         11D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :         -6D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -2D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=671,678 ) /
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        2D0,     0D0,     0D0 /

*
*  Planetary argument multipliers
*    :         L   L'  F   D   Om  Me  Ve  E  Ma  Ju  Sa  Ur  Ne  pre

      DATA ( ( NAPL(I,J), I=1,14 ), J=  1, 10 ) /
     :         0,  0,  0,  0,  0,  0,  0,  8,-16,  4,  5,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -8, 16, -4, -5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  8,-16,  4,  5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  2,  2,
     :         0,  0,  0,  0,  0,  0,  0, -4,  8, -1, -5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  4, -8,  3,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  0,  3, -8,  3,  0,  0,  0,  0,
     :        -1,  0,  0,  0,  0,  0, 10, -3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0, -2,  6, -3,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  4, -8,  3,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J= 11, 20 ) /
     :         0,  0,  1, -1,  1,  0,  0, -5,  8, -3,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -4,  8, -3,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  4, -8,  1,  5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -5,  6,  4,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -5,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  2, -5,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -5,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0, -2,  5,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0, -2,  5,  0,  0,  1 /
      DATA ( ( NAPL(I,J), I=1,14 ), J= 21, 30 ) /
     :         0,  0,  0,  0,  0,  0,  0,  0,  0, -2,  5,  0,  0,  2,
     :         2,  0, -1, -1,  0,  0,  0,  3, -7,  0,  0,  0,  0,  0,
     :         1,  0,  0, -2,  0,  0, 19,-21,  3,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  2, -4,  0, -3,  0,  0,  0,  0,
     :         1,  0,  0, -1,  1,  0,  0, -1,  0,  2,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0, -4, 10,  0,  0,  0,
     :        -2,  0,  0,  2,  1,  0,  0,  2,  0,  0, -5,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  3, -7,  4,  0,  0,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  0,  1,  0,  1, -1,  0,  0,  0,
     :        -2,  0,  0,  2,  1,  0,  0,  2,  0, -2,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J= 31, 40 ) /
     :        -1,  0,  0,  0,  0,  0, 18,-16,  0,  0,  0,  0,  0,  0,
     :        -2,  0,  1,  1,  2,  0,  0,  1,  0, -2,  0,  0,  0,  0,
     :        -1,  0,  1, -1,  1,  0, 18,-17,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  0,  1,  1,  0,  0,  2, -2,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -8, 13,  0,  0,  0,  0,  0,  2,
     :         0,  0,  2, -2,  2,  0, -8, 11,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -8, 13,  0,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0, -8, 12,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  8,-13,  0,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  8,-14,  0,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J= 41, 50 ) /
     :         0,  0,  0,  0,  0,  0,  8,-13,  0,  0,  0,  0,  0,  1,
     :        -2,  0,  0,  2,  1,  0,  0,  2,  0, -4,  5,  0,  0,  0,
     :        -2,  0,  0,  2,  2,  0,  3, -3,  0,  0,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  2,  0, -3,  1,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  3, -5,  0,  2,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  2,  0, -4,  3,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0, -1,  2,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  2,  0,  0, -2,  2,  0,  0,  0,  0,  0,
     :        -1,  0,  1,  0,  1,  0,  3, -5,  0,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J= 51, 60 ) /
     :        -1,  0,  0,  1,  0,  0,  3, -4,  0,  0,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  2,  0, -2, -2,  0,  0,  0,
     :        -2,  0,  2,  0,  2,  0,  0, -5,  9,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  0, -1,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  0,  0,  2,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  1,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  2,
     :        -1,  0,  0,  1,  0,  0,  0,  3, -4,  0,  0,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  0,  1,  0,  0,  2,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J= 61, 70 ) /
     :         0,  0,  1, -1,  2,  0,  0, -1,  0,  0,  2,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0, -9, 17,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  2,  0, -3,  5,  0,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0, -1,  2,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0,
     :         1,  0,  0, -2,  0,  0, 17,-16,  0, -2,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  1, -3,  0,  0,  0,
     :        -2,  0,  0,  2,  1,  0,  0,  5, -6,  0,  0,  0,  0,  0,
     :         0,  0, -2,  2,  0,  0,  0,  9,-13,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  2,  0,  0, -1,  0,  0,  1,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J= 71, 80 ) /
     :         0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  1,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  0,  1,  0,  0,  1,  0,  0,  0,
     :         0,  0, -2,  2,  0,  0,  5, -6,  0,  0,  0,  0,  0,  0,
     :         0,  0, -1,  1,  1,  0,  5, -7,  0,  0,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  6, -8,  0,  0,  0,  0,  0,  0,
     :         2,  0,  1, -3,  1,  0, -6,  7,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  2,  0,  0,  0,  0,  1,  0,  0,  0,  0,
     :         0,  0, -1,  1,  1,  0,  0,  1,  0,  1,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  0,  2,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  1 /
      DATA ( ( NAPL(I,J), I=1,14 ), J= 81, 90 ) /
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -8, 15,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -8, 15,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  0, -9, 15,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  8,-15,  0,  0,  0,  0,  0,
     :         1,  0, -1, -1,  0,  0,  0,  8,-15,  0,  0,  0,  0,  0,
     :         2,  0,  0, -2,  0,  0,  2, -5,  0,  0,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  2,  0, -5,  5,  0,  0,  0,
     :         2,  0,  0, -2,  1,  0,  0, -6,  8,  0,  0,  0,  0,  0,
     :         2,  0,  0, -2,  1,  0,  0, -2,  0,  3,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J= 91,100 ) /
     :        -2,  0,  1,  1,  0,  0,  0,  1,  0, -3,  0,  0,  0,  0,
     :        -2,  0,  1,  1,  1,  0,  0,  1,  0, -3,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  2,  0, -3,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  6, -8,  0,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  2,  0, -1, -5,  0,  0,  0,
     :        -1,  0,  0,  1,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0,
     :        -1,  0,  1,  1,  1,  0,-20, 20,  0,  0,  0,  0,  0,  0,
     :         1,  0,  0, -2,  0,  0, 20,-21,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  8,-15,  0,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0,  0,-10, 15,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=101,110 ) /
     :         0,  0, -1,  1,  0,  0,  0,  1,  0,  1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  0,  0,  1,  0,  0,  0,  0,
     :         0,  0,  1, -1,  2,  0,  0, -1,  0,  1,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0, -2,  4,  0,  0,  0,
     :         2,  0,  0, -2,  1,  0, -6,  8,  0,  0,  0,  0,  0,  0,
     :         0,  0, -2,  2,  1,  0,  5, -6,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  0, -1,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  1,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=111,120 ) /
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  2,
     :         0,  0,  2, -2,  1,  0,  0, -9, 13,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  7,-13,  0,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  5, -6,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  9,-17,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -9, 17,  0,  0,  0,  0,  2,
     :         1,  0,  0, -1,  1,  0,  0, -3,  4,  0,  0,  0,  0,  0,
     :         1,  0,  0, -1,  1,  0, -3,  4,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  2,  0,  0, -1,  2,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=121,130 ) /
     :         0,  0, -1,  1,  1,  0,  0,  0,  2,  0,  0,  0,  0,  0,
     :         0,  0, -2,  2,  0,  1,  0, -2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  3, -5,  0,  2,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  1,  0,  0,  2,  0, -3,  1,  0,  0,  0,
     :        -2,  0,  0,  2,  1,  0,  3, -3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  8,-13,  0,  0,  0,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  8,-12,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0, -8, 11,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  0,  1,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0,
     :        -1,  0,  0,  0,  1,  0, 18,-16,  0,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=131,140 ) /
     :         0,  0,  1, -1,  1,  0,  0, -1,  0, -1,  1,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  3, -7,  4,  0,  0,  0,  0,  0,
     :        -2,  0,  1,  1,  1,  0,  0, -3,  7,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  2,  0,  0, -1,  0, -2,  5,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  0,  0, -2,  5,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0, -4,  8, -3,  0,  0,  0,  0,
     :         1,  0,  0,  0,  1,  0,-10,  3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0,  0, -2,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  0,  0,  1,  0, 10, -3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  4, -8,  3,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=141,150 ) /
     :         0,  0,  0,  0,  1,  0,  0,  0,  0,  2, -5,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  0,  1,  0,  2, -5,  0,  0,  0,
     :         2,  0, -1, -1,  1,  0,  0,  3, -7,  0,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  2,  0,  0, -5,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0, -3,  7, -4,  0,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,
     :         1,  0,  0,  0,  1,  0,-18, 16,  0,  0,  0,  0,  0,  0,
     :        -2,  0,  1,  1,  1,  0,  0,  1,  0, -2,  0,  0,  0,  0,
     :         0,  0,  1, -1,  2,  0, -8, 12,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0, -8, 13,  0,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=151,160 ) /
     :         0,  0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  0,  0, -2,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -2,  2,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -1,  2,  0,  0,  0,  0,  1,
     :        -1,  0,  0,  1,  1,  0,  3, -4,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  0,  1,  1,  0,  0,  3, -4,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  0, -2,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  2,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  1 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=161,170 ) /
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  2,
     :         0,  0,  1, -1,  0,  0,  3, -6,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0, -3,  5,  0,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  2,  0, -3,  4,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0, -2,  4,  0,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0, -5,  6,  0,  0,  0,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  5, -7,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  5, -8,  0,  0,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  1,  0,  6, -8,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0, -8, 15,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=171,180 ) /
     :        -2,  0,  0,  2,  1,  0,  0,  2,  0, -3,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  1,  0,  0,  6, -8,  0,  0,  0,  0,  0,
     :         1,  0,  0, -1,  1,  0,  0, -1,  0,  1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  3, -5,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0, -1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  1 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=181,190 ) /
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  2,
     :         0,  0,  1, -1,  2,  0,  0, -1,  0,  0, -1,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  0,  0,  0, -1,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  0,  1,  0,  0, -1,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -7, 13,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  7,-13,  0,  0,  0,  0,  0,
     :         2,  0,  0, -2,  1,  0,  0, -5,  6,  0,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0,  0, -8, 11,  0,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1, -1,  0,  2,  0,  0,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  4, -4,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=191,200 ) /
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  3,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  2,
     :        -2,  0,  0,  2,  0,  0,  3, -3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  2,  0,  0, -4,  8, -3,  0,  0,  0,  0,
     :         0,  0,  0,  0,  2,  0,  0,  4, -8,  3,  0,  0,  0,  0,
     :         2,  0,  0, -2,  1,  0,  0, -2,  0,  2,  0,  0,  0,  0,
     :         0,  0,  1, -1,  2,  0,  0, -1,  0,  2,  0,  0,  0,  0,
     :         0,  0,  1, -1,  2,  0,  0,  0, -2,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=201,210 ) /
     :         0,  0,  0,  0,  1,  0,  0,  1, -2,  0,  0,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  0,  1,  0,  0, -2,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0,  0, -2,  0,  0,  2,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  3, -6,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  3, -5,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  3, -5,  0,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0, -3,  4,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -3,  5,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0, -3,  5,  0,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=211,220 ) /
     :         0,  0,  2, -2,  2,  0, -3,  3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -3,  5,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2, -4,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  0,  1, -4,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  2, -4,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -2,  4,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  0, -3,  4,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -2,  4,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0, -2,  4,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -5,  8,  0,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=221,230 ) /
     :         0,  0,  2, -2,  2,  0, -5,  6,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -5,  8,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -5,  8,  0,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0, -5,  7,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -5,  8,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  5, -8,  0,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  2,  0,  0, -1,  0, -1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  0,  0, -1,  0,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0,  0, -2,  0,  1,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=231,240 ) /
     :         0,  0,  0,  0,  0,  0,  0, -6, 11,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  6,-11,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0, -1,  0,  4,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  1,  0, -4,  0,  0,  0,  0,  0,  0,
     :         2,  0,  0, -2,  1,  0, -3,  3,  0,  0,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  2,  0,  0, -2,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0,  0, -7,  9,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  4, -5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  1 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=241,250 ) /
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  2,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  2,
     :         0,  0,  2, -2,  2,  0,  0, -2,  0,  2,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  5,  0,  0,  2,
     :         0,  0,  0,  0,  1,  0,  3, -5,  0,  0,  0,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  3, -4,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0, -3,  3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  2, -4,  0,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0,  0, -4,  4,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=251,260 ) /
     :         0,  0,  1, -1,  2,  0, -5,  7,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  3, -6,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -3,  6,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  0, -4,  6,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -3,  6,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0, -3,  6,  0,  0,  0,  0,  2,
     :         0,  0, -1,  1,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  2, -3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -5,  9,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -5,  9,  0,  0,  0,  0,  1 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=261,270 ) /
     :         0,  0,  0,  0,  0,  0,  0,  5, -9,  0,  0,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  0,  1,  0, -2,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0,  0, -2,  0,  2,  0,  0,  0,  0,
     :        -2,  0,  1,  1,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :         0,  0, -2,  2,  0,  0,  3, -3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -6, 10,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0, -6, 10,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -2,  3,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -2,  3,  0,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0, -2,  2,  0,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=271,280 ) /
     :         0,  0,  0,  0,  0,  0,  2, -3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  2, -3,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  3,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  4, -8,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -4,  8,  0,  0,  0,  0,  2,
     :         0,  0, -2,  2,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -4,  7,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=281,290 ) /
     :         0,  0,  0,  0,  0,  0,  0, -4,  7,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  4, -7,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0, -2,  3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0,  0, -2,  0,  3,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -5, 10,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  1,  0, -1,  2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  4,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -3,  5,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -3,  5,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  3, -5,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=291,300 ) /
     :         0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  1, -3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -1,  2,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0, -1,  2,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -7, 11,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -7, 11,  0,  0,  0,  0,  0,  1,
     :         0,  0, -2,  2,  0,  0,  4, -4,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  2, -3,  0,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0, -4,  4,  0,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=301,310 ) /
     :         0,  0, -1,  1,  0,  0,  4, -5,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -4,  7,  0,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0, -4,  6,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -4,  7,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -4,  6,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -4,  6,  0,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0, -4,  5,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -4,  6,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  4, -6,  0,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=311,320 ) /
     :        -2,  0,  0,  2,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  1, -1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -1,  0,  5,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  1, -3,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -1,  3,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -7, 12,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -1,  1,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -1,  1,  0,  0,  0,  0,  0,  1 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=321,330 ) /
     :         0,  0,  1, -1,  1,  0, -1,  0,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  1, -2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -2,  5,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -1,  0,  4,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0, -4,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0, -1,  1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -6, 10,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -6, 10,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=331,340 ) /
     :         0,  0,  2, -2,  1,  0,  0, -3,  0,  3,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -3,  7,  0,  0,  0,  0,  2,
     :        -2,  0,  0,  2,  0,  0,  4, -4,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -5,  8,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  5, -8,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -1,  0,  3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -1,  0,  3,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0, -3,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  2, -4,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -2,  4,  0,  0,  0,  0,  0,  1 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=341,350 ) /
     :         0,  0,  1, -1,  1,  0, -2,  3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -2,  4,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -6,  9,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -6,  9,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  6, -9,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  1,  0, -2,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0, -2,  2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -4,  6,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  4, -6,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  3, -4,  0,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=351,360 ) /
     :         0,  0,  0,  0,  0,  0,  0, -1,  0,  2,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0, -2,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  1,  0, -1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -5,  9,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  3, -4,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -3,  4,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -3,  4,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  3, -4,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  3, -4,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  1,  0,  0,  2, -2,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=361,370 ) /
     :         0,  0,  0,  0,  1,  0,  0, -1,  0,  2,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  0, -3,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  1, -5,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -1,  0,  1,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0, -1,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0, -3,  5,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0, -3,  4,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  0, -2,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=371,380 ) /
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  0, -1,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0, -1,  0,  1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0, -2,  2,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -8, 14,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  2, -5,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  5, -8,  3,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  5, -8,  3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  3, -8,  3,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=381,390 ) /
     :         0,  0,  0,  0,  0,  0,  0, -3,  8, -3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0, -2,  5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -8, 12,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -8, 12,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  1, -2,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  1,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  2,  0,  0,  2,
     :         0,  0,  2, -2,  1,  0, -5,  5,  0,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=391,400 ) /
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  3, -6,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -3,  6,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0, -3,  6,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -1,  4,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -5,  7,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -5,  7,  0,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0, -5,  6,  0,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=401,410 ) /
     :         0,  0,  0,  0,  0,  0,  5, -7,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0,  0, -1,  0,  1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -1,  0,  1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0, -1,  0,  3,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  2,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -2,  6,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  1,  0,  2, -2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -6,  9,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  6, -9,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -2,  2,  0,  0,  0,  0,  0,  1 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=411,420 ) /
     :         0,  0,  1, -1,  1,  0, -2,  1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -5,  7,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  5, -7,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0, -2,  2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  4, -5,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  1, -3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -1,  3,  0,  0,  0,  0,  0,  1 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=421,430 ) /
     :         0,  0,  1, -1,  1,  0, -1,  2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -1,  3,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -7, 10,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -7, 10,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  3, -3,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -4,  8,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -4,  5,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -4,  5,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  4, -5,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=431,440 ) /
     :         0,  0,  0,  0,  0,  0,  0, -2,  0,  5,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -9, 13,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -1,  5,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -2,  0,  4,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0, -4,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -2,  7,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0, -3,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=441,450 ) /
     :         0,  0,  0,  0,  0,  0, -2,  5,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0, -2,  5,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -6,  8,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -6,  8,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  6, -8,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  2,  0, -2,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -3,  9,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  5, -6,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  5, -6,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=451,460 ) /
     :         0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -5, 10,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  4, -4,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  4, -4,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -3,  3,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  3, -3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  3, -3,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  3, -3,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  0, -3,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=461,470 ) /
     :         0,  0,  0,  0,  0,  0,  0, -5, 13,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0, -1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0, -1,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  0, -2,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  0, -2,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  3, -2,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  3, -2,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  0, -1,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -6, 15,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -8, 15,  0,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=471,480 ) /
     :         0,  0,  0,  0,  0,  0, -3,  9, -4,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  2, -5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -2,  8, -1, -5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  6, -8,  3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=481,490 ) /
     :         0,  0,  0,  0,  0,  0,  0, -6, 16, -4, -5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -2,  8, -3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -2,  8, -3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  6, -8,  1,  5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  3, -5,  4,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -8, 11,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -8, 11,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0, -8, 11,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, 11,  0,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=491,500 ) /
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  1,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  3, -3,  0,  2,  0,  0,  0,  2,
     :         0,  0,  2, -2,  1,  0,  0,  4, -8,  3,  0,  0,  0,  0,
     :         0,  0,  1, -1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0,  0, -4,  8, -3,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  1,  2,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  1,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -3,  7,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  0,  4,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -5,  6,  0,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=501,510 ) /
     :         0,  0,  0,  0,  0,  0, -5,  6,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  5, -6,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  5, -6,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  2,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -1,  6,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  7, -9,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  2, -1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  2, -1,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  6, -7,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  5, -5,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=511,520 ) /
     :         0,  0,  0,  0,  0,  0, -1,  4,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0, -1,  4,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -7,  9,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -7,  9,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  4, -3,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  3, -1,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -4,  4,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  4, -4,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  4, -4,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  4, -4,  0,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=521,530 ) /
     :         0,  0,  0,  0,  0,  0,  0,  2,  1,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -3,  0,  5,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -9, 12,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  3,  0, -4,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0,  1, -1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  7, -8,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  3,  0, -3,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=531,540 ) /
     :         0,  0,  0,  0,  0,  0,  0,  3,  0, -3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -2,  6,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -6,  7,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  6, -7,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  6, -6,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  3,  0, -2,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  3,  0, -2,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  5, -4,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  3, -2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  3, -2,  0,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=541,550 ) /
     :         0,  0,  0,  0,  0,  0,  0,  3,  0, -1,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  3,  0, -1,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  3,  0,  0, -2,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  4, -2,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  3,  0,  0, -1,  0,  0,  2,
     :         0,  0,  2, -2,  1,  0,  0,  1,  0, -1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -8, 16,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  3,  0,  2, -5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  7, -8,  3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -5, 16, -4, -5,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=551,560 ) /
     :         0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -1,  8, -3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -8, 10,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -8, 10,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0, -8, 10,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  2,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  3,  0,  1,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -3,  8,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -5,  5,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  5, -5,  0,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=561,570 ) /
     :         0,  0,  0,  0,  0,  0,  5, -5,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  5, -5,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  7, -7,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  7, -7,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  6, -5,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  7, -8,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  5, -3,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=571,580 ) /
     :         0,  0,  0,  0,  0,  0,  4, -3,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  1,  2,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -9, 11,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -9, 11,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  4,  0, -4,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  4,  0, -3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -6,  6,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  6, -6,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  6, -6,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  4,  0, -2,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=581,590 ) /
     :         0,  0,  0,  0,  0,  0,  0,  6, -4,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  3, -1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  3, -1,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  3, -1,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  4,  0, -1,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  4,  0,  0, -2,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  5, -2,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  4,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  8, -9,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  5, -4,  0,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=591,600 ) /
     :         0,  0,  0,  0,  0,  0,  2,  1,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  2,  1,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  2,  1,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0, -7,  7,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  7, -7,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  4, -2,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  4, -2,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  4, -2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  4, -2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  5,  0, -4,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=601,610 ) /
     :         0,  0,  0,  0,  0,  0,  0,  5,  0, -3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  5,  0, -2,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -8,  8,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  8, -8,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  5, -3,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  5, -3,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -9,  9,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0, -9,  9,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0, -9,  9,  0,  0,  0,  0,  0,  1 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=611,620 ) /
     :         0,  0,  0,  0,  0,  0,  9, -9,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  6, -4,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=621,630 ) /
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,
     :         1,  0,  0, -2,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,
     :         1,  0,  0, -2,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0,
     :         1,  0,  0, -2,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0,
     :         1,  0,  0, -2,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  0,  0,  0,  0,  3, -3,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,
     :        -1,  0,  0,  2,  0,  0,  0,  4, -8,  3,  0,  0,  0,  0,
     :         1,  0,  0, -2,  0,  0,  0,  4, -8,  3,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  4, -8,  3,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=631,640 ) /
     :        -1,  0,  0,  0,  0,  0,  0,  2,  0, -3,  0,  0,  0,  0,
     :        -1,  0,  0,  0,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0,
     :        -1,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  0,  2,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0,
     :         1,  0, -1,  1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  0,  2,  0,  0,  0,  2,  0, -3,  0,  0,  0,  0,
     :        -2,  0,  0,  0,  0,  0,  0,  2,  0, -3,  0,  0,  0,  0,
     :         1,  0,  0,  0,  0,  0,  0,  4, -8,  3,  0,  0,  0,  0,
     :        -1,  0,  1, -1,  1,  0,  0, -1,  0,  0,  0,  0,  0,  0,
     :         1,  0,  1, -1,  1,  0,  0, -1,  0,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=641,650 ) /
     :        -1,  0,  0,  0,  0,  0,  0,  4, -8,  3,  0,  0,  0,  0,
     :        -1,  0,  0,  2,  1,  0,  0,  2,  0, -2,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,
     :        -1,  0,  0,  2,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,
     :        -1,  0,  0,  2,  0,  0,  3, -3,  0,  0,  0,  0,  0,  0,
     :         1,  0,  0, -2,  1,  0,  0, -2,  0,  2,  0,  0,  0,  0,
     :         1,  0,  2, -2,  2,  0, -3,  3,  0,  0,  0,  0,  0,  0,
     :         1,  0,  2, -2,  2,  0,  0, -2,  0,  2,  0,  0,  0,  0,
     :         1,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0,
     :         1,  0,  0,  0,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=651,660 ) /
     :         0,  0,  0, -2,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0, -2,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0,
     :         0,  0,  2,  0,  2,  0, -2,  2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2,  0,  2,  0,  0, -1,  0,  1,  0,  0,  0,  0,
     :         0,  0,  2,  0,  2,  0, -1,  1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2,  0,  2,  0, -2,  3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  2,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,
     :         0,  0,  1,  1,  2,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :         1,  0,  2,  0,  2,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  2,  0,  2,  0, 10, -3,  0,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=661,670 ) /
     :         0,  0,  1,  1,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :         1,  0,  2,  0,  2,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2,  0,  2,  0,  0,  4, -8,  3,  0,  0,  0,  0,
     :         0,  0,  2,  0,  2,  0,  0, -4,  8, -3,  0,  0,  0,  0,
     :        -1,  0,  2,  0,  2,  0,  0, -4,  8, -3,  0,  0,  0,  0,
     :         2,  0,  2, -2,  2,  0,  0, -2,  0,  3,  0,  0,  0,  0,
     :         1,  0,  2,  0,  1,  0,  0, -2,  0,  3,  0,  0,  0,  0,
     :         0,  0,  1,  1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  2,  0,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :        -2,  0,  2,  2,  2,  0,  0,  2,  0, -2,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=671,680 ) /
     :         0,  0,  2,  0,  2,  0,  2, -3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2,  0,  2,  0,  1, -1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2,  0,  2,  0,  0,  1,  0, -1,  0,  0,  0,  0,
     :         0,  0,  2,  0,  2,  0,  2, -2,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  2,  2,  2,  0,  0, -1,  0,  1,  0,  0,  0,  0,
     :         1,  0,  2,  0,  2,  0, -1,  1,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  2,  2,  2,  0,  0,  2,  0, -3,  0,  0,  0,  0,
     :         2,  0,  2,  0,  2,  0,  0,  2,  0, -3,  0,  0,  0,  0,
     :         1,  0,  2,  0,  2,  0,  0, -4,  8, -3,  0,  0,  0,  0,
     :         1,  0,  2,  0,  2,  0,  0,  4, -8,  3,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=681,687 ) /
     :         1,  0,  1,  1,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2,  0,  2,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :         2,  0,  2,  0,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  2,  2,  2,  0,  0,  2,  0, -2,  0,  0,  0,  0,
     :        -1,  0,  2,  2,  2,  0,  3, -3,  0,  0,  0,  0,  0,  0,
     :         1,  0,  2,  0,  2,  0,  1, -1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2,  2,  2,  0,  0,  2,  0, -2,  0,  0,  0,  0 /

*
*  Planetary nutation coefficients, unit 1e-7 arcsec
*  longitude (sin, cos), obliquity (sin, cos)
*

      DATA ( ( ICPL(I,J), I=1,4 ), J=  1, 10 ) /
     :       1440,          0,          0,          0,
     :         56,       -117,        -42,        -40,
     :        125,        -43,          0,        -54,
     :          0,          5,          0,          0,
     :          3,         -7,         -3,          0,
     :          3,          0,          0,         -2,
     :       -114,          0,          0,         61,
     :       -219,         89,          0,          0,
     :         -3,          0,          0,          0,
     :       -462,       1604,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J= 11, 20 ) /
     :         99,          0,          0,        -53,
     :         -3,          0,          0,          2,
     :          0,          6,          2,          0,
     :          3,          0,          0,          0,
     :        -12,          0,          0,          0,
     :         14,       -218,        117,          8,
     :         31,       -481,       -257,        -17,
     :       -491,        128,          0,          0,
     :      -3084,       5123,       2735,       1647,
     :      -1444,       2409,      -1286,       -771 /
      DATA ( ( ICPL(I,J), I=1,4 ), J= 21, 30 ) /
     :         11,        -24,        -11,         -9,
     :         26,         -9,          0,          0,
     :        103,        -60,          0,          0,
     :          0,        -13,         -7,          0,
     :        -26,        -29,        -16,         14,
     :          9,        -27,        -14,         -5,
     :         12,          0,          0,         -6,
     :         -7,          0,          0,          0,
     :          0,         24,          0,          0,
     :        284,          0,          0,       -151 /
      DATA ( ( ICPL(I,J), I=1,4 ), J= 31, 40 ) /
     :        226,        101,          0,          0,
     :          0,         -8,         -2,          0,
     :          0,         -6,         -3,          0,
     :          5,          0,          0,         -3,
     :        -41,        175,         76,         17,
     :          0,         15,          6,          0,
     :        425,        212,       -133,        269,
     :       1200,        598,        319,       -641,
     :        235,        334,          0,          0,
     :         11,        -12,         -7,         -6 /
      DATA ( ( ICPL(I,J), I=1,4 ), J= 41, 50 ) /
     :          5,         -6,          3,          3,
     :         -5,          0,          0,          3,
     :          6,          0,          0,         -3,
     :         15,          0,          0,          0,
     :         13,          0,          0,         -7,
     :         -6,         -9,          0,          0,
     :        266,        -78,          0,          0,
     :       -460,       -435,       -232,        246,
     :          0,         15,          7,          0,
     :         -3,          0,          0,          2 /
      DATA ( ( ICPL(I,J), I=1,4 ), J= 51, 60 ) /
     :          0,        131,          0,          0,
     :          4,          0,          0,          0,
     :          0,          3,          0,          0,
     :          0,          4,          2,          0,
     :          0,          3,          0,          0,
     :        -17,        -19,        -10,          9,
     :         -9,        -11,          6,         -5,
     :         -6,          0,          0,          3,
     :        -16,          8,          0,          0,
     :          0,          3,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J= 61, 70 ) /
     :         11,         24,         11,         -5,
     :         -3,         -4,         -2,          1,
     :          3,          0,          0,         -1,
     :          0,         -8,         -4,          0,
     :          0,          3,          0,          0,
     :          0,          5,          0,          0,
     :          0,          3,          2,          0,
     :         -6,          4,          2,          3,
     :         -3,         -5,          0,          0,
     :         -5,          0,          0,          2 /
      DATA ( ( ICPL(I,J), I=1,4 ), J= 71, 80 ) /
     :          4,         24,         13,         -2,
     :        -42,         20,          0,          0,
     :        -10,        233,          0,          0,
     :         -3,          0,          0,          1,
     :         78,        -18,          0,          0,
     :          0,          3,          1,          0,
     :          0,         -3,         -1,          0,
     :          0,         -4,         -2,          1,
     :          0,         -8,         -4,         -1,
     :          0,         -5,          3,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J= 81, 90 ) /
     :         -7,          0,          0,          3,
     :        -14,          8,          3,          6,
     :          0,          8,         -4,          0,
     :          0,         19,         10,          0,
     :         45,        -22,          0,          0,
     :         -3,          0,          0,          0,
     :          0,         -3,          0,          0,
     :          0,          3,          0,          0,
     :          3,          5,          3,         -2,
     :         89,        -16,         -9,        -48 /
      DATA ( ( ICPL(I,J), I=1,4 ), J= 91,100 ) /
     :          0,          3,          0,          0,
     :         -3,          7,          4,          2,
     :       -349,        -62,          0,          0,
     :        -15,         22,          0,          0,
     :         -3,          0,          0,          0,
     :        -53,          0,          0,          0,
     :          5,          0,          0,         -3,
     :          0,         -8,          0,          0,
     :         15,         -7,         -4,         -8,
     :         -3,          0,          0,          1 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=101,110 ) /
     :        -21,        -78,          0,          0,
     :         20,        -70,        -37,        -11,
     :          0,          6,          3,          0,
     :          5,          3,          2,         -2,
     :        -17,         -4,         -2,          9,
     :          0,          6,          3,          0,
     :         32,         15,         -8,         17,
     :        174,         84,         45,        -93,
     :         11,         56,          0,          0,
     :        -66,        -12,         -6,         35 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=111,120 ) /
     :         47,          8,          4,        -25,
     :          0,          8,          4,          0,
     :         10,        -22,        -12,         -5,
     :         -3,          0,          0,          2,
     :        -24,         12,          0,          0,
     :          5,         -6,          0,          0,
     :          3,          0,          0,         -2,
     :          4,          3,          1,         -2,
     :          0,         29,         15,          0,
     :         -5,         -4,         -2,          2 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=121,130 ) /
     :          8,         -3,         -1,         -5,
     :          0,         -3,          0,          0,
     :         10,          0,          0,          0,
     :          3,          0,          0,         -2,
     :         -5,          0,          0,          3,
     :         46,         66,         35,        -25,
     :        -14,          7,          0,          0,
     :          0,          3,          2,          0,
     :         -5,          0,          0,          0,
     :        -68,        -34,        -18,         36 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=131,140 ) /
     :          0,         14,          7,          0,
     :         10,         -6,         -3,         -5,
     :         -5,         -4,         -2,          3,
     :         -3,          5,          2,          1,
     :         76,         17,          9,        -41,
     :         84,        298,        159,        -45,
     :          3,          0,          0,         -1,
     :         -3,          0,          0,          2,
     :         -3,          0,          0,          1,
     :        -82,        292,        156,         44 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=141,150 ) /
     :        -73,         17,          9,         39,
     :         -9,        -16,          0,          0,
     :          3,          0,         -1,         -2,
     :         -3,          0,          0,          0,
     :         -9,         -5,         -3,          5,
     :       -439,          0,          0,          0,
     :         57,        -28,        -15,        -30,
     :          0,         -6,         -3,          0,
     :         -4,          0,          0,          2,
     :        -40,         57,         30,         21 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=151,160 ) /
     :         23,          7,          3,        -13,
     :        273,         80,         43,       -146,
     :       -449,        430,          0,          0,
     :         -8,        -47,        -25,          4,
     :          6,         47,         25,         -3,
     :          0,         23,         13,          0,
     :         -3,          0,          0,          2,
     :          3,         -4,         -2,         -2,
     :        -48,       -110,        -59,         26,
     :         51,        114,         61,        -27 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=161,170 ) /
     :       -133,          0,          0,         57,
     :          0,          4,          0,          0,
     :        -21,         -6,         -3,         11,
     :          0,         -3,         -1,          0,
     :        -11,        -21,        -11,          6,
     :        -18,       -436,       -233,          9,
     :         35,         -7,          0,          0,
     :          0,          5,          3,          0,
     :         11,         -3,         -1,         -6,
     :         -5,         -3,         -1,          3 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=171,180 ) /
     :        -53,         -9,         -5,         28,
     :          0,          3,          2,          1,
     :          4,          0,          0,         -2,
     :          0,         -4,          0,          0,
     :        -50,        194,        103,         27,
     :        -13,         52,         28,          7,
     :        -91,        248,          0,          0,
     :          6,         49,         26,         -3,
     :         -6,        -47,        -25,          3,
     :          0,          5,          3,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=181,190 ) /
     :         52,         23,         10,        -23,
     :         -3,          0,          0,          1,
     :          0,          5,          3,          0,
     :         -4,          0,          0,          0,
     :         -4,          8,          3,          2,
     :         10,          0,          0,          0,
     :          3,          0,          0,         -2,
     :          0,          8,          4,          0,
     :          0,          8,          4,          1,
     :         -4,          0,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=191,200 ) /
     :         -4,          0,          0,          0,
     :         -8,          4,          2,          4,
     :          8,         -4,         -2,         -4,
     :          0,         15,          7,          0,
     :       -138,          0,          0,          0,
     :          0,         -7,         -3,          0,
     :          0,         -7,         -3,          0,
     :         54,          0,          0,        -29,
     :          0,         10,          4,          0,
     :         -7,          0,          0,          3 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=201,210 ) /
     :        -37,         35,         19,         20,
     :          0,          4,          0,          0,
     :         -4,          9,          0,          0,
     :          8,          0,          0,         -4,
     :         -9,        -14,         -8,          5,
     :         -3,         -9,         -5,          3,
     :       -145,         47,          0,          0,
     :        -10,         40,         21,          5,
     :         11,        -49,        -26,         -7,
     :      -2150,          0,          0,        932 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=211,220 ) /
     :        -12,          0,          0,          5,
     :         85,          0,          0,        -37,
     :          4,          0,          0,         -2,
     :          3,          0,          0,         -2,
     :        -86,        153,          0,          0,
     :         -6,          9,          5,          3,
     :          9,        -13,         -7,         -5,
     :         -8,         12,          6,          4,
     :        -51,          0,          0,         22,
     :        -11,       -268,       -116,          5 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=221,230 ) /
     :          0,         12,          5,          0,
     :          0,          7,          3,          0,
     :         31,          6,          3,        -17,
     :        140,         27,         14,        -75,
     :         57,         11,          6,        -30,
     :        -14,        -39,          0,          0,
     :          0,         -6,         -2,          0,
     :          4,         15,          8,         -2,
     :          0,          4,          0,          0,
     :         -3,          0,          0,          1 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=231,240 ) /
     :          0,         11,          5,          0,
     :          9,          6,          0,          0,
     :         -4,         10,          4,          2,
     :          5,          3,          0,          0,
     :         16,          0,          0,         -9,
     :         -3,          0,          0,          0,
     :          0,          3,          2,         -1,
     :          7,          0,          0,         -3,
     :        -25,         22,          0,          0,
     :         42,        223,        119,        -22 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=241,250 ) /
     :        -27,       -143,        -77,         14,
     :          9,         49,         26,         -5,
     :      -1166,          0,          0,        505,
     :         -5,          0,          0,          2,
     :         -6,          0,          0,          3,
     :         -8,          0,          1,          4,
     :          0,         -4,          0,          0,
     :        117,          0,          0,        -63,
     :         -4,          8,          4,          2,
     :          3,          0,          0,         -2 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=251,260 ) /
     :         -5,          0,          0,          2,
     :          0,         31,          0,          0,
     :         -5,          0,          1,          3,
     :          4,          0,          0,         -2,
     :         -4,          0,          0,          2,
     :        -24,        -13,         -6,         10,
     :          3,          0,          0,          0,
     :          0,        -32,        -17,          0,
     :          8,         12,          5,         -3,
     :          3,          0,          0,         -1 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=261,270 ) /
     :          7,         13,          0,          0,
     :         -3,         16,          0,          0,
     :         50,          0,          0,        -27,
     :          0,         -5,         -3,          0,
     :         13,          0,          0,          0,
     :          0,          5,          3,          1,
     :         24,          5,          2,        -11,
     :          5,        -11,         -5,         -2,
     :         30,         -3,         -2,        -16,
     :         18,          0,          0,         -9 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=271,280 ) /
     :          8,        614,          0,          0,
     :          3,         -3,         -1,         -2,
     :          6,         17,          9,         -3,
     :         -3,         -9,         -5,          2,
     :          0,          6,          3,         -1,
     :       -127,         21,          9,         55,
     :          3,          5,          0,          0,
     :         -6,        -10,         -4,          3,
     :          5,          0,          0,          0,
     :         16,          9,          4,         -7 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=281,290 ) /
     :          3,          0,          0,         -2,
     :          0,         22,          0,          0,
     :          0,         19,         10,          0,
     :          7,          0,          0,         -4,
     :          0,         -5,         -2,          0,
     :          0,          3,          1,          0,
     :         -9,          3,          1,          4,
     :         17,          0,          0,         -7,
     :          0,         -3,         -2,         -1,
     :        -20,         34,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=291,300 ) /
     :        -10,          0,          1,          5,
     :         -4,          0,          0,          2,
     :         22,        -87,          0,          0,
     :         -4,          0,          0,          2,
     :         -3,         -6,         -2,          1,
     :        -16,         -3,         -1,          7,
     :          0,         -3,         -2,          0,
     :          4,          0,          0,          0,
     :        -68,         39,          0,          0,
     :         27,          0,          0,        -14 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=301,310 ) /
     :          0,         -4,          0,          0,
     :        -25,          0,          0,          0,
     :        -12,         -3,         -2,          6,
     :          3,          0,          0,         -1,
     :          3,         66,         29,         -1,
     :        490,          0,          0,       -213,
     :        -22,         93,         49,         12,
     :         -7,         28,         15,          4,
     :         -3,         13,          7,          2,
     :        -46,         14,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=311,320 ) /
     :         -5,          0,          0,          0,
     :          2,          1,          0,          0,
     :          0,         -3,          0,          0,
     :        -28,          0,          0,         15,
     :          5,          0,          0,         -2,
     :          0,          3,          0,          0,
     :        -11,          0,          0,          5,
     :          0,          3,          1,          0,
     :         -3,          0,          0,          1,
     :         25,        106,         57,        -13 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=321,330 ) /
     :          5,         21,         11,         -3,
     :       1485,          0,          0,          0,
     :         -7,        -32,        -17,          4,
     :          0,          5,          3,          0,
     :         -6,         -3,         -2,          3,
     :         30,         -6,         -2,        -13,
     :         -4,          4,          0,          0,
     :        -19,          0,          0,         10,
     :          0,          4,          2,         -1,
     :          0,          3,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=331,340 ) /
     :          4,          0,          0,         -2,
     :          0,         -3,         -1,          0,
     :         -3,          0,          0,          0,
     :          5,          3,          1,         -2,
     :          0,         11,          0,          0,
     :        118,          0,          0,        -52,
     :          0,         -5,         -3,          0,
     :        -28,         36,          0,          0,
     :          5,         -5,          0,          0,
     :         14,        -59,        -31,         -8 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=341,350 ) /
     :          0,          9,          5,          1,
     :       -458,          0,          0,        198,
     :          0,        -45,        -20,          0,
     :          9,          0,          0,         -5,
     :          0,         -3,          0,          0,
     :          0,         -4,         -2,         -1,
     :         11,          0,          0,         -6,
     :          6,          0,          0,         -2,
     :        -16,         23,          0,          0,
     :          0,         -4,         -2,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=351,360 ) /
     :         -5,          0,          0,          2,
     :       -166,        269,          0,          0,
     :         15,          0,          0,         -8,
     :         10,          0,          0,         -4,
     :        -78,         45,          0,          0,
     :          0,         -5,         -2,          0,
     :          7,          0,          0,         -4,
     :         -5,        328,          0,          0,
     :          3,          0,          0,         -2,
     :          5,          0,          0,         -2 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=361,370 ) /
     :          0,          3,          1,          0,
     :         -3,          0,          0,          0,
     :         -3,          0,          0,          0,
     :          0,         -4,         -2,          0,
     :      -1223,        -26,          0,          0,
     :          0,          7,          3,          0,
     :          3,          0,          0,          0,
     :          0,          3,          2,          0,
     :         -6,         20,          0,          0,
     :       -368,          0,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=371,380 ) /
     :        -75,          0,          0,          0,
     :         11,          0,          0,         -6,
     :          3,          0,          0,         -2,
     :         -3,          0,          0,          1,
     :        -13,        -30,          0,          0,
     :         21,          3,          0,          0,
     :         -3,          0,          0,          1,
     :         -4,          0,          0,          2,
     :          8,        -27,          0,          0,
     :        -19,        -11,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=381,390 ) /
     :         -4,          0,          0,          2,
     :          0,          5,          2,          0,
     :         -6,          0,          0,          2,
     :         -8,          0,          0,          0,
     :         -1,          0,          0,          0,
     :        -14,          0,          0,          6,
     :          6,          0,          0,          0,
     :        -74,          0,          0,         32,
     :          0,         -3,         -1,          0,
     :          4,          0,          0,         -2 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=391,400 ) /
     :          8,         11,          0,          0,
     :          0,          3,          2,          0,
     :       -262,          0,          0,        114,
     :          0,         -4,          0,          0,
     :         -7,          0,          0,          4,
     :          0,        -27,        -12,          0,
     :        -19,         -8,         -4,          8,
     :        202,          0,          0,        -87,
     :         -8,         35,         19,          5,
     :          0,          4,          2,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=401,410 ) /
     :         16,         -5,          0,          0,
     :          5,          0,          0,         -3,
     :          0,         -3,          0,          0,
     :          1,          0,          0,          0,
     :        -35,        -48,        -21,         15,
     :         -3,         -5,         -2,          1,
     :          6,          0,          0,         -3,
     :          3,          0,          0,         -1,
     :          0,         -5,          0,          0,
     :         12,         55,         29,         -6 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=411,420 ) /
     :          0,          5,          3,          0,
     :       -598,          0,          0,          0,
     :         -3,        -13,         -7,          1,
     :         -5,         -7,         -3,          2,
     :          3,          0,          0,         -1,
     :          5,         -7,          0,          0,
     :          4,          0,          0,         -2,
     :         16,         -6,          0,          0,
     :          8,         -3,          0,          0,
     :          8,        -31,        -16,         -4 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=421,430 ) /
     :          0,          3,          1,          0,
     :        113,          0,          0,        -49,
     :          0,        -24,        -10,          0,
     :          4,          0,          0,         -2,
     :         27,          0,          0,          0,
     :         -3,          0,          0,          1,
     :          0,         -4,         -2,          0,
     :          5,          0,          0,         -2,
     :          0,         -3,          0,          0,
     :        -13,          0,          0,          6 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=431,440 ) /
     :          5,          0,          0,         -2,
     :        -18,        -10,         -4,          8,
     :         -4,        -28,          0,          0,
     :         -5,          6,          3,          2,
     :         -3,          0,          0,          1,
     :         -5,         -9,         -4,          2,
     :         17,          0,          0,         -7,
     :         11,          4,          0,          0,
     :          0,         -6,         -2,          0,
     :         83,         15,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=441,450 ) /
     :         -4,          0,          0,          2,
     :          0,       -114,        -49,          0,
     :        117,          0,          0,        -51,
     :         -5,         19,         10,          2,
     :         -3,          0,          0,          0,
     :         -3,          0,          0,          2,
     :          0,         -3,         -1,          0,
     :          3,          0,          0,          0,
     :          0,         -6,         -2,          0,
     :        393,          3,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=451,460 ) /
     :         -4,         21,         11,          2,
     :         -6,          0,         -1,          3,
     :         -3,          8,          4,          1,
     :          8,          0,          0,          0,
     :         18,        -29,        -13,         -8,
     :          8,         34,         18,         -4,
     :         89,          0,          0,          0,
     :          3,         12,          6,         -1,
     :         54,        -15,         -7,        -24,
     :          0,          3,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=461,470 ) /
     :          3,          0,          0,         -1,
     :          0,         35,          0,          0,
     :       -154,        -30,        -13,         67,
     :         15,          0,          0,          0,
     :          0,          4,          2,          0,
     :          0,          9,          0,          0,
     :         80,        -71,        -31,        -35,
     :          0,        -20,         -9,          0,
     :         11,          5,          2,         -5,
     :         61,        -96,        -42,        -27 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=471,480 ) /
     :         14,          9,          4,         -6,
     :        -11,         -6,         -3,          5,
     :          0,         -3,         -1,          0,
     :        123,       -415,       -180,        -53,
     :          0,          0,          0,        -35,
     :         -5,          0,          0,          0,
     :          7,        -32,        -17,         -4,
     :          0,         -9,         -5,          0,
     :          0,         -4,          2,          0,
     :        -89,          0,          0,         38 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=481,490 ) /
     :          0,        -86,        -19,         -6,
     :          0,          0,        -19,          6,
     :       -123,       -416,       -180,         53,
     :          0,         -3,         -1,          0,
     :         12,         -6,         -3,         -5,
     :        -13,          9,          4,          6,
     :          0,        -15,         -7,          0,
     :          3,          0,          0,         -1,
     :        -62,        -97,        -42,         27,
     :        -11,          5,          2,          5 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=491,500 ) /
     :          0,        -19,         -8,          0,
     :         -3,          0,          0,          1,
     :          0,          4,          2,          0,
     :          0,          3,          0,          0,
     :          0,          4,          2,          0,
     :        -85,        -70,        -31,         37,
     :        163,        -12,         -5,        -72,
     :        -63,        -16,         -7,         28,
     :        -21,        -32,        -14,          9,
     :          0,         -3,         -1,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=501,510 ) /
     :          3,          0,          0,         -2,
     :          0,          8,          0,          0,
     :          3,         10,          4,         -1,
     :          3,          0,          0,         -1,
     :          0,         -7,         -3,          0,
     :          0,         -4,         -2,          0,
     :          6,         19,          0,          0,
     :          5,       -173,        -75,         -2,
     :          0,         -7,         -3,          0,
     :          7,        -12,         -5,         -3 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=511,520 ) /
     :         -3,          0,          0,          2,
     :          3,         -4,         -2,         -1,
     :         74,          0,          0,        -32,
     :         -3,         12,          6,          2,
     :         26,        -14,         -6,        -11,
     :         19,          0,          0,         -8,
     :          6,         24,         13,         -3,
     :         83,          0,          0,          0,
     :          0,        -10,         -5,          0,
     :         11,         -3,         -1,         -5 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=521,530 ) /
     :          3,          0,          1,         -1,
     :          3,          0,          0,         -1,
     :         -4,          0,          0,          0,
     :          5,        -23,        -12,         -3,
     :       -339,          0,          0,        147,
     :          0,        -10,         -5,          0,
     :          5,          0,          0,          0,
     :          3,          0,          0,         -1,
     :          0,         -4,         -2,          0,
     :         18,         -3,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=531,540 ) /
     :          9,        -11,         -5,         -4,
     :         -8,          0,          0,          4,
     :          3,          0,          0,         -1,
     :          0,          9,          0,          0,
     :          6,         -9,         -4,         -2,
     :         -4,        -12,          0,          0,
     :         67,        -91,        -39,        -29,
     :         30,        -18,         -8,        -13,
     :          0,          0,          0,          0,
     :          0,       -114,        -50,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=541,550 ) /
     :          0,          0,          0,         23,
     :        517,         16,          7,       -224,
     :          0,         -7,         -3,          0,
     :        143,         -3,         -1,        -62,
     :         29,          0,          0,        -13,
     :         -4,          0,          0,          2,
     :         -6,          0,          0,          3,
     :          5,         12,          5,         -2,
     :        -25,          0,          0,         11,
     :         -3,          0,          0,          1 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=551,560 ) /
     :          0,          4,          2,          0,
     :        -22,         12,          5,         10,
     :         50,          0,          0,        -22,
     :          0,          7,          4,          0,
     :          0,          3,          1,          0,
     :         -4,          4,          2,          2,
     :         -5,        -11,         -5,          2,
     :          0,          4,          2,          0,
     :          4,         17,          9,         -2,
     :         59,          0,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=561,570 ) /
     :          0,         -4,         -2,          0,
     :         -8,          0,          0,          4,
     :         -3,          0,          0,          0,
     :          4,        -15,         -8,         -2,
     :        370,         -8,          0,       -160,
     :          0,          0,         -3,          0,
     :          0,          3,          1,          0,
     :         -6,          3,          1,          3,
     :          0,          6,          0,          0,
     :        -10,          0,          0,          4 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=571,580 ) /
     :          0,          9,          4,          0,
     :          4,         17,          7,         -2,
     :         34,          0,          0,        -15,
     :          0,          5,          3,          0,
     :         -5,          0,          0,          2,
     :        -37,         -7,         -3,         16,
     :          3,         13,          7,         -2,
     :         40,          0,          0,          0,
     :          0,         -3,         -2,          0,
     :       -184,         -3,         -1,         80 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=581,590 ) /
     :         -3,          0,          0,          1,
     :         -3,          0,          0,          0,
     :          0,        -10,         -6,         -1,
     :         31,         -6,          0,        -13,
     :         -3,        -32,        -14,          1,
     :         -7,          0,          0,          3,
     :          0,         -8,         -4,          0,
     :          3,         -4,          0,          0,
     :          0,          4,          0,          0,
     :          0,          3,          1,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=591,600 ) /
     :         19,        -23,        -10,          2,
     :          0,          0,          0,        -10,
     :          0,          3,          2,          0,
     :          0,          9,          5,         -1,
     :         28,          0,          0,          0,
     :          0,         -7,         -4,          0,
     :          8,         -4,          0,         -4,
     :          0,          0,         -2,          0,
     :          0,          3,          0,          0,
     :         -3,          0,          0,          1 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=601,610 ) /
     :         -9,          0,          1,          4,
     :          3,         12,          5,         -1,
     :         17,         -3,         -1,          0,
     :          0,          7,          4,          0,
     :         19,          0,          0,          0,
     :          0,         -5,         -3,          0,
     :         14,         -3,          0,         -1,
     :          0,          0,         -1,          0,
     :          0,          0,          0,         -5,
     :          0,          5,          3,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=611,620 ) /
     :         13,          0,          0,          0,
     :          0,         -3,         -2,          0,
     :          2,          9,          4,          3,
     :          0,          0,          0,         -4,
     :          8,          0,          0,          0,
     :          0,          4,          2,          0,
     :          6,          0,          0,         -3,
     :          6,          0,          0,          0,
     :          0,          3,          1,          0,
     :          5,          0,          0,         -2 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=621,630 ) /
     :          3,          0,          0,         -1,
     :         -3,          0,          0,          0,
     :          6,          0,          0,          0,
     :          7,          0,          0,          0,
     :         -4,          0,          0,          0,
     :          4,          0,          0,          0,
     :          6,          0,          0,          0,
     :          0,         -4,          0,          0,
     :          0,         -4,          0,          0,
     :          5,          0,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=631,640 ) /
     :         -3,          0,          0,          0,
     :          4,          0,          0,          0,
     :         -5,          0,          0,          0,
     :          4,          0,          0,          0,
     :          0,          3,          0,          0,
     :         13,          0,          0,          0,
     :         21,         11,          0,          0,
     :          0,         -5,          0,          0,
     :          0,         -5,         -2,          0,
     :          0,          5,          3,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=641,650 ) /
     :          0,         -5,          0,          0,
     :         -3,          0,          0,          2,
     :         20,         10,          0,          0,
     :        -34,          0,          0,          0,
     :        -19,          0,          0,          0,
     :          3,          0,          0,         -2,
     :         -3,          0,          0,          1,
     :         -6,          0,          0,          3,
     :         -4,          0,          0,          0,
     :          3,          0,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=651,660 ) /
     :          3,          0,          0,          0,
     :          4,          0,          0,          0,
     :          3,          0,          0,         -1,
     :          6,          0,          0,         -3,
     :         -8,          0,          0,          3,
     :          0,          3,          1,          0,
     :         -3,          0,          0,          0,
     :          0,         -3,         -2,          0,
     :        126,        -63,        -27,        -55,
     :         -5,          0,          1,          2 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=661,670 ) /
     :         -3,         28,         15,          2,
     :          5,          0,          1,         -2,
     :          0,          9,          4,          1,
     :          0,          9,          4,         -1,
     :       -126,        -63,        -27,         55,
     :          3,          0,          0,         -1,
     :         21,        -11,         -6,        -11,
     :          0,         -4,          0,          0,
     :        -21,        -11,         -6,         11,
     :         -3,          0,          0,          1 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=671,680 ) /
     :          0,          3,          1,          0,
     :          8,          0,          0,         -4,
     :         -6,          0,          0,          3,
     :         -3,          0,          0,          1,
     :          3,          0,          0,         -1,
     :         -3,          0,          0,          1,
     :         -5,          0,          0,          2,
     :         24,        -12,         -5,        -11,
     :          0,          3,          1,          0,
     :          0,          3,          1,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=681,687 ) /
     :          0,          3,          2,          0,
     :        -24,        -12,         -5,         10,
     :          4,          0,         -1,         -2,
     :         13,          0,          0,         -6,
     :          7,          0,          0,         -3,
     :          3,          0,          0,         -1,
     :          3,          0,          0,         -1 /

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Interval between fundamental date J2000.0 and given date (JC).
      T = ( ( DATE1-DJ00 ) + DATE2 ) / DJC

*  -------------------
*  LUNI-SOLAR NUTATION
*  -------------------

*
*  Fundamental (Delaunay) arguments
*

*  Mean anomaly of the Moon (IERS 2003).
      EL = iau_FAL03 ( T )

*  Mean anomaly of the Sun (MHB2000).
      ELP = MOD (       1287104.79305D0 +
     :            T*( 129596581.0481D0 +
     :            T*(       - 0.5532D0 +
     :            T*(         0.000136D0 +
     :            T*(       - 0.00001149D0 )))), TURNAS ) * DAS2R

*  Mean longitude of the Moon minus that of the ascending node
*  (IERS 2003.
      F = iau_FAF03 ( T )

*  Mean elongation of the Moon from the Sun (MHB2000).
      D = MOD (        1072260.70369D0 +
     :          T*( 1602961601.2090D0 +
     :          T*(        - 6.3706D0 +
     :          T*(          0.006593D0 +
     :          T*(        - 0.00003169D0 )))), TURNAS ) * DAS2R

*  Mean longitude of the ascending node of the Moon (IERS 2003).
      OM = iau_FAOM03 ( T )

*  Initialize the nutation values.
      DP = 0D0
      DE = 0D0

*  Summation of luni-solar nutation series (in reverse order).
      DO 100 I = NLS, 1, -1

*     Argument and functions.
         ARG = MOD ( DBLE ( NALS(1,I) ) * EL  +
     :               DBLE ( NALS(2,I) ) * ELP +
     :               DBLE ( NALS(3,I) ) * F   +
     :               DBLE ( NALS(4,I) ) * D   +
     :               DBLE ( NALS(5,I) ) * OM, D2PI )
         SARG = SIN(ARG)
         CARG = COS(ARG)

*     Term.
         DP = DP + ( CLS(1,I) + CLS(2,I) * T ) * SARG
     :           +   CLS(3,I)                  * CARG
         DE = DE + ( CLS(4,I) + CLS(5,I) * T ) * CARG
     :           +   CLS(6,I)                  * SARG

 100  CONTINUE

*  Convert from 0.1 microarcsec units to radians.
      DPSILS = DP * U2R
      DEPSLS = DE * U2R

*  ------------------
*  PLANETARY NUTATION
*  ------------------

*  n.b.  The MHB2000 code computes the luni-solar and planetary nutation
*        in different routines, using slightly different Delaunay
*        arguments in the two cases.  This behaviour is faithfully
*        reproduced here.  Use of the IERS 2003 expressions for both
*        cases leads to negligible changes, well below
*        0.1 microarcsecond.

*  Mean anomaly of the Moon (MHB2000).
      AL = MOD ( 2.35555598D0 + 8328.6914269554D0 * T, D2PI )

*  Mean anomaly of the Sun (MHB2000).
      ALSU = MOD ( 6.24006013D0 + 628.301955D0 * T, D2PI )

*  Mean longitude of the Moon minus that of the ascending node
* (MHB2000).
      AF = MOD ( 1.627905234D0 + 8433.466158131D0 * T, D2PI )

*  Mean elongation of the Moon from the Sun (MHB2000).
      AD = MOD ( 5.198466741D0 + 7771.3771468121D0 * T, D2PI )

*  Mean longitude of the ascending node of the Moon (MHB2000).
      AOM = MOD ( 2.18243920D0 - 33.757045D0 * T, D2PI )

*  General accumulated precession in longitude (IERS 2003).
      APA = iau_FAPA03 ( T )

*  Planetary longitudes, Mercury through Uranus (IERS 2003).
      ALME = iau_FAME03 ( T )
      ALVE = iau_FAVE03 ( T )
      ALEA = iau_FAE03 ( T )
      ALMA = iau_FAMA03 ( T )
      ALJU = iau_FAJU03 ( T )
      ALSA = iau_FASA03 ( T )
      ALUR = iau_FAUR03 ( T )

*  Neptune longitude (MHB2000).
      ALNE = MOD ( 5.321159000D0 + 3.8127774000D0 * T, D2PI )

*  Initialize the nutation values.
      DP = 0D0
      DE = 0D0

*  Summation of planetary nutation series (in reverse order).
      DO 200 I = NPL, 1, -1

*     Argument and functions.
         ARG = MOD ( DBLE ( NAPL( 1,I) ) * AL   +
     :               DBLE ( NAPL( 2,I) ) * ALSU +
     :               DBLE ( NAPL( 3,I) ) * AF   +
     :               DBLE ( NAPL( 4,I) ) * AD   +
     :               DBLE ( NAPL( 5,I) ) * AOM  +
     :               DBLE ( NAPL( 6,I) ) * ALME +
     :               DBLE ( NAPL( 7,I) ) * ALVE +
     :               DBLE ( NAPL( 8,I) ) * ALEA +
     :               DBLE ( NAPL( 9,I) ) * ALMA +
     :               DBLE ( NAPL(10,I) ) * ALJU +
     :               DBLE ( NAPL(11,I) ) * ALSA +
     :               DBLE ( NAPL(12,I) ) * ALUR +
     :               DBLE ( NAPL(13,I) ) * ALNE +
     :               DBLE ( NAPL(14,I) ) * APA, D2PI )
         SARG = SIN(ARG)
         CARG = COS(ARG)

*     Term.
         DP = DP + DBLE( ICPL(1,I)) * SARG + DBLE( ICPL(2,I)) * CARG
         DE = DE + DBLE( ICPL(3,I)) * SARG + DBLE( ICPL(4,I)) * CARG

 200  CONTINUE

*  Convert from 0.1 microarcsec units to radians.
      DPSIPL = DP * U2R
      DEPSPL = DE * U2R

*  -------
*  RESULTS
*  -------

*  Add luni-solar and planetary components.
      DPSI = DPSILS + DPSIPL
      DEPS = DEPSLS + DEPSPL

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
