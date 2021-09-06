      SUBROUTINE RG_ZONT2 ( T, DUT, DLOD, DOMEGA )
*+
*  - - - - - - - - - - -
*   R G _ Z O N T 2
*  - - - - - - - - - - -
*
*  This routine is part of the International Earth Rotation and
*  Reference Systems Service (IERS) Conventions software collection.
*
*  This subroutine evaluates the effects of zonal Earth tides on the
*  rotation of the Earth.  The model used is a combination of Yoder
*  et al. (1981) elastic body tide, Wahr and Bergen (1986) inelastic
*  body tide, and Kantha et al. (1998) ocean tide models 
*  as recommended by the IERS Conventions (2010).  Refer to
*  Chapter 8 pp. 105 - 106.  The latest version of the model is located
*  at http://tai.bipm.org/iers/convupdt/convupdt_c8.html.
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
*     T           d      TT, Julian centuries since J2000 (Note 1)
*
*  Returned:
*     DUT         d      Effect on UT1 (Note 2)
*     DLOD        d      Effect on excess length of day (LOD) (Note 3)
*     DOMEGA      d      Effect on rotational speed (Note 4)
*
*  Notes:
*
*  1) Though T is strictly TDB, it is usually more convenient to use
*     TT, which makes no significant difference.  Julian centuries since
*     J2000 is (JD - 2451545.0)/36525.
*
*  2) The expression used is as adopted in IERS Conventions (2010).
*     DUT is expressed in seconds and is double precision.
*
*  3) The expression used is as adopted in IERS Conventions (2010).
*     DLOD is the excess in LOD and is expressed in seconds per day
*     and is double precision.  The phrase 'per day' is generally
*     understood, so it has been omitted commonly in speech and
*     literature.  
*     See: Stephenson, F. R., Morrison, L. V., Whitrow, G. J., 1984,
*     "Long-Term Changes in the Rotation of the Earth: 700 B. C. to
*     A. D. 1980 [and Discussion]", Phil. Trans. Roy. Soc. of London.
*     Series A, 313, pp. 47 - 70.
* 
*  4) The expression used is as adopted in IERS Conventions (2010).
*     Rotational speed is expressed in radians per second and is
*     double precision.
*  
*  Called:
*     FUNDARG      Computation of the fundamental lunisolar arguments
*
*  Test case:
*     given input: T = .07995893223819302 Julian centuries since J2000
*                  (MJD = 54465)
*     expected output: DUT    =  7.983287678576557467E-002 seconds
*                      DLOD   =  5.035331113978199288E-005 seconds / day
*                      DOMEGA = -4.249711616463017E-014 radians / second
*
*  References:
*
*     Yoder, C. F., Williams, J. G., and Parke, M. E., (1981),
*     "Tidal Variations of Earth Rotation," J. Geophys. Res., 86,
*     pp. 881 - 891.
*
*     Wahr, J. and Bergen, Z., (1986), "The effects of mantle 
*     anelasticity on nutations, Earth tides, and tidal variations
*     in rotation rate," Geophys. J. Roy. astr. Soc., 87, pp. 633 - 668.
*
*     Kantha, L. H., Stewart, J. S., and Desai, S. D., (1998), "Long-
*     period lunar fortnightly and monthly ocean tides," J. Geophys.
*     Res., 103, pp. 12639 - 12647.
*
*     Gross, R. S., (2009), "Ocean tidal effects on Earth rotation,"
*     J. Geodyn., 48(3-5), pp. 219 - 225.
* 
*     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
*     IERS Technical Note No. 36, BKG (2010)
*
*  Revisions:  
*  2008 January 18 B.E. Stetzler  Initial changes to header
*               and used 2PI instead of PI as parameter
*  2008 January 25 B.E. Stetzler Additional changes to header
*  2008 February 21 B.E. Stetzler Definition of (excess) LOD clarified
*  2008 March   12 B.E. Stetzler Applied changes to wording of notes.
*  2008 March   14 B.E. Stetzler Further changes applied to code.
*  2008 April   03 B.E. Stetzler Provided example test case
*  2009 February 11 B.E. Stetzler Updated test case due to changes made
*                                 to FUNDARG.F subroutine
*  2009 April   10 B.E. Stetzler DLOD corrected to say it is expressed
*                                in seconds per day
*  2009 May     04 B.E. Stetzler Code formatting changes based on 
*                                client recommendations
*  2009 May     07 B.E. Stetzler Updated test case due to above changes
*  2010 February 19 B.E. Stetzler Replaced Conventions 2003 recommended
*                                 model with Conventions 2010 model
*  2010 February 22 B.E. Stetzler Provided example test case
*  2010 February 23 B.E. Stetzler Updated values to two decimal places
*  2010 February 23 B.E. Stetzler Split fundamental arguments and
*                                 coefficients for four decimal place
*                                 precision
*  2010 February 25 B.E. Stetzler Recalculation of fundamental arguments
*  2010 March    01 B.E. Stetzler Updated table values to four decimal
*                                 places and double precision
*  2010 March    12 B.E. Stetzler Applied changes to wording of notes.
*  2010 March    22 B.E. Stetzler Corrected DOMEGA output for test case
*  2011 December 20 B.E. Stetzler Corrected typo in the sine coefficient
*                                 of the 365.26 period for delta LOD.
*                                 Noted by J. Lundberg. Test case output
*                                 for DLOD updated.
*-----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER I, J
      DOUBLE PRECISION T, DUT, DLOD, DOMEGA, L, LP, F, D, OM, ARG, D2PI

*  Arcseconds to radians
      DOUBLE PRECISION DAS2R
      PARAMETER ( DAS2R = 4.848136811095359935899141D-6 )

*  Arcseconds in a full circle
      DOUBLE PRECISION TURNAS
      PARAMETER ( TURNAS = 1296000D0 )

*  2Pi
      PARAMETER (D2PI= 6.283185307179586476925287D0)

*  ----------------------
*  Zonal Earth tide model
*  ----------------------

*  Number of terms in the zonal Earth tide model  
      INTEGER NZONT
      PARAMETER ( NZONT = 62 )

*  Coefficients for the fundamental arguments
      INTEGER NFUND(5,NZONT)

*  Zonal tide term coefficients
      DOUBLE PRECISION TIDE(6,NZONT)

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*  --------------------------------------------------
*  Tables of multiples of arguments and coefficients
*  --------------------------------------------------
*  Luni-Solar argument multipliers
*      l   l'  F   D OMEGA
      DATA ( ( NFUND(I,J), I=1,5 ), J= 1, 20 ) /
     : 1,  0,  2,  2,  2, 
     : 2,  0,  2,  0,  1, 
     : 2,  0,  2,  0,  2, 
     : 0,  0,  2,  2,  1, 
     : 0,  0,  2,  2,  2, 
     : 1,  0,  2,  0,  0, 
     : 1,  0,  2,  0,  1, 
     : 1,  0,  2,  0,  2, 
     : 3,  0,  0,  0,  0, 
     :-1,  0,  2,  2,  1, 
     :-1,  0,  2,  2,  2, 
     : 1,  0,  0,  2,  0, 
     : 2,  0,  2, -2,  2, 
     : 0,  1,  2,  0,  2, 
     : 0,  0,  2,  0,  0, 
     : 0,  0,  2,  0,  1, 
     : 0,  0,  2,  0,  2, 
     : 2,  0,  0,  0, -1, 
     : 2,  0,  0,  0,  0, 
     : 2,  0,  0,  0,  1 / 
                                         
      DATA ( ( NFUND(I,J), I=1,5 ), J= 21, 40 ) /
     : 0, -1,  2,  0,  2, 
     : 0,  0,  0,  2, -1, 
     : 0,  0,  0,  2,  0, 
     : 0,  0,  0,  2,  1, 
     : 0, -1,  0,  2,  0, 
     : 1,  0,  2, -2,  1, 
     : 1,  0,  2, -2,  2, 
     : 1,  1,  0,  0,  0, 
     :-1,  0,  2,  0,  0, 
     :-1,  0,  2,  0,  1, 
     :-1,  0,  2,  0,  2, 
     : 1,  0,  0,  0, -1, 
     : 1,  0,  0,  0,  0, 
     : 1,  0,  0,  0,  1, 
     : 0,  0,  0,  1,  0, 
     : 1, -1,  0,  0,  0, 
     :-1,  0,  0,  2, -1, 
     :-1,  0,  0,  2,  0, 
     :-1,  0,  0,  2,  1, 
     : 1,  0, -2,  2, -1 / 

      DATA ( ( NFUND(I,J), I=1,5 ), J= 41, 60 ) /
     :-1, -1,  0,  2,  0, 
     : 0,  2,  2, -2,  2, 
     : 0,  1,  2, -2,  1, 
     : 0,  1,  2, -2,  2, 
     : 0,  0,  2, -2,  0, 
     : 0,  0,  2, -2,  1, 
     : 0,  0,  2, -2,  2, 
     : 0,  2,  0,  0,  0, 
     : 2,  0,  0, -2, -1, 
     : 2,  0,  0, -2,  0, 
     : 2,  0,  0, -2,  1, 
     : 0, -1,  2, -2,  1, 
     : 0,  1,  0,  0, -1, 
     : 0, -1,  2, -2,  2, 
     : 0,  1,  0,  0,  0, 
     : 0,  1,  0,  0,  1, 
     : 1,  0,  0, -1,  0, 
     : 2,  0, -2,  0,  0, 
     :-2,  0,  2,  0,  1, 
     :-1,  1,  0,  1,  0  /

      DATA ( ( NFUND(I,J), I=1,5 ), J= 61, 62 ) /
     : 0,  0,  0,  0,  2, 
     : 0,  0,  0,  0,  1 /

*    Multiple of     DUT          DLOD              DOMEGA
*            sin     cos      cos      sin       cos      sin
      DATA ( ( TIDE(I,J), I=1,6 ), J = 1,20 ) /
     :    -0.0235D0,0.0000D0, 0.2617D0, 0.0000D0, -0.2209D0, 0.0000D0,
     :    -0.0404D0,0.0000D0, 0.3706D0, 0.0000D0, -0.3128D0, 0.0000D0,
     :    -0.0987D0,0.0000D0, 0.9041D0, 0.0000D0, -0.7630D0, 0.0000D0,
     :    -0.0508D0,0.0000D0, 0.4499D0, 0.0000D0, -0.3797D0, 0.0000D0,
     :    -0.1231D0,0.0000D0, 1.0904D0, 0.0000D0, -0.9203D0, 0.0000D0,
     :    -0.0385D0,0.0000D0, 0.2659D0, 0.0000D0, -0.2244D0, 0.0000D0,
     :    -0.4108D0,0.0000D0, 2.8298D0, 0.0000D0, -2.3884D0, 0.0000D0,
     :    -0.9926D0,0.0000D0, 6.8291D0, 0.0000D0, -5.7637D0, 0.0000D0,
     :    -0.0179D0,0.0000D0, 0.1222D0, 0.0000D0, -0.1031D0, 0.0000D0,
     :    -0.0818D0,0.0000D0, 0.5384D0, 0.0000D0, -0.4544D0, 0.0000D0,
     :    -0.1974D0,0.0000D0, 1.2978D0, 0.0000D0, -1.0953D0, 0.0000D0,
     :    -0.0761D0,0.0000D0, 0.4976D0, 0.0000D0, -0.4200D0, 0.0000D0,
     :     0.0216D0,0.0000D0,-0.1060D0, 0.0000D0,  0.0895D0, 0.0000D0,
     :     0.0254D0,0.0000D0,-0.1211D0, 0.0000D0,  0.1022D0, 0.0000D0,
     :    -0.2989D0,0.0000D0, 1.3804D0, 0.0000D0, -1.1650D0, 0.0000D0,
     :    -3.1873D0,0.2010D0,14.6890D0, 0.9266D0,-12.3974D0,-0.7820D0,
     :    -7.8468D0,0.5320D0,36.0910D0, 2.4469D0,-30.4606D0,-2.0652D0,
     :     0.0216D0,0.0000D0,-0.0988D0, 0.0000D0,  0.0834D0, 0.0000D0,
     :    -0.3384D0,0.0000D0, 1.5433D0, 0.0000D0, -1.3025D0, 0.0000D0,
     :     0.0179D0,0.0000D0,-0.0813D0, 0.0000D0,  0.0686D0, 0.0000D0 /

      DATA ( ( TIDE(I,J), I=1,6 ), J = 21,40 ) /
     :    -0.0244D0,0.0000D0, 0.1082D0, 0.0000D0, -0.0913D0, 0.0000D0,
     :     0.0470D0,0.0000D0,-0.2004D0, 0.0000D0,  0.1692D0, 0.0000D0,
     :    -0.7341D0,0.0000D0, 3.1240D0, 0.0000D0, -2.6367D0, 0.0000D0,
     :    -0.0526D0,0.0000D0, 0.2235D0, 0.0000D0, -0.1886D0, 0.0000D0,
     :    -0.0508D0,0.0000D0, 0.2073D0, 0.0000D0, -0.1749D0, 0.0000D0,
     :     0.0498D0,0.0000D0,-0.1312D0, 0.0000D0,  0.1107D0, 0.0000D0,
     :     0.1006D0,0.0000D0,-0.2640D0, 0.0000D0,  0.2228D0, 0.0000D0,
     :     0.0395D0,0.0000D0,-0.0968D0, 0.0000D0,  0.0817D0, 0.0000D0,
     :     0.0470D0,0.0000D0,-0.1099D0, 0.0000D0,  0.0927D0, 0.0000D0,
     :     0.1767D0,0.0000D0,-0.4115D0, 0.0000D0,  0.3473D0, 0.0000D0,
     :     0.4352D0,0.0000D0,-1.0093D0, 0.0000D0,  0.8519D0, 0.0000D0,
     :     0.5339D0,0.0000D0,-1.2224D0, 0.0000D0,  1.0317D0, 0.0000D0,
     :    -8.4046D0,0.2500D0,19.1647D0, 0.5701D0,-16.1749D0,-0.4811D0,
     :     0.5443D0,0.0000D0,-1.2360D0, 0.0000D0,  1.0432D0, 0.0000D0,
     :     0.0470D0,0.0000D0,-0.1000D0, 0.0000D0,  0.0844D0, 0.0000D0,
     :    -0.0555D0,0.0000D0, 0.1169D0, 0.0000D0, -0.0987D0, 0.0000D0,
     :     0.1175D0,0.0000D0,-0.2332D0, 0.0000D0,  0.1968D0, 0.0000D0,
     :    -1.8236D0,0.0000D0, 3.6018D0, 0.0000D0, -3.0399D0, 0.0000D0,
     :     0.1316D0,0.0000D0,-0.2587D0, 0.0000D0,  0.2183D0, 0.0000D0,
     :     0.0179D0,0.0000D0,-0.0344D0, 0.0000D0,  0.0290D0, 0.0000D0 /

      DATA ( ( TIDE(I,J), I=1,6 ), J = 41,60 ) /
     :    -0.0855D0,0.0000D0, 0.1542D0, 0.0000D0, -0.1302D0, 0.0000D0,
     :    -0.0573D0,0.0000D0, 0.0395D0, 0.0000D0, -0.0333D0, 0.0000D0,
     :     0.0329D0,0.0000D0,-0.0173D0, 0.0000D0,  0.0146D0, 0.0000D0,
     :    -1.8847D0,0.0000D0, 0.9726D0, 0.0000D0, -0.8209D0, 0.0000D0,
     :     0.2510D0,0.0000D0,-0.0910D0, 0.0000D0,  0.0768D0, 0.0000D0,
     :     1.1703D0,0.0000D0,-0.4135D0, 0.0000D0,  0.3490D0, 0.0000D0,
     :   -49.7174D0,0.4330D0,17.1056D0, 0.1490D0,-14.4370D0,-0.1257D0,
     :    -0.1936D0,0.0000D0, 0.0666D0, 0.0000D0, -0.0562D0, 0.0000D0,
     :     0.0489D0,0.0000D0,-0.0154D0, 0.0000D0,  0.0130D0, 0.0000D0,
     :    -0.5471D0,0.0000D0, 0.1670D0, 0.0000D0, -0.1409D0, 0.0000D0,
     :     0.0367D0,0.0000D0,-0.0108D0, 0.0000D0,  0.0092D0, 0.0000D0,
     :    -0.0451D0,0.0000D0, 0.0082D0, 0.0000D0, -0.0069D0, 0.0000D0,
     :     0.0921D0,0.0000D0,-0.0167D0, 0.0000D0,  0.0141D0, 0.0000D0,
     :     0.8281D0,0.0000D0,-0.1425D0, 0.0000D0,  0.1202D0, 0.0000D0,
     :   -15.8887D0,0.1530D0, 2.7332D0, 0.0263D0, -2.3068D0,-0.0222D0,
     :    -0.1382D0,0.0000D0, 0.0225D0, 0.0000D0, -0.0190D0, 0.0000D0,
     :     0.0348D0,0.0000D0,-0.0053D0, 0.0000D0,  0.0045D0, 0.0000D0,
     :    -0.1372D0,0.0000D0,-0.0079D0, 0.0000D0,  0.0066D0, 0.0000D0,
     :     0.4211D0,0.0000D0,-0.0203D0, 0.0000D0,  0.0171D0, 0.0000D0,
     :    -0.0404D0,0.0000D0, 0.0008D0, 0.0000D0, -0.0007D0, 0.0000D0 /

      DATA ( ( TIDE(I,J), I=1,6 ), J = 61,62 ) /
     :     7.8998D0,0.0000D0, 0.1460D0, 0.0000D0, -0.1232D0, 0.0000D0,
     : -1617.2681D0,0.0000D0,-14.9471D0,0.0000D0, 12.6153D0, 0.0000D0 /

*  -------------------------------------
*   Computation of fundamental arguments
*  -------------------------------------
      CALL FUNDARG(T,L,LP,F,D,OM)

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*  Set initial values to zero.
      DUT    = 0.0D0
      DLOD   = 0.0D0
      DOMEGA = 0.0D0

*  Sum zonal tide terms.
      DO 10 I = 1, NZONT 
*     Formation of multiples of arguments.
         ARG =      MOD ( DBLE ( NFUND(1,I) ) * L 
     .       +            DBLE ( NFUND(2,I) ) * LP 
     .       +            DBLE ( NFUND(3,I) ) * F 
     .       +            DBLE ( NFUND(4,I) ) * D 
     .       +            DBLE ( NFUND(5,I) ) * OM, D2PI )

         IF (ARG.LT.0D0) ARG = ARG + D2PI

*     Evaluate zonal tidal terms.
         DUT    = DUT    + TIDE(1,I) *DSIN(ARG) + TIDE(2,I) *DCOS(ARG)
         DLOD   = DLOD   + TIDE(3,I) *DCOS(ARG) + TIDE(4,I) *DSIN(ARG)
         DOMEGA = DOMEGA + TIDE(5,I) *DCOS(ARG) + TIDE(6,I) *DSIN(ARG)
10    CONTINUE

*  Rescale corrections so that they are in units involving seconds.

      DUT    = DUT    * 1.0D-4
      DLOD   = DLOD   * 1.0D-5
      DOMEGA = DOMEGA * 1.0D-14

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
