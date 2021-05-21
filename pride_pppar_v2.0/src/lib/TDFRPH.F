      SUBROUTINE TDFRPH (jd,sod,IDOOD,FREQ,PHASE)
*+
*  - - - - - - - - - - -
*   T D F R P H 
*  - - - - - - - - - - -
*
*  This routine is part of the International Earth Rotation and
*  Reference Systems Service (IERS) Conventions software collection.
*
*  This subroutine returns the frequency and phase of a tidal
*  constituent when its Doodson number is given as input. 
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
*     idood       i      Doodson number of a tidal constituent
*
*  Returned:
*     freq        d      Frequency of a tidal constituent
*     phase       d      Phase of a tidal constituent (Note 1)
*
*  Notes:
*
*  1) The phases must be decreased by 90 degrees if the sum of the order 
*     and the species number is odd (as for the 2nd degree diurnals, and 
*     3rd degree low frequency and semidiurnals).
*     
*     These phases may need further adjustment to allow for the spherical
*     harmonic normalization used; e.g. for that used for the potential
*     by Cartwright and Tayler, 180 degrees must be added for (species,
*     order) = (1,2), (1,3), or (3,3). 
*
*  Called:
*     TOYMD     Converts year-day of year to year-month-day format 
*     LEAP      Returns true if a given year is a leap year
*     JULDAT    Converts Gregorian date to Julian date 
*     ETUTC     Returns difference of Epheremis Time (ET) and 
*               Coordinated Universal Time (UTC) 
*
*  Test case:
*     given input: For June 25, 2009 (DOY = 176) 0 Hr 0 Min, 0 Sec, M2 tide
*                  INTEGER ITM = / 2009, 176, 0, 0, 0 /
*                  INTEGER IDOOD = / 2, 0, 0, 0, 0, 0 /  
*
*     expected output: FREQ = 1.93227361605688D0
*                      PHASE = 303.343338720297D0
*
*  Test case:
*     given input: For June 25, 2009 (DOY = 176) 12 Hr 1 Min, 45 Sec, M2 tide
*                  INTEGER ITM = / 2009, 176, 12, 1, 45 /
*                  INTEGER IDOOD = / 2, 0, 0, 0, 0, 0 /  
*
*     expected output: FREQ = 1.93227361605689D0
*                      PHASE = 291.997959322689D0
*
*  References:
*
*     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
*     IERS Technical Note No. 36, BKG (2010)
*
*  Revisions:  
*  2009 June   15 B.E.Stetzler  Initial changes to code 
*  2009 August 19 B.E.Stetzler  Capitalized all variables for FORTRAN
*                               77 compatibility
*  2010 March  19 B.E.Stetzler  Provided test case
*  2013 September 11 MSC        Fix error in test case.
*  2015 April  29 M.A. Davis    Corrected error in seconds per day in 
*                               variable DAYFR (from 84600.0D0 to 
*                               86400.0D0)
*  2015 May    21 M.A. Davis    Clarified, corrected existing test case
*                               Added new test case for non-zero hhmmss
*-----------------------------------------------------------------------

      IMPLICIT NONE
      INTEGER I,IDOOD,INITIAL,ITM,ITM2,ITMSAVE,JD
      integer jdsav
      DOUBLE PRECISION YEAR,DELTA,FREQ,PHASE,sod,sodsav,
     .                 D,DAYFR,DD,DJD,F1,F2,F3,F4,F5,
     .                 FD1,FD2,FD3,FD4,FD5,T,timdif,taiutc
      DIMENSION IDOOD(6),ITM2(6),ITMSAVE(5),D(6),DD(6)
      data jdsav,sodsav/0,0.d0/
      SAVE jdsav,sodsav,D,DD

* Common block 'date' stores time information in Universal Time (UT)

c     COMMON/DATE/ITM(5)
c     DATA ITMSAVE/5*0/
*------------------------------------------------------------------------
*  Test to see if time has changed; if so, set the phases and frequencies
*  for each of the Doodson arguments
*------------------------------------------------------------------------
c     INITIAL=0
c     DO I=1,5
c       IF(ITM(I).NE.ITMSAVE(I)) INITIAL=1
c     ENDDO

c     IF(INITIAL.EQ.1) THEN
c       DO I=1,5
c          ITMSAVE(I) = ITM(I)
c       ENDDO
      if(dabs(timdif(jd,sod,jdsav,sodsav)).gt.1.d-3) then
        jdsav=jd
        sodsav=sod

* Convert times to Julian days (UT) then to Julian centuries from J2000.0
*   (ET)

c        CALL TOYMD(ITM,ITM2)
c        JD = JULDAT(ITM2)
c        DAYFR=  ITM(3)/24.0D0 + ITM(4)/1440.0D0 + ITM(5)/86400.0D0
c        YEAR=ITM(1)+(ITM(2)+DAYFR)/(365.0D0+LEAP(ITM(1)))
c        CALL ETUTC(YEAR,DELTA)
c        DJD= JD - 0.5D0 + DAYFR
c        T = (DJD - 2451545.0D0 + DELTA/86400.0D0)/36525.0D0
        t=(jd-51544.5d0+(taiutc(jd)+32.184d0+sod)/86400.d0)/36525.d0


* IERS expressions for the Delaunay arguments, in degrees

        F1 =     134.9634025100D0 +
     .    T*( 477198.8675605000D0 +
     .    T*(      0.0088553333D0 +
     .    T*(      0.0000143431D0 +
     .    T*(     -0.0000000680D0 ))))
        F2 =     357.5291091806D0 +
     .    T*(  35999.0502911389D0 +
     .    T*(     -0.0001536667D0 +
     .    T*(      0.0000000378D0 +
     .    T*(     -0.0000000032D0 ))))
        F3 =      93.2720906200D0 +
     .    T*( 483202.0174577222D0 +
     .    T*(     -0.0035420000D0 +
     .    T*(     -0.0000002881D0 +
     .    T*(      0.0000000012D0 ))))
        F4 =     297.8501954694D0 +
     .    T*( 445267.1114469445D0 +
     .    T*(     -0.0017696111D0 +
     .    T*(      0.0000018314D0 +
     .    T*(     -0.0000000088D0 ))))
        F5 =     125.0445550100D0 +
     .    T*(  -1934.1362619722D0 +
     .    T*(      0.0020756111D0 +
     .    T*(      0.0000021394D0 +
     .    T*(     -0.0000000165D0 ))))

*  Convert to Doodson (Darwin) variables

        D(1) = sod/86400.d0*360.d0 - F4
        D(2) = F3 + F5
        D(3) = D(2) - F4
        D(4) = D(2) - F1
        D(5) = -F5
        D(6) = D(3) - F2

*  Find frequencies of Delauney variables (in cycles/day), and from these
*  the same for the Doodson arguments

        FD1 =  0.0362916471D0 + 0.0000000013D0*T
        FD2 =  0.0027377786D0
        FD3 =  0.0367481951D0 - 0.0000000005D0*T
        FD4 =  0.0338631920D0 - 0.0000000003D0*T
        FD5 = -0.0001470938D0 + 0.0000000003D0*T
        DD(1) = 1.0D0 - FD4
        DD(2) = FD3 + FD5
        DD(3) = DD(2) - FD4
        DD(4) = DD(2) - FD1
        DD(5) = -FD5
        DD(6) = DD(3) - FD2
      ENDIF

*  End of intialization (likely to be called only once)

*  Compute phase and frequency of the given tidal constituent

      FREQ=0.0D0
      PHASE=0.0D0
      DO I = 1,6
         FREQ =  FREQ + IDOOD(I)*DD(I)
         PHASE = PHASE + IDOOD(I)*D(I)
      ENDDO

* Adjust phases so that they fall in the positive range 0 to 360
      PHASE = DMOD(PHASE,360.0D0)
      IF(PHASE.LT.0.0D0) PHASE = PHASE + 360.0D0

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
