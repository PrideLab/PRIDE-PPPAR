!! reference and acknowledgement to RTKLIB
!! modified by Kunlun Zhang
! constants -----------------------------------------------------------------

module rtklib_h_
implicit none

character(5), parameter :: VER_RTKLIB       = "2.4.3"  ! library version
character(3), parameter :: PATCH_LEVEL      = "b33"    ! patch level

real*8, parameter :: PI     = 3.1415926535897932d0     ! pi 
real*8, parameter :: D2R    = (PI/180.d0)              ! deg to rad 
real*8, parameter :: R2D    = (180.d0/PI)              ! rad to deg 
real*8, parameter :: CLIGHT = 299792458.d0             ! speed of light (m/s) 
real*8, parameter :: SC2RAD = 3.1415926535898d0        ! semi-circle to radian (IS-GPS) 
real*8, parameter :: AU     = 149597870691.d0          ! 1 AU (m) 
real*8, parameter :: AS2R   = (D2R/3600.d0)            ! arc sec to radian 

real*8, parameter :: OMGE   = 7.2921151467d-5          ! earth angular velocity (IS-GPS) (rad/s)
real*8, parameter :: RE_WGS84 = 6378137.d0             ! earth semimajor axis (WGS84) (m) 
real*8, parameter :: FE_WGS84 = (1.d0/298.257223563d0) ! earth flattening (WGS84) 

real*8, parameter :: HION     = 350000.d0              ! ionosphere height (m) 
integer*4, parameter :: MAXFREQ  = 7                   ! max NFREQ 

real*8, parameter :: FREQ1     = 1.57542d9             ! L1/E1  frequency (Hz)
real*8, parameter :: FREQ2     = 1.22760d9             ! L2     frequency (Hz)
real*8, parameter :: FREQ5     = 1.17645d9             ! L5/E5a frequency (Hz)
real*8, parameter :: FREQ6     = 1.27875d9             ! E6/LEX frequency (Hz)
real*8, parameter :: FREQ7     = 1.20714d9             ! E5b    frequency (Hz)
real*8, parameter :: FREQ8     = 1.191795d9            ! E5a+b  frequency (Hz)
real*8, parameter :: FREQ9     = 2.492028d9            ! S      frequency (Hz)
real*8, parameter :: FREQ1_GLO = 1.60200d9             ! GLONASS G1 base frequency (Hz)
real*8, parameter :: DFRQ1_GLO = 0.56250d6             ! GLONASS G1 bias frequency (Hz/n)
real*8, parameter :: FREQ2_GLO = 1.24600d9             ! GLONASS G2 base frequency (Hz)
real*8, parameter :: DFRQ2_GLO = 0.43750d6             ! GLONASS G2 bias frequency (Hz/n)
real*8, parameter :: FREQ3_GLO = 1.202025d9            ! GLONASS G3 frequency (Hz)
real*8, parameter :: FREQ1_CMP = 1.561098d9            ! BeiDou B1 frequency (Hz)
real*8, parameter :: FREQ2_CMP = 1.20714d9             ! BeiDou B2 frequency (Hz)
real*8, parameter :: FREQ3_CMP = 1.26852d9             ! BeiDou B3 frequency (Hz)

real*8, parameter :: EFACT_GPS = 1.d0                  ! error factor: GPS
real*8, parameter :: EFACT_GLO = 1.5d0                 ! error factor: GLONASS
real*8, parameter :: EFACT_GAL = 1.d0                  ! error factor: Galileo
real*8, parameter :: EFACT_QZS = 1.d0                  ! error factor: QZSS
real*8, parameter :: EFACT_CMP = 1.d0                  ! error factor: BeiDou
real*8, parameter :: EFACT_IRN = 1.5d0                 ! error factor: IRNSS
real*8, parameter :: EFACT_SBS = 3.d0                  ! error factor: SBAS

integer*4, parameter :: SYS_NONE = 0   !0x00           ! navigation system: none 
integer*4, parameter :: SYS_GPS  = 1   !0x01           ! navigation system: GPS 
integer*4, parameter :: SYS_SBS  = 2   !0x02           ! navigation system: SBAS 
integer*4, parameter :: SYS_GLO  = 4   !0x04           ! navigation system: GLONASS 
integer*4, parameter :: SYS_GAL  = 8   !0x08           ! navigation system: Galileo 
integer*4, parameter :: SYS_QZS  = 16  !0x10           ! navigation system: QZSS 
integer*4, parameter :: SYS_CMP  = 32  !0x20           ! navigation system: BeiDou 
integer*4, parameter :: SYS_IRN  = 64  !0x40           ! navigation system: IRNS 
integer*4, parameter :: SYS_LEO  = 128 !0x80           ! navigation system: LEO 
integer*4, parameter :: SYS_ALL  = 255 !0xFF           ! navigation system: all 

integer*4, parameter :: TSYS_GPS = 0                   ! time system: GPS time 
integer*4, parameter :: TSYS_UTC = 1                   ! time system: UTC 
integer*4, parameter :: TSYS_GLO = 2                   ! time system: GLONASS time 
integer*4, parameter :: TSYS_GAL = 3                   ! time system: Galileo time 
integer*4, parameter :: TSYS_QZS = 4                   ! time system: QZSS time 
integer*4, parameter :: TSYS_CMP = 5                   ! time system: BeiDou time 
integer*4, parameter :: TSYS_IRN = 6                   ! time system: IRNSS time 


!#ifndef NFREQ
integer*4, parameter :: NFREQ    = 3                   ! number of carrier frequencies 
!#endif
integer*4, parameter :: NFREQGLO = 2                   ! number of carrier frequencies of GLONASS 

!#ifndef NEXOBS
integer*4, parameter :: NEXOBS   = 3                   ! number of extended obs codes 
!#endif

integer*4, parameter :: MINPRNGPS = 1                  ! min satellite PRN number of GPS 
integer*4, parameter :: MAXPRNGPS = 32                 ! max satellite PRN number of GPS 
integer*4, parameter :: NSATGPS   = (MAXPRNGPS-MINPRNGPS+1) ! number of GPS satellites 
integer*4, parameter :: NSYSGPS   = 1

integer*4, parameter :: MINPRNGLO = 1                  ! min satellite slot number of GLONASS 
integer*4, parameter :: MAXPRNGLO = 27                 ! max satellite slot number of GLONASS 
integer*4, parameter :: NSATGLO   = (MAXPRNGLO-MINPRNGLO+1) ! number of GLONASS satellites 
integer*4, parameter :: NSYSGLO   = 1

integer*4, parameter :: MINPRNGAL = 1                  ! min satellite PRN number of Galileo 
integer*4, parameter :: MAXPRNGAL = 36                 ! max satellite PRN number of Galileo 
integer*4, parameter :: NSATGAL   = (MAXPRNGAL-MINPRNGAL+1) ! number of Galileo satellites 
integer*4, parameter :: NSYSGAL   = 1

integer*4, parameter :: MINPRNQZS = 193                ! min satellite PRN number of QZSS 
integer*4, parameter :: MAXPRNQZS = 202                ! max satellite PRN number of QZSS 
integer*4, parameter :: MINPRNQZS_S = 183              ! min satellite PRN number of QZSS SAIF 
integer*4, parameter :: MAXPRNQZS_S = 191              ! max satellite PRN number of QZSS SAIF 
integer*4, parameter :: NSATQZS   = (MAXPRNQZS-MINPRNQZS+1) ! number of QZSS satellites 
integer*4, parameter :: NSYSQZS   = 1

integer*4, parameter :: MINPRNCMP = 1                  ! min satellite sat number of BeiDou 
integer*4, parameter :: MAXPRNCMP = 37                 ! max satellite sat number of BeiDou 
integer*4, parameter :: NSATCMP   = (MAXPRNCMP-MINPRNCMP+1) ! number of BeiDou satellites 
integer*4, parameter :: NSYSCMP   = 1

integer*4, parameter :: MINPRNIRN = 1                  ! min satellite sat number of IRNSS 
integer*4, parameter :: MAXPRNIRN = 7                  ! max satellite sat number of IRNSS 
integer*4, parameter :: NSATIRN   = (MAXPRNIRN-MINPRNIRN+1) ! number of IRNSS satellites 
integer*4, parameter :: NSYSIRN   = 1

integer*4, parameter :: MINPRNLEO = 1                  ! min satellite sat number of LEO 
integer*4, parameter :: MAXPRNLEO = 10                 ! max satellite sat number of LEO 
integer*4, parameter :: NSATLEO   = (MAXPRNLEO-MINPRNLEO+1) ! number of LEO satellites 
integer*4, parameter :: NSYSLEO   = 1

integer*4, parameter :: NSYS      = (NSYSGPS+NSYSGLO+NSYSGAL+NSYSQZS+NSYSCMP+NSYSIRN+NSYSLEO) ! number of systems

integer*4, parameter :: MINPRNSBS = 120                ! min satellite PRN number of SBAS 
integer*4, parameter :: MAXPRNSBS = 142                ! max satellite PRN number of SBAS 
integer*4, parameter :: NSATSBS   = (MAXPRNSBS-MINPRNSBS+1) ! number of SBAS satellites 
integer*4, parameter :: MAXSAT    = (NSATGPS+NSATGLO+NSATGAL+NSATQZS+NSATCMP+NSATIRN+NSATSBS+NSATLEO)
                                                       ! max satellite number (1 to MAXSAT) 
integer*4, parameter :: MAXSTA    = 255

integer*4, parameter :: MAXOBS    = 128                ! max number of obs in an epoch
integer*4, parameter :: MAXRCV    = 64                 ! max receiver number (1 to MAXRCV)
integer*4, parameter :: MAXOBSTYPE = 64                ! max number of obs type in RINEX

real*8, parameter :: DTTOL = 0.015d0                   ! tolerance of time difference (s) 

real*8, parameter :: MAXDTOE     = 7200.d0             ! max time difference to GPS Toe (s) 
real*8, parameter :: MAXDTOE_QZS = 7200.d0             ! max time difference to QZSS Toe (s) 
real*8, parameter :: MAXDTOE_GAL = 14400.d0            ! max time difference to Galileo Toe (s) 
real*8, parameter :: MAXDTOE_CMP = 21600.d0            ! max time difference to BeiDou Toe (s) 
real*8, parameter :: MAXDTOE_GLO = 1800.d0             ! max time difference to GLONASS Toe (s) 
real*8, parameter :: MAXDTOE_SBS = 360.d0              ! max time difference to SBAS Toe (s) 
real*8, parameter :: MAXDTOE_S   = 86400.d0            ! max time difference to ephem toe (s) for other 
real*8, parameter :: MAXGDOP     = 300.d0              ! max GDOP 

real*8, parameter :: INT_SWAP_TRAC = 86400.d0          ! swap interval of trace file (s) 
real*8, parameter :: INT_SWAP_STAT = 86400.d0          ! swap interval of solution status file (s) 

integer*4, parameter :: MAXEXFILE  = 1024              ! max number of expanded files 
real*8,    parameter :: MAXSBSAGEF = 30.d0             ! max age of SBAS fast correction (s) 
real*8,    parameter :: MAXSBSAGEL = 1800.d0           ! max age of SBAS long term corr (s) 
integer*4, parameter :: MAXSBSURA  = 8                 ! max URA of SBAS satellite 
integer*4, parameter :: MAXBAND    = 10                ! max SBAS band of IGP 
integer*4, parameter :: MAXNIGP    = 201               ! max number of IGP in SBAS band 
integer*4, parameter :: MAXNGEO    = 4                 ! max number of GEO satellites 
integer*4, parameter :: MAXCOMMENT = 10                ! max number of RINEX comments 
integer*4, parameter :: MAXSTRPATH = 1024              ! max length of stream path 
integer*4, parameter :: MAXSTRMSG  = 1024              ! max length of stream message 
integer*4, parameter :: MAXSTRRTK  = 8                 ! max number of stream in RTK server 
integer*4, parameter :: MAXSBSMSG  = 32                ! max number of SBAS msg in RTK server 
integer*4, parameter :: MAXSOLMSG  = 8191              ! max length of solution message 
integer*4, parameter :: MAXRAWLEN  = 4096              ! max length of receiver raw message 
integer*4, parameter :: MAXERRMSG  = 4096              ! max length of error/warning message 
integer*4, parameter :: MAXANT     = 64                ! max length of station name/antenna type 
integer*4, parameter :: MAXSOLBUF  = 256               ! max number of solution buffer 
integer*4, parameter :: MAXSOLNUM  = 86400             ! max number of solution in program
integer*4, parameter :: MAXOBSBUF  = 128               ! max number of observation data buffer 
integer*4, parameter :: MAXNRPOS   = 16                ! max number of reference positions 
integer*4, parameter :: MAXLEAPS   = 64                ! max number of leap seconds table 
integer*4, parameter :: MAXGISLAYER= 32                ! max number of GIS data layers 
integer*4, parameter :: MAXRCVCMD  = 4096              ! max length of receiver commands 

real*8, parameter :: RNX2VER = 2.10d0                  ! RINEX ver.2 default output version 
real*8, parameter :: RNX3VER = 3.00d0                  ! RINEX ver.3 default output version 

integer*4, parameter :: OBSTYPE_PR  = 1   !#01         ! observation type: pseudorange 
integer*4, parameter :: OBSTYPE_CP  = 2   !#02         ! observation type: carrier-phase 
integer*4, parameter :: OBSTYPE_DOP = 4   !#04         ! observation type: doppler-freq 
integer*4, parameter :: OBSTYPE_SNR = 8   !#08         ! observation type: SNR 
integer*4, parameter :: OBSTYPE_ALL = 255 !#FF         ! observation type: all 

integer*4, parameter :: FREQTYPE_L1 = 1    !0x01       ! frequency type: L1/E1 
integer*4, parameter :: FREQTYPE_L2 = 2    !0x02       ! frequency type: L2/B1 
integer*4, parameter :: FREQTYPE_L5 = 4    !0x04       ! frequency type: L5/E5a/L3 
integer*4, parameter :: FREQTYPE_L6 = 8    !0x08       ! frequency type: E6/LEX/B3 
integer*4, parameter :: FREQTYPE_L7 = 16   !0x10       ! frequency type: E5b/B2 
integer*4, parameter :: FREQTYPE_L8 = 32   !0x20       ! frequency type: E5(a+b) 
integer*4, parameter :: FREQTYPE_L9 = 64   !0x40       ! frequency type: S 
integer*4, parameter :: FREQTYPE_ALL = 255 !0xFF       ! frequency type: all 

integer*4, parameter :: CODE_NONE = 0                  ! obs code: none or unknown 
integer*4, parameter :: CODE_L1C  = 1                  ! obs code: L1C/A,G1C/A,E1C (GPS,GLO,GAL,QZS,SBS) 
integer*4, parameter :: CODE_L1P  = 2                  ! obs code: L1P,G1P    (GPS,GLO) 
integer*4, parameter :: CODE_L1W  = 3                  ! obs code: L1 Z-track (GPS) 
integer*4, parameter :: CODE_L1Y  = 4                  ! obs code: L1Y        (GPS) 
integer*4, parameter :: CODE_L1M  = 5                  ! obs code: L1M        (GPS) 
integer*4, parameter :: CODE_L1N  = 6                  ! obs code: L1codeless (GPS) 
integer*4, parameter :: CODE_L1S  = 7                  ! obs code: L1C(D)     (GPS,QZS) 
integer*4, parameter :: CODE_L1L  = 8                  ! obs code: L1C(P)     (GPS,QZS) 
integer*4, parameter :: CODE_L1E  = 9                  ! (not used) 
integer*4, parameter :: CODE_L1A  = 10                 ! obs code: E1A        (GAL) 
integer*4, parameter :: CODE_L1B  = 11                 ! obs code: E1B        (GAL) 
integer*4, parameter :: CODE_L1X  = 12                 ! obs code: E1B+C,L1C(D+P) (GAL,QZS) 
integer*4, parameter :: CODE_L1Z  = 13                 ! obs code: E1A+B+C,L1SAIF (GAL,QZS) 
integer*4, parameter :: CODE_L2C  = 14                 ! obs code: L2C/A,G1C/A (GPS,GLO) 
integer*4, parameter :: CODE_L2D  = 15                 ! obs code: L2 L1C/A-(P2-P1) (GPS) 
integer*4, parameter :: CODE_L2S  = 16                 ! obs code: L2C(M)     (GPS,QZS) 
integer*4, parameter :: CODE_L2L  = 17                 ! obs code: L2C(L)     (GPS,QZS) 
integer*4, parameter :: CODE_L2X  = 18                 ! obs code: L2C(M+L),B1I+Q (GPS,QZS,CMP) 
integer*4, parameter :: CODE_L2P  = 19                 ! obs code: L2P,G2P    (GPS,GLO) 
integer*4, parameter :: CODE_L2W  = 20                 ! obs code: L2 Z-track (GPS) 
integer*4, parameter :: CODE_L2Y  = 21                 ! obs code: L2Y        (GPS) 
integer*4, parameter :: CODE_L2M  = 22                 ! obs code: L2M        (GPS) 
integer*4, parameter :: CODE_L2N  = 23                 ! obs code: L2codeless (GPS) 
integer*4, parameter :: CODE_L5I  = 24                 ! obs code: L5/E5aI    (GPS,GAL,QZS,SBS) 
integer*4, parameter :: CODE_L5Q  = 25                 ! obs code: L5/E5aQ    (GPS,GAL,QZS,SBS) 
integer*4, parameter :: CODE_L5X  = 26                 ! obs code: L5/E5aI+Q/L5B+C (GPS,GAL,QZS,IRN,SBS) 
integer*4, parameter :: CODE_L7I  = 27                 ! obs code: E5bI,B2I   (GAL,CMP) 
integer*4, parameter :: CODE_L7Q  = 28                 ! obs code: E5bQ,B2Q   (GAL,CMP) 
integer*4, parameter :: CODE_L7X  = 29                 ! obs code: E5bI+Q,B2I+Q (GAL,CMP) 
integer*4, parameter :: CODE_L6A  = 30                 ! obs code: E6A        (GAL) 
integer*4, parameter :: CODE_L6B  = 31                 ! obs code: E6B        (GAL) 
integer*4, parameter :: CODE_L6C  = 32                 ! obs code: E6C        (GAL) 
integer*4, parameter :: CODE_L6X  = 33                 ! obs code: E6B+C,LEXS+L,B3I+Q (GAL,QZS,CMP) 
integer*4, parameter :: CODE_L6Z  = 34                 ! obs code: E6A+B+C    (GAL) 
integer*4, parameter :: CODE_L6S  = 35                 ! obs code: LEXS       (QZS) 
integer*4, parameter :: CODE_L6L  = 36                 ! obs code: LEXL       (QZS) 
integer*4, parameter :: CODE_L8I  = 37                 ! obs code: E5(a+b)I   (GAL) 
integer*4, parameter :: CODE_L8Q  = 38                 ! obs code: E5(a+b)Q   (GAL) 
integer*4, parameter :: CODE_L8X  = 39                 ! obs code: E5(a+b)I+Q (GAL) 
integer*4, parameter :: CODE_L2I  = 40                 ! obs code: B1I        (BDS) 
integer*4, parameter :: CODE_L2Q  = 41                 ! obs code: B1Q        (BDS) 
integer*4, parameter :: CODE_L6I  = 42                 ! obs code: B3I        (BDS) 
integer*4, parameter :: CODE_L6Q  = 43                 ! obs code: B3Q        (BDS) 
integer*4, parameter :: CODE_L3I  = 44                 ! obs code: G3I        (GLO) 
integer*4, parameter :: CODE_L3Q  = 45                 ! obs code: G3Q        (GLO) 
integer*4, parameter :: CODE_L3X  = 46                 ! obs code: G3I+Q      (GLO) 
integer*4, parameter :: CODE_L1I  = 47                 ! obs code: B1I        (BDS) 
integer*4, parameter :: CODE_L1Q  = 48                 ! obs code: B1Q        (BDS) 
integer*4, parameter :: CODE_L5A  = 49                 ! obs code: L5A SPS    (IRN) 
integer*4, parameter :: CODE_L5B  = 50                 ! obs code: L5B RS(D)  (IRN) 
integer*4, parameter :: CODE_L5C  = 51                 ! obs code: L5C RS(P)  (IRN) 
integer*4, parameter :: CODE_L9A  = 52                 ! obs code: SA SPS     (IRN) 
integer*4, parameter :: CODE_L9B  = 53                 ! obs code: SB RS(D)   (IRN) 
integer*4, parameter :: CODE_L9C  = 54                 ! obs code: SC RS(P)   (IRN) 
integer*4, parameter :: CODE_L9X  = 55                 ! obs code: SB+C       (IRN) 
integer*4, parameter :: MAXCODE   = 55                 ! max number of obs code 

integer*4, parameter :: PMODE_SINGLE = 0               ! positioning mode: single 
integer*4, parameter :: PMODE_DGPS   = 1               ! positioning mode: DGPS/DGNSS 
integer*4, parameter :: PMODE_KINEMA = 2               ! positioning mode: kinematic 
integer*4, parameter :: PMODE_STATIC = 3               ! positioning mode: static 
integer*4, parameter :: PMODE_MOVEB  = 4               ! positioning mode: moving-base 
integer*4, parameter :: PMODE_FIXED  = 5               ! positioning mode: fixed 
integer*4, parameter :: PMODE_PPP_KINEMA = 6           ! positioning mode: PPP-kinemaric 
integer*4, parameter :: PMODE_PPP_STATIC = 7           ! positioning mode: PPP-static 
integer*4, parameter :: PMODE_PPP_FIXED = 8            ! positioning mode: PPP-fixed 

integer*4, parameter :: SOLF_LLH  = 0                  ! solution format: lat/lon/height 
integer*4, parameter :: SOLF_XYZ  = 1                  ! solution format: x/y/z-ecef 
integer*4, parameter :: SOLF_ENU  = 2                  ! solution format: e/n/u-baseline 
integer*4, parameter :: SOLF_NMEA = 3                  ! solution format: NMEA-183 
integer*4, parameter :: SOLF_STAT = 4                  ! solution format: solution status 
integer*4, parameter :: SOLF_GSIF = 5                  ! solution format: GSI F1/F2 

integer*4, parameter :: SOLQ_NONE  = 0                 ! solution status: no solution 
integer*4, parameter :: SOLQ_FIX   = 1                 ! solution status: fix 
integer*4, parameter :: SOLQ_FLOAT = 2                 ! solution status: float 
integer*4, parameter :: SOLQ_SBAS  = 3                 ! solution status: SBAS 
integer*4, parameter :: SOLQ_DGPS  = 4                 ! solution status: DGPS/DGNSS 
integer*4, parameter :: SOLQ_SINGLE = 5                ! solution status: single 
integer*4, parameter :: SOLQ_PPP   = 6                 ! solution status: PPP 
integer*4, parameter :: SOLQ_DR    = 7                 ! solution status: dead reconing 
integer*4, parameter :: MAXSOLQ    = 7                 ! max number of solution status 

integer*4, parameter :: TIMES_GPST = 0                 ! time system: gps time 
integer*4, parameter :: TIMES_UTC  = 1                 ! time system: utc 
integer*4, parameter :: TIMES_JST  = 2                 ! time system: jst 

integer*4, parameter :: IONOOPT_OFF  = 0               ! ionosphere option: correction off 
integer*4, parameter :: IONOOPT_BRDC = 1               ! ionosphere option: broadcast model 
integer*4, parameter :: IONOOPT_SBAS = 2               ! ionosphere option: SBAS model 
integer*4, parameter :: IONOOPT_IFLC = 3               ! ionosphere option: L1/L2 or L1/L5 iono-free LC 
integer*4, parameter :: IONOOPT_EST  = 4               ! ionosphere option: estimation 
integer*4, parameter :: IONOOPT_TEC  = 5               ! ionosphere option: IONEX TEC model 
integer*4, parameter :: IONOOPT_QZS  = 6               ! ionosphere option: QZSS broadcast model 
integer*4, parameter :: IONOOPT_LEX  = 7               ! ionosphere option: QZSS LEX ionospehre 
integer*4, parameter :: IONOOPT_STEC = 8               ! ionosphere option: SLANT TEC model 

integer*4, parameter :: TROPOPT_OFF  = 0               ! troposphere option: correction off 
integer*4, parameter :: TROPOPT_SAAS = 1               ! troposphere option: Saastamoinen model 
integer*4, parameter :: TROPOPT_SBAS = 2               ! troposphere option: SBAS model 
integer*4, parameter :: TROPOPT_EST  = 3               ! troposphere option: ZTD estimation 
integer*4, parameter :: TROPOPT_ESTG = 4               ! troposphere option: ZTD+grad estimation 
integer*4, parameter :: TROPOPT_ZTD  = 5               ! troposphere option: ZTD correction 

integer*4, parameter :: EPHOPT_BRDC = 0                ! ephemeris option: broadcast ephemeris 
integer*4, parameter :: EPHOPT_PREC = 1                ! ephemeris option: precise ephemeris 
integer*4, parameter :: EPHOPT_SBAS = 2                ! ephemeris option: broadcast + SBAS 
integer*4, parameter :: EPHOPT_SSRAPC = 3              ! ephemeris option: broadcast + SSR_APC 
integer*4, parameter :: EPHOPT_SSRCOM = 4              ! ephemeris option: broadcast + SSR_COM 
integer*4, parameter :: EPHOPT_LEX  = 5                ! ephemeris option: QZSS LEX ephemeris 

integer*4, parameter :: ARMODE_OFF  = 0                ! AR mode: off 
integer*4, parameter :: ARMODE_CONT = 1                ! AR mode: continuous 
integer*4, parameter :: ARMODE_INST = 2                ! AR mode: instantaneous 
integer*4, parameter :: ARMODE_FIXHOLD = 3             ! AR mode: fix and hold 
integer*4, parameter :: ARMODE_WLNL = 4                ! AR mode: wide lane/narrow lane 
integer*4, parameter :: ARMODE_TCAR = 5                ! AR mode: triple carrier ar 

integer*4, parameter :: SBSOPT_LCORR = 1               ! SBAS option: long term correction 
integer*4, parameter :: SBSOPT_FCORR = 2               ! SBAS option: fast correction 
integer*4, parameter :: SBSOPT_ICORR = 4               ! SBAS option: ionosphere correction 
integer*4, parameter :: SBSOPT_RANGE = 8               ! SBAS option: ranging 

integer*4, parameter :: POSOPT_POS   = 0               ! pos option: LLH/XYZ 
integer*4, parameter :: POSOPT_SINGLE = 1              ! pos option: average of single pos 
integer*4, parameter :: POSOPT_FILE  = 2               ! pos option: read from pos file 
integer*4, parameter :: POSOPT_RINEX = 3               ! pos option: rinex header pos 
integer*4, parameter :: POSOPT_RTCM  = 4               ! pos option: rtcm station pos 
integer*4, parameter :: POSOPT_RAW   = 5               ! pos option: raw station pos 

integer*4, parameter :: STR_NONE     = 0               ! stream type: none 
integer*4, parameter :: STR_SERIAL   = 1               ! stream type: serial 
integer*4, parameter :: STR_FILE     = 2               ! stream type: file 
integer*4, parameter :: STR_TCPSVR   = 3               ! stream type: TCP server 
integer*4, parameter :: STR_TCPCLI   = 4               ! stream type: TCP client 
integer*4, parameter :: STR_NTRIPSVR = 6               ! stream type: NTRIP server 
integer*4, parameter :: STR_NTRIPCLI = 7               ! stream type: NTRIP client 
integer*4, parameter :: STR_FTP      = 8               ! stream type: ftp 
integer*4, parameter :: STR_HTTP     = 9               ! stream type: http 
integer*4, parameter :: STR_NTRIPC_S = 10              ! stream type: NTRIP caster server 
integer*4, parameter :: STR_NTRIPC_C = 11              ! stream type: NTRIP caster client 
integer*4, parameter :: STR_UDPSVR   = 12              ! stream type: UDP server 
integer*4, parameter :: STR_UDPCLI   = 13              ! stream type: UDP server 
integer*4, parameter :: STR_MEMBUF   = 14              ! stream type: memory buffer 

integer*4, parameter :: STRFMT_RTCM2 = 0               ! stream format: RTCM 2 
integer*4, parameter :: STRFMT_RTCM3 = 1               ! stream format: RTCM 3 
integer*4, parameter :: STRFMT_OEM4  = 2               ! stream format: NovAtel OEMV/4 
integer*4, parameter :: STRFMT_OEM3  = 3               ! stream format: NovAtel OEM3 
integer*4, parameter :: STRFMT_UBX   = 4               ! stream format: u-blox LEA-*T 
integer*4, parameter :: STRFMT_SS2   = 5               ! stream format: NovAtel Superstar II 
integer*4, parameter :: STRFMT_CRES  = 6               ! stream format: Hemisphere 
integer*4, parameter :: STRFMT_STQ   = 7               ! stream format: SkyTraq S1315F 
integer*4, parameter :: STRFMT_GW10  = 8               ! stream format: Furuno GW10 
integer*4, parameter :: STRFMT_JAVAD = 9               ! stream format: JAVAD GRIL/GREIS 
integer*4, parameter :: STRFMT_NVS   = 10              ! stream format: NVS NVC08C 
integer*4, parameter :: STRFMT_BINEX = 11              ! stream format: BINEX 
integer*4, parameter :: STRFMT_RT17  = 12              ! stream format: Trimble RT17 
integer*4, parameter :: STRFMT_SEPT  = 13              ! stream format: Septentrio 
integer*4, parameter :: STRFMT_CMR   = 14              ! stream format: CMR/CMR+ 
integer*4, parameter :: STRFMT_TERSUS = 15             ! stream format: TERSUS 
integer*4, parameter :: STRFMT_LEXR  = 16              ! stream format: Furuno LPY-10000 
integer*4, parameter :: STRFMT_RINEX = 17              ! stream format: RINEX 
integer*4, parameter :: STRFMT_SP3   = 18              ! stream format: SP3 
integer*4, parameter :: STRFMT_RNXCLK = 19             ! stream format: RINEX CLK 
integer*4, parameter :: STRFMT_SBAS  = 20              ! stream format: SBAS messages 
integer*4, parameter :: STRFMT_NMEA  = 21              ! stream format: NMEA 0183 
integer*4, parameter :: MAXRCVFMT    = 15              ! max number of receiver format 

integer*4, parameter :: STR_MODE_R   = 1               ! stream mode: read 
integer*4, parameter :: STR_MODE_W   = 2               ! stream mode: write 
integer*4, parameter :: STR_MODE_RW  = 3               ! stream mode: read/write 

integer*4, parameter :: GEOID_EMBEDDED    = 0          ! geoid model: embedded geoid 
integer*4, parameter :: GEOID_EGM96_M150  = 1          ! geoid model: EGM96 15x15' 
integer*4, parameter :: GEOID_EGM2008_M25 = 2          ! geoid model: EGM2008 2.5x2.5' 
integer*4, parameter :: GEOID_EGM2008_M10 = 3          ! geoid model: EGM2008 1.d0x1.d0' 
integer*4, parameter :: GEOID_GSI2000_M15 = 4          ! geoid model: GSI geoid 2000 1.d0x1.5' 
integer*4, parameter :: GEOID_RAF09       = 5          ! geoid model: IGN RAF09 for France 1.5'x2' 

character(1), parameter :: COMMENTH = '%'              ! comment line indicator for solution 
character(12), parameter :: MSG_DISCONN = '$_DISCONNECT' ! disconnect message

integer*4, parameter :: DLOPT_FORCE   = 1              ! download option: force download existing
integer*4, parameter :: DLOPT_KEEPCMP = 2              ! download option: keep compressed file
integer*4, parameter :: DLOPT_HOLDERR = 4              ! download option: hold on error file
integer*4, parameter :: DLOPT_HOLDLST = 8              ! download option: hold on listing file

integer*4, parameter :: LLI_SLIP   = 1                 ! LLI: cycle-slip
integer*4, parameter :: LLI_HALFC  = 2                 ! LLI: half-cycle not resovled
integer*4, parameter :: LLI_BOCTRK = 4                 ! LLI: boc tracking of mboc signal
integer*4, parameter :: LLI_HALFA  = 64                ! LLI: half-cycle added
integer*4, parameter :: LLI_HALFS  = 128               ! LLI: half-cycle subtracted

integer*4, parameter :: IMUFMT_KVH = 1                 ! imu data format KVH

real*8, parameter :: P2_5  = 0.03125d0                 ! 2^-5
real*8, parameter :: P2_6  = 0.015625d0                ! 2^-6
real*8, parameter :: P2_11 = 4.882812500000000d-04     ! 2^-11
real*8, parameter :: P2_15 = 3.051757812500000d-05     ! 2^-15
real*8, parameter :: P2_17 = 7.629394531250000d-06     ! 2^-17
real*8, parameter :: P2_19 = 1.907348632812500d-06     ! 2^-19
real*8, parameter :: P2_20 = 9.536743164062500d-07     ! 2^-20
real*8, parameter :: P2_21 = 4.768371582031250d-07     ! 2^-21
real*8, parameter :: P2_23 = 1.192092895507810d-07     ! 2^-23
real*8, parameter :: P2_24 = 5.960464477539063d-08     ! 2^-24
real*8, parameter :: P2_27 = 7.450580596923828d-09     ! 2^-27
real*8, parameter :: P2_29 = 1.862645149230957d-09     ! 2^-29
real*8, parameter :: P2_30 = 9.313225746154785d-10     ! 2^-30
real*8, parameter :: P2_31 = 4.656612873077393d-10     ! 2^-31
real*8, parameter :: P2_32 = 2.328306436538696d-10     ! 2^-32
real*8, parameter :: P2_33 = 1.164153218269348d-10     ! 2^-33
real*8, parameter :: P2_35 = 2.910383045673370d-11     ! 2^-35
real*8, parameter :: P2_38 = 3.637978807091710d-12     ! 2^-38
real*8, parameter :: P2_39 = 1.818989403545856d-12     ! 2^-39
real*8, parameter :: P2_40 = 9.094947017729280d-13     ! 2^-40
real*8, parameter :: P2_43 = 1.136868377216160d-13     ! 2^-43
real*8, parameter :: P2_48 = 3.552713678800501d-15     ! 2^-48
real*8, parameter :: P2_50 = 8.881784197001252d-16     ! 2^-50
real*8, parameter :: P2_55 = 2.775557561562891d-17     ! 2^-55

integer*4, parameter :: FPOUT  = 10                    ! solution files  (except for 1/2/5/6)
integer*4, parameter :: FPREAD = 11                    ! obs/nav  files
integer*4, parameter :: FPOPT  = 12                    ! option   files
integer*4, parameter :: FPADD  = 13                    ! additional files
integer*4, parameter :: FPCLK  = 14                    ! additional files (sat clock bias)
integer*4, parameter :: FPPOS  = 15                    ! additional files (sat pos)

! type definitions ----------------------------------------------------------
type gtime_t                 ! time struct 
    !time_t time             ! time (s) expressed by standard time_t
    integer*4 time
    real*8 sec               ! fraction of second under 1 s 
end type

type obsd_t                  ! observation data record 
    type(gtime_t) time       ! receiver sampling time (GPST) 
    integer*4 sat,rcv        ! satellite/receiver number 
    integer*4 SNR  (NFREQ+NEXOBS) ! signal strength (0.25 dBHz) 
    integer*4 LLI  (NFREQ+NEXOBS) ! loss of lock indicator 
    integer*4 code (NFREQ+NEXOBS) ! code indicator (CODE_???) 
    integer*4 qualL(NFREQ+NEXOBS) ! quality of carrier phase measurement 
    integer*4 qualP(NFREQ+NEXOBS) ! quality of pseudorange measurement 
    real*8 L(NFREQ+NEXOBS)   ! observation data carrier-phase (cycle) 
    real*8 P(NFREQ+NEXOBS)   ! observation data pseudorange (m) 
    real*8 D(NFREQ+NEXOBS)   ! observation data doppler frequency (Hz) 
end type

type obs_t                   ! observation data 
    integer*4 n,nmax         ! number of obervation data/allocated 
    real*8 tint              ! time interval
    type(obsd_t), pointer :: mydata(:)  ! observation data records 
end type

type alm_t                   ! almanac type 
    integer*4 sat            ! satellite number 
    integer*4 svh            ! sv health (0:ok) 
    integer*4 svconf         ! as and sv config 
    integer*4 week           ! GPS/QZS: gps week, GAL: galileo week 
    type(gtime_t) toa        ! Toa 
                             ! SV orbit parameters 
    real*8 A,e,i0,OMG0,omg,M0,OMGd
    real*8 toas              ! Toa (s) in week 
    real*8 f0,f1             ! SV clock parameters (af0,af1) 
end type

type eph_t                   ! GPS/QZS/GAL broadcast ephemeris type 
    integer*4 sat            ! satellite number 
    integer*4 iode,iodc      ! IODE,IODC 
    integer*4 sva            ! SV accuracy (URA index) 
    integer*4 svh            ! SV health (0:ok) 
    integer*4 week           ! GPS/QZS: gps week, GAL: galileo week 
    integer*4 code           ! GPS/QZS: code on L2, GAL/CMP: data sources 
    integer*4 flag           ! GPS/QZS: L2 P data flag, CMP: nav type 
    type(gtime_t) toe,toc,ttr ! Toe,Toc,T_trans 
                             ! SV orbit parameters 
    real*8 A,e,i0,OMG0,omg,M0,deln,OMGd,idot
    real*8 crc,crs,cuc,cus,cic,cis
    real*8 toes              ! Toe (s) in week 
    real*8 fit               ! fit interval (h) 
    real*8 f0,f1,f2          ! SV clock parameters (af0,af1,af2) 
    real*8 tgd(4)            ! group delay parameters 
                             ! GPS/QZS:tgd(0)=TGD 
                             ! GAL    :tgd(0)=BGD E5a/E1,tgd(1)=BGD E5b/E1 
                             ! CMP    :tgd(0)=BGD1,tgd(1)=BGD2 
    real*8 Adot,ndot         ! Adot,ndot for CNAV 
end type

type geph_t                  ! GLONASS broadcast ephemeris type 
    integer*4 sat            ! satellite number 
    integer*4 iode           ! IODE (0-6 bit of tb field) 
    integer*4 frq            ! satellite frequency number 
    integer*4 svh,sva,age    ! satellite health, accuracy, age of operation 
    type(gtime_t) toe        ! epoch of epherides (gpst) 
    type(gtime_t) tof        ! message frame time (gpst) 
    real*8 pos(3)            ! satellite position (ecef) (m) 
    real*8 vel(3)            ! satellite velocity (ecef) (m/s) 
    real*8 acc(3)            ! satellite acceleration (ecef) (m/s^2) 
    real*8 taun,gamn         ! SV clock bias (s)/relative freq bias 
    real*8 dtaun             ! delay between L1 and L2 (s) 
end type

type trop_t                  ! trop data type 
    type(gtime_t) time       ! time (GPST) 
    real*8 trp(3)            ! zenith tropos delay/gradient (m) 
    real*8 std(3)            ! std-dev (m) 
end type

type nav_t                   ! navigation data type 
    integer*4 n,nmax         ! number of broadcast ephemeris 
    integer*4 ng,ngmax       ! number of glonass ephemeris 
    integer*4 ns,nsmax       ! number of sbas ephemeris 
    integer*4 ne,nemax       ! number of precise ephemeris 
    integer*4 nc,ncmax       ! number of precise clock 
    integer*4 na,namax       ! number of almanac data 
    integer*4 nt,ntmax       ! number of tec grid data 
    integer*4 nf,nfmax       ! number of satellite fcb data 
    type(eph_t), pointer :: eph(:)   ! GPS/QZS/GAL ephemeris 
    type(geph_t), pointer :: geph(:) ! GLONASS ephemeris 
    type(alm_t), pointer :: alm(:)   ! almanac data 
    real*8 utc_gps(4)        ! GPS delta-UTC parameters {A0,A1,T,W} 
    real*8 utc_glo(4)        ! GLONASS UTC GPS time parameters 
    real*8 utc_gal(4)        ! Galileo UTC GPS time parameters 
    real*8 utc_qzs(4)        ! QZS UTC GPS time parameters 
    real*8 utc_cmp(4)        ! BeiDou UTC parameters 
    real*8 utc_irn(4)        ! IRNSS UTC parameters 
    real*8 utc_sbs(4)        ! SBAS UTC parameters 
    real*8 ion_gps(8)        ! GPS iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} 
    real*8 ion_gal(4)        ! Galileo iono model parameters {ai0,ai1,ai2,0} 
    real*8 ion_qzs(8)        ! QZSS iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} 
    real*8 ion_cmp(8)        ! BeiDou iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} 
    real*8 ion_irn(8)        ! IRNSS iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} 
    integer*4 leaps          ! leap seconds (s) 
    real*8 lam(MAXSAT,NFREQ) ! carrier wave lengths (m) 
    real*8 cbias(MAXSAT,3)   ! satellite dcb (0:p1-p2,1:p1-c1,2:p2-c2) (m) 
    real*8 rbias(MAXRCV,2,3) ! receiver dcb (0:p1-p2,1:p1-c1,2:p2-c2) (m) 
    real*8 wlbias(MAXSAT)    ! wide-lane bias (cycle) 
    real*8 glo_cpbias(4)     ! glonass code-phase bias {1C,1P,2C,2P} (m) 
    character(MAXPRNGLO+1) glo_fcn  ! glonass frequency channel number + 8 
end type

type sta_t                   ! station parameter type 
    character(MAXANT) name   ! marker name 
    character(MAXANT) marker ! marker number 
    character(MAXANT) antdes ! antenna descriptor 
    character(MAXANT) antsno ! antenna serial number 
    character(MAXANT) rectype ! receiver type descriptor 
    character(MAXANT) recver ! receiver firmware version 
    character(MAXANT) recsno ! receiver serial number 
    integer*4 antsetup       ! antenna setup id 
    integer*4 itrf           ! ITRF realization year 
    integer*4 deltype        ! antenna delta type (0:enu,1:xyz) 
    real*8 pos(3)            ! station position (ecef) (m) 
    real*8 del(3)            ! antenna position delta (e/n/u or x/y/z) (m) 
    real*8 hgt               ! antenna height (m) 
end type

type anyvar
    integer*4 mytype         ! 1-int; 2-double; 3-char*; 0-null
    integer*4 ivalue
    real*8    rvalue
    character(1024) cvalue
end type

type opt_t                   ! option type
    character(1024) name     ! option name 
    integer*4 myformat       ! option format (0:integer*4,1:real*8,2:string,3:enum) 
    type(anyvar) var
    !void *var               ! pointer to option variable 
    character(1024) comment  ! option comment/enum labels/unit 
end type

type snrmask_t               ! SNR mask type 
    integer*4 ena(2)         ! enable flag {rover,base} 
    real*8 mask(NFREQ,9)     ! mask (dBHz) at 5,10,...85 deg 
end type

type prcopt_t                ! processing options type 
    integer*4 mode           ! positioning mode (PMODE_???) 
    integer*4 soltype        ! solution type (0:forward,1:backward,2:combined) 
    integer*4 nf             ! number of frequencies (1:L1,2:L1+L2,3:L1+L2+L5) 
    integer*4 navsys         ! navigation system 
    real*8 elmin             ! elevation mask angle (rad) 
    integer*4 sateph         ! satellite ephemeris/clock (EPHOPT_???) 
    integer*4 ionoopt        ! ionosphere option (IONOOPT_???) 
    integer*4 tropopt        ! troposphere option (TROPOPT_???) 
    integer*4 niter          ! number of filter iteration 
    real*8 maxinno           ! reject threshold of innovation (m) 
    real*8 maxgdop           ! reject threshold of gdop 
    integer(4) exsats(MAXSAT) ! excluded satellites (1:excluded,2:included)
    character(256) rnxopt(2) ! rinex options {rover,base}
    integer*4 posopt(6)      ! positioning options
end type

type sol_t                   ! solution type 
    type(gtime_t) time       ! time (GPST) 
    type(gtime_t) time0      ! time (GPST) of observation
    real*8 rr(6)             ! position/velocity (m|m/s) 
                             ! {x,y,z,vx,vy,vz} or {e,n,u,ve,vn,vu} 
    real*8 qr(6)             ! position variance/covariance (m^2) 
                             ! {c_xx,c_yy,c_zz,c_xy,c_yz,c_zx} or 
                             ! {c_ee,c_nn,c_uu,c_en,c_nu,c_ue} 
    real*8 dtr(6)            ! receiver clock bias to time systems (s) 
    integer*4 mytype         ! type (0:xyz-ecef,1:enu-baseline)
    integer*4 stat           ! solution status (SOLQ_???) 
    integer*4 ns             ! number of valid satellites 
    real*8 dop(4)            ! dilution of precision (gdop, pdop, hdop, vdop)
end type

type solopt_t                ! solution options type 
    integer*4 posf           ! solution format (SOLF_???) 
    integer*4 times          ! time system (TIMES_???) 
    integer*4 timef          ! time format (0:sssss.s,1:yyyy/mm/dd hh:mm:ss.s) 
    integer*4 timeu          ! time digits under decimal point 
    integer*4 degf           ! latitude/longitude format (0:ddd.ddd,1:ddd mm ss) 
    integer*4 outhead        ! output header (0:no,1:yes) 
    integer*4 outopt         ! output processing options (0:no,1:yes) 
    integer*4 outvel         ! output velocity options (0:no,1:yes) 
    integer*4 datum          ! datum (0:WGS84,1:Tokyo) 
    integer*4 height         ! height (0:ellipsoidal,1:geodetic) 
    integer*4 geoid          ! geoid model (0:EGM96,1:JGD2000) 
    integer*4 solstatic      ! solution of static mode (0:all,1:single) 
    integer*4 sstat          ! solution statistics level (0:off,1:states,2:residuals) 
    integer*4 trace          ! debug trace level (0:off,1-5:debug) 
    integer*4 issingle       ! single solution (0:no,1:yes)
    real*8 nmeaintv(2)       ! nmea output interval (s) ( .lt. 0:no,0:all) 
                             ! nmeaintv(0):gprmc,gpgga,nmeaintv(1):gpgsv 
    character(64) sep        ! field separator 
    character(64) prog       ! program name 
    real*8 maxsolstd         ! max std-dev for solution output (m) (0:all) 
end type

type rnxopt_t                ! RINEX options type 
    type(gtime_t) ts,te      ! time start/end 
    real*8 tint              ! time interval (s) 
    real*8 ttol              ! time tolerance (s) 
    real*8 tunit             ! time unit for multiple-session (s) 
    real*8 rnxver            ! RINEX version 
    integer*4 navsys         ! navigation system 
    integer*4 obstype        ! observation type 
    integer*4 freqtype       ! frequency type 
    character(64) mask(7)    ! code mask {GPS,GLO,GAL,QZS,SBS,CMP,IRN} 
    character(32) staid      ! station id for rinex file name 
    character(32) prog       ! program 
    character(32) runby      ! run-by 
    character(64) marker     ! marker name 
    character(32) markerno   ! marker number 
    character(32) markertype ! marker type (ver.3) 
    character(32) name(2)    ! observer/agency 
    character(32) rec (3)    ! receiver #/type/vers 
    character(32) ant (3)    ! antenna #/type 
    integer*4 exsats(MAXSAT) ! excluded satellites 
    type(gtime_t) tstart     ! first obs time 
    type(gtime_t) tend       ! last obs time 
    character(4) tobs(7, MAXOBSTYPE) ! obs types {GPS,GLO,GAL,QZS,SBS,CMP,IRN} 
    integer*4 nobs(7)        ! number of obs types {GPS,GLO,GAL,QZS,SBS,CMP,IRN} 
end type

type ssat_t                  ! satellite status type 
    integer*4 sys            ! navigation system 
    integer*4 vs             ! valid satellite flag single 
    real*8 azel(2)           ! azimuth/elevation angles {az,el} (rad) 
    real*8 resp(NFREQ)       ! residuals of pseudorange (m) 
    real*8 resc(NFREQ)       ! residuals of carrier-phase (m) 
    real*8 icbias(NFREQ)     ! glonass IC bias (cycles) 
    integer*4 vsat(NFREQ)    ! valid satellite flag 
    integer*4 snr (NFREQ)    ! signal strength (0.25 dBHz) 
    integer*4 lock (NFREQ)   ! lock counter of phase 
end type

type rtk_t                   ! RTK control/result type 
    type(sol_t) sol          ! RTK solution 
    real*8 rb(6)             ! base position/velocity (ecef) (m|m/s) 
    integer*4 nx,na          ! number of float states/fixed states 
    real*8 tt                ! time difference between current and previous (s) 
    real*8, pointer :: x(:), P(:,:)   ! float states and their covariance 
    real*8, pointer :: xa(:),Pa(:,:)  ! fixed states and their covariance 
    integer*4 nfix           ! number of continuous fixes of ambiguity 
    integer*4 excsat         ! index of next satellite to be excluded for partial ambiguity resolution 
    real*8 com_bias          ! phase bias common between all sats (used to be distributed to all sats 
    type(ssat_t) ssat(MAXSAT) ! satellite status 
    integer*4 neb            ! bytes in error message buffer 
    type(prcopt_t) opt       ! processing options 
    integer*4 initial_mode   ! initial positioning mode 
end type

! global variables ----------------------------------------------------------
real*8, save :: chisqr(100)               ! chi-sqr(n) table (alpha=0.001) 
real*8, save :: ura_value(15)             ! ura max values
type(prcopt_t), save :: prcopt_default    ! defaults processing options
type(solopt_t), save :: solopt_default    ! defaults solution output options
type(opt_t), save :: sysopts(38)          ! system options table 
! type(sol_t), save :: allsol_(MAXSOLNUM) ! all solutions
type(sol_t), pointer, save :: allsol_(:)  ! all solutions
integer*4, save :: solindex_              ! solution index
integer*4, save :: headwritten_           ! output the header (0-none, 1-written)
character(1024), save :: clkfile_

end module rtklib_h_
