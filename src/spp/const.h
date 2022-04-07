
! Definition of constant ----------------------------------------------------
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
integer*4, parameter :: MAXFREQ = 7                    ! max NFREQ 

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

integer*4, parameter :: NFREQ    = 3                   ! number of carrier frequencies 
integer*4, parameter :: NFREQGLO = 2                   ! number of carrier frequencies of GLONASS 
integer*4, parameter :: NEXOBS   = 3                   ! number of extended obs codes 

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

real*8, parameter :: DTTOL       = 0.015d0             ! tolerance of time difference (s) 
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

integer*4, parameter :: OBSTYPE_PR  = 1    !#01        ! observation type: pseudorange 
integer*4, parameter :: OBSTYPE_CP  = 2    !#02        ! observation type: carrier-phase 
integer*4, parameter :: OBSTYPE_DOP = 4    !#04        ! observation type: doppler-freq 
integer*4, parameter :: OBSTYPE_SNR = 8    !#08        ! observation type: SNR 
integer*4, parameter :: OBSTYPE_ALL = 255  !#FF        ! observation type: all 

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

! rtkcmn.f90 ----------------------------------------------------------------
real*8, parameter :: MAX_VAR_EPH = 9d4                 ! max variance eph to reject satellite (m^2) 
real*8, parameter :: gpst0(6)=(/1980,1, 6,0,0,0/)      ! gps time reference 
real*8, parameter :: gst0 (6)=(/1999,8,22,0,0,0/)      ! galileo system time reference 
real*8, parameter :: bdt0 (6)=(/2006,1, 1,0,0,0/)      ! beidou time reference 

real*8, parameter :: lam_carr(MAXFREQ)=(/&             ! carrier wave length (m)
    CLIGHT/FREQ1,CLIGHT/FREQ2,CLIGHT/FREQ5,CLIGHT/FREQ6,CLIGHT/FREQ7,&
    CLIGHT/FREQ8,CLIGHT/FREQ9/)

! leap seconds (y,m,d,h,m,s,utc-gpst) 
real*8, parameter :: leaps(19,7)=reshape((/&
    2017,1,1,0,0,0,-18,&
    2015,7,1,0,0,0,-17,&
    2012,7,1,0,0,0,-16,&
    2009,1,1,0,0,0,-15,&
    2006,1,1,0,0,0,-14,&
    1999,1,1,0,0,0,-13,&
    1997,7,1,0,0,0,-12,&
    1996,1,1,0,0,0,-11,&
    1994,7,1,0,0,0,-10,&
    1993,7,1,0,0,0, -9,&
    1992,7,1,0,0,0, -8,&
    1991,1,1,0,0,0, -7,&
    1990,1,1,0,0,0, -6,&
    1988,1,1,0,0,0, -5,&
    1985,7,1,0,0,0, -4,&
    1983,7,1,0,0,0, -3,&
    1982,7,1,0,0,0, -2,&
    1981,7,1,0,0,0, -1,&
    0   ,0,0,0,0,0,  0 &
/),shape(leaps),order=(/2,1/))

! observation code strings 
character(2), parameter :: obscodes(60)=(/&
    '  ','1C','1P','1W','1Y', '1M','1N','1S','1L','1E',&  !  0- 9 
    '1A','1B','1X','1Z','2C', '2D','2S','2L','2X','2P',&  ! 10-19 
    '2W','2Y','2M','2N','5I', '5Q','5X','7I','7Q','7X',&  ! 20-29 
    '6A','6B','6C','6X','6Z', '6S','6L','8L','8Q','8X',&  ! 30-39 
    '2I','2Q','6I','6Q','3I', '3Q','3X','1I','1Q','5A',&  ! 40-49 
    '5B','5C','9A','9B','9C', '9X','  ','  ','  ','  ' &  ! 50-59
/)

! 1:L1/E1, 2:L2/B1, 3:L5/E5a/L3, 4:L6/LEX/B3, 5:E5b/B2, 6:E5(a+b), 7:S
integer(4), parameter :: obsfreqs(60)=(/&
    0, 1, 1, 1, 1,  1, 1, 1, 1, 1,&  !  0- 9
    1, 1, 1, 1, 2,  2, 2, 2, 2, 2,&  ! 10-19
    2, 2, 2, 2, 3,  3, 3, 5, 5, 5,&  ! 20-29
    4, 4, 4, 4, 4,  4, 4, 6, 6, 6,&  ! 30-39
    2, 2, 4, 4, 3,  3, 3, 1, 1, 3,&  ! 40-49
    3, 3, 7, 7, 7,  7, 0, 0, 0, 0 &  ! 50-59
/)

! code priority table 
! L1/E1; L2/B1; L5/E5a/L3; L6/LEX/B3; E5b/B2; E5(a+b); S
character(16), parameter :: codepris(7,MAXFREQ)=reshape((/&
    'CPYWMNSL  ','PYWCMNDSLX','IQX       ','          ','          ','          ','          ',&  ! GPS
    'PC        ','PC        ','IQX       ','          ','          ','          ','          ',&  ! GLO
    'CABXZ     ','          ','IQX       ','ABCXZ     ','IQX       ','IQX       ','          ',&  ! GAL
    'CSLXZ     ','SLX       ','IQX       ','SLX       ','          ','          ','          ',&  ! QZS
    'C         ','          ','IQX       ','          ','          ','          ','          ',&  ! SBS
    'IQX       ','IQX       ','IQX       ','IQX       ','IQX       ','          ','          ',&  ! BDS
    '          ','          ','ABCX      ','          ','          ','          ','ABCX      ' &  ! IRN
/),shape(codepris),order=(/2,1/))

! ephemeris.f90 -------------------------------------------------------------
real*8, parameter :: RE_GLO = 6378136.d0               ! radius of earth (m)            ref (2) 
real*8, parameter :: MU_GPS = 3.9860050d14             ! gravitational constant         ref (1) 
real*8, parameter :: MU_GLO = 3.9860044d14             ! gravitational constant         ref (2) 
real*8, parameter :: MU_GAL = 3.986004418d14           ! earth gravitational constant   ref (7) 
real*8, parameter :: MU_CMP = 3.986004418d14           ! earth gravitational constant   ref (9) 
real*8, parameter :: J2_GLO = 1.0826257d-3             ! 2nd zonal harmonic of geopot   ref (2) 

real*8, parameter :: OMGE_GLO = 7.292115d-5            ! earth angular velocity (rad/s) ref (2) 
real*8, parameter :: OMGE_GAL = 7.2921151467d-5        ! earth angular velocity (rad/s) ref (7) 
real*8, parameter :: OMGE_CMP = 7.292115d-5            ! earth angular velocity (rad/s) ref (9) 
real*8, parameter :: SIN_5 = -0.0871557427476582d0     ! dsin(-5.d0 deg) 
real*8, parameter :: COS_5 =  0.9961946980917456d0     ! dcos(-5.d0 deg) 

real*8, parameter :: ERREPH_GLO  = 5.d0                ! error of glonass ephemeris (m) 
real*8, parameter :: TSTEP       = 60.d0               ! integration step glonass ephemeris (s) 
real*8, parameter :: RTOL_KEPLER = 1d-13               ! relative tolerance for Kepler equation 

real*8, parameter :: DEFURASSR  = 0.15d0               ! default accurary of ssr corr (m) 
real*8, parameter :: MAXECORSSR = 10.d0                ! max orbit correction of ssr (m) 
real*8, parameter :: MAXCCORSSR = (1d-6*CLIGHT)        ! max clock correction of ssr (m) 
real*8, parameter :: MAXAGESSR  = 90.d0                ! max age of ssr orbit and clock (s) 
real*8, parameter :: MAXAGESSR_HRCLK = 10.d0           ! max age of ssr high-rate clock (s) 
real*8, parameter :: STD_BRDCCLK = 30.d0               ! error of broadcast clock (m) 
real*8, parameter :: STD_GAL_NAPA = 500.d0             ! error of galileo ephemeris for NAPA (m) 

integer*4, parameter :: MAX_ITER_KEPLER = 30           ! max number of iteration of Kelpler 
integer*4, parameter :: eph_sel(6)=(/0,0,1,0,0,0/)     ! GPS,GLO,GAL,QZS,BDS,SBS (ephemeris selections)

! pntpos.f90 ----------------------------------------------------------------
integer*4, parameter :: NX     = (4+3)                 ! # of estimated parameters 
integer*4, parameter :: MAXITR = 10                    ! max number of iteration for point pos 
real*8, parameter :: ERR_ION   = 5.d0                  ! ionospheric delay std (m) 
real*8, parameter :: ERR_TROP  = 3.d0                  ! tropspheric delay std (m) 
real*8, parameter :: ERR_SAAS  = 0.3d0                 ! saastamoinen model error std (m) 
real*8, parameter :: ERR_BRDCI = 0.5d0                 ! broadcast iono model error factor 
real*8, parameter :: ERR_CBIAS = 0.3d0                 ! code bias error std (m) 
real*8, parameter :: REL_HUMI  = 0.7d0                 ! relative humidity for saastamoinen model 

! rinex.f90 -----------------------------------------------------------------
integer*4, parameter :: NUMSYS      = 7                ! number of systems 
integer*4, parameter :: MAXRNXLEN   = (16*MAXOBSTYPE+4) ! max rinex record length 
integer*4, parameter :: MAXPOSHEAD  = 1024             ! max head line position 
integer*4, parameter :: MINFREQ_GLO = -7               ! min frequency number glonass 
integer*4, parameter :: MAXFREQ_GLO = 13               ! max frequency number glonass 
integer*4, parameter :: NINCOBS     = 262144/2         ! inclimental number of obs data 

character(7), parameter :: systyps  = 'GREJSCI'        ! satellite system codes 
character(4), parameter :: obstyps  = 'CLDS'           ! obs type codes 
character(7), parameter :: frqtyps  = '1256789'        ! frequency codes 
real*8, parameter :: ura_eph(16) = (/&                 ! ura values (ref [3] 20.3.3.3.1.1) 
    2.4,3.4,4.85,6.85,9.65,13.65,24.0,48.0,96.0,&
    192.0,384.0,768.0,1536.0,3072.0,6144.0,0.00/)
real*8, parameter :: ura_nominal(16) = (/&             ! ura nominal values 
    2.0,2.8,4.0,5.7,8.0,11.3,16.0,32.0,64.0,&
    128.0,256.0,512.0,1024.0,2048.0,4096.0,8192.0/)

! postpos.f90 ---------------------------------------------------------------
integer*4, parameter :: MAXINFILE = 1000               ! max number of input files 

! rnx2rtkp.f90 --------------------------------------------------------------
integer*4, parameter :: MAXFILE  = 16                  ! max number of input files 
integer*4, parameter :: FLNUMLIM = 200                 ! max number of rnxlist
