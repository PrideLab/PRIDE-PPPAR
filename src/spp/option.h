
! Struct of options ---------------------------------------------------------
type prcopt_t                ! processing options type 
    sequence
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

type solopt_t                ! solution options type 
    sequence
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
    sequence
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
