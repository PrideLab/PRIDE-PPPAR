
! Struct of navigation ------------------------------------------------------
type eph_t                   ! GPS/QZS/GAL broadcast ephemeris type 
    sequence
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
    sequence
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

type alm_t                   ! almanac type 
    sequence
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

type nav_t                   ! navigation data type 
    sequence
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
