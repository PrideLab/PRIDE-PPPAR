!! reference and acknowledgement to RTKLIB
!! modified by Kunlun Zhang
! constants/macros ----------------------------------------------------------

module rtkpos_f90_
use pntpos_f90_
implicit none

real*8, parameter :: VAR_POS    = 30.d0**2   ! initial variance of receiver pos (m^2) 
real*8, parameter :: VAR_VEL    = 10.d0**2   ! initial variance of receiver vel ((m/s)^2) 
real*8, parameter :: VAR_ACC    = 10.d0**2   ! initial variance of receiver acc ((m/ss)^2) 
real*8, parameter :: VAR_HWBIAS = 1.d0**2    ! initial variance of h/w bias ((m/MHz)^2) 
real*8, parameter :: VAR_GRA    = 1d-3**2    ! initial variance of gradient (m^2) 
real*8, parameter :: INIT_ZWD   = 0.15d0     ! initial zwd (m) 

real*8, parameter :: PRN_HWBIAS = 1d-6       ! process noise of h/w bias (m/MHz/dsqrt(s)) 
integer*4, parameter :: GAP_RESION = 120     ! gap to reset ionosphere parameters (epochs) 
real*8, parameter :: MAXACC     = 30.d0      ! max accel for doppler slip detection (m/s^2) 
real*8, parameter :: VAR_HOLDAMB = 1d-3      ! constraint to hold ambiguity (cycle^2) 
real*8, parameter :: TTOL_MOVEB = (1.d0+2*DTTOL) ! time sync tolerance for moving-baseline (s) 

contains
!! number of parameters (pos,ionos,tropos,hw-bias,phase-bias,real,estimated)
!define NF(opt)     ((opt)%ionoopt==IONOOPT_IFLC?1:(opt)%nf)
!define NP(opt)     ((opt)%dynamics==0?3:9)
!define NI(opt)     ((opt)%ionoopt/=IONOOPT_EST?0:MAXSAT)
!define NT(opt)     ((opt)%tropopt<TROPOPT_EST?0:((opt)%tropopt<TROPOPT_ESTG?2:6))
!define NL(opt)     ((opt)%glomodear/=GLO_ARMODE_AUTOCAL?0:NFREQGLO)
!define NB(opt)     ((opt)%mode<=PMODE_DGPS?0:MAXSAT*NF(opt))
!define NR(opt)     (NP(opt)+NI(opt)+NT(opt)+NL(opt))
!define NX(opt)     (NR(opt)+NB(opt))
!! state variable index 
!define II(s,opt)   (NP(opt)+(s)-1)                 ! ionos (s:satellite no) 
!define IT(r,opt)   (NP(opt)+NI(opt)+NT(opt)/2*(r)) ! tropos (r:0=rov,1:ref) 
!define IL(f,opt)   (NP(opt)+NI(opt)+NT(opt)+(f))   ! receiver h/w bias 
!define IB(s,f,opt) (NR(opt)+MAXSAT*(f)+(s)-1) ! phase bias (s:satno,f:freq) 

! open solution status file ---------------------------------------------------
! open solution status file and set output level
! args   : char     *file   I   rtk status file
!          integer*4      level   I   rtk status level (0: off)
! return : status (1:ok,0:error)
! notes  : file can constain time keywords (%Y,%y,%m...) defined in reppath().
!          The time to replace keywords is based on UTC of CPU time.
! output : solution status file record format
!
!   $POS,week,tow,stat,posx,posy,posz,posxf,posyf,poszf
!          week/tow : gps week no/time of week (s)
!          stat     : solution status
!          posx/posy/posz    : position x/y/z ecef (m) float
!          posxf/posyf/poszf : position x/y/z ecef (m) fixed
!
!   $VELACC,week,tow,stat,vele,veln,velu,acce,accn,accu,velef,velnf,veluf,accef,accnf,accuf
!          week/tow : gps week no/time of week (s)
!          stat     : solution status
!          vele/veln/velu    : velocity e/n/u (m/s) float
!          acce/accn/accu    : acceleration e/n/u (m/s^2) float
!          velef/velnf/veluf : velocity e/n/u (m/s) fixed
!          accef/accnf/accuf : acceleration e/n/u (m/s^2) fixed
!
!   $CLK,week,tow,stat,clk1,clk2,clk3,clk4
!          week/tow : gps week no/time of week (s)
!          stat     : solution status
!          clk1     : receiver clock bias GPS (ns)
!          clk2     : receiver clock bias GLO-GPS (ns)
!          clk3     : receiver clock bias GAL-GPS (ns)
!          clk4     : receiver clock bias BDS-GPS (ns)
!          cmn_bias : common phase bias removed from all states
!
!   $ION,week,tow,stat,sat,az,el,ion,ion-fixed
!          week/tow : gps week no/time of week (s)
!          stat     : solution status
!          sat      : satellite id
!          az/el    : azimuth/elevation angle(deg)
!          ion      : vertical ionospheric delay L1 (m) float
!          ion-fixed: vertical ionospheric delay L1 (m) fixed
!
!   $TROP,week,tow,stat,rcv,ztd,ztdf
!          week/tow : gps week no/time of week (s)
!          stat     : solution status
!          rcv      : receiver (1:rover,2:base station)
!          ztd      : zenith total delay (m) float
!          ztdf     : zenith total delay (m) fixed
!
!   $HWBIAS,week,tow,stat,frq,bias,biasf
!          week/tow : gps week no/time of week (s)
!          stat     : solution status
!          frq      : frequency (1:L1,2:L2,...)
!          bias     : h/w bias coefficient (m/MHz) float
!          biasf    : h/w bias coefficient (m/MHz) fixed
!
!   $SAT,week,tow,sat,frq,az,el,resp,resc,vsat,snr,fix,slip,lock,outc,slipc,rejc
!          week/tow : gps week no/time of week (s)
!          sat/frq  : satellite id/frequency (1:L1,2:L2,...)
!          az/el    : azimuth/elevation angle (deg)
!          resp     : pseudorange residual (m)
!          resc     : carrier-phase residual (m)
!          vsat     : valid data flag (0:invalid,1:valid)
!          snr      : signal strength (dbHz)
!          fix      : ambiguity flag  (0:no data,1:float,2:fixed,3:hold,4:ppp)
!          slip     : cycle-slip flag (bit1:slip,bit2:parity unknown)
!          lock     : carrier-lock count
!          outc     : data outage count
!          slipc    : cycle-slip count
!          rejc     : data reject (outlier) count
!          icbias   : interchannel bias (GLONASS)
!          bias     : phase bias 
!          bias_var : variance of phase bias
!          lambda   : wavelength
!
!-----------------------------------------------------------------------------

! initialize rtk control -----------------------------------------------------
! initialize rtk control struct
! args   : rtk_t    *rtk    IO  rtk control/result struct
!          prcopt_t *opt    I   positioning options (see rtklib.h)
! return : none
!-----------------------------------------------------------------------------
subroutine rtkinit(rtk, opt)
implicit none
type(rtk_t), intent(out) :: rtk
type(prcopt_t), intent(in) :: opt
type(sol_t) :: sol0=sol_t(gtime_t(0,0.d0),0.d0,0.d0,0.d0,0,0,0)
type(ssat_t) :: ssat0=ssat_t(0,0,0.d0,0.d0,0.d0,0.d0,0,0,0)
integer*4 i

rtk%sol=sol0; rtk%rb=0.d0; rtk%tt=0.d0
rtk%nx=0; rtk%na=0
allocate(rtk%x(rtk%nx))
allocate(rtk%P(rtk%nx,rtk%nx))
allocate(rtk%xa(rtk%na))
allocate(rtk%Pa(rtk%na,rtk%na))
rtk%x=0.d0; rtk%P=0.d0; rtk%xa=0.d0; rtk%Pa=0.d0
rtk%nfix=0;     rtk%neb=0
rtk%ssat=ssat0; rtk%excsat=0
rtk%opt=opt;    rtk%initial_mode=opt%mode
end subroutine

! free rtk control ------------------------------------------------------------
! free memory for rtk control struct
! args   : rtk_t    *rtk    IO  rtk control/result struct
! return : none
!-----------------------------------------------------------------------------
subroutine rtkfree(rtk)
implicit none
type(rtk_t), intent(inout) :: rtk
rtk%nx=0; rtk%na=0
deallocate(rtk%x ); nullify(rtk%x)
deallocate(rtk%P ); nullify(rtk%P)
deallocate(rtk%xa); nullify(rtk%xa)
deallocate(rtk%Pa); nullify(rtk%Pa)
end subroutine

! precise positioning ---------------------------------------------------------
! input observation data and navigation message, compute rover position by 
! precise positioning
! args   : rtk_t *rtk       IO  rtk control/result struct
!            rtk%sol        IO  solution
!                .time      O   solution time
!                .rr()      IO  rover position/velocity
!                               (I:fixed mode,O:single mode)
!                .dtr(0)    O   receiver clock bias (s)
!                .dtr(1)    O   receiver glonass-gps time offset (s)
!                .Qr()      O   rover position covarinace
!                .stat      O   solution status (SOLQ_???)
!                .ns        O   number of valid satellites
!                .age       O   age of differential (s)
!                .ratio     O   ratio factor for ambiguity validation
!            rtk%rb()       IO  base station position/velocity
!                               (I:relative mode,O:moving-base mode)
!            rtk%nx         I   number of all states
!            rtk%na         I   number of integer states
!            rtk%ns         O   number of valid satellite
!            rtk%tt         O   time difference between current and previous (s)
!            rtk%x()        IO  float states pre-filter and post-filter
!            rtk%P()        IO  float covariance pre-filter and post-filter
!            rtk%xa()       O   fixed states after AR
!            rtk%Pa()       O   fixed covariance after AR
!            rtk%ssat(s)    IO  sat(s+1) status
!                .sys       O   system (SYS_???)
!                .az   (r)  O   azimuth angle   (rad) (r=0:rover,1:base)
!                .el   (r)  O   elevation angle (rad) (r=0:rover,1:base)
!                .vs   (r)  O   data valid single     (r=0:rover,1:base)
!                .resp (f)  O   freq(f+1) pseudorange residual (m)
!                .resc (f)  O   freq(f+1) carrier-phase residual (m)
!                .vsat (f)  O   freq(f+1) data vaild (0:invalid,1:valid)
!                .fix  (f)  O   freq(f+1) ambiguity flag
!                               (0:nodata,1:float,2:fix,3:hold)
!                .slip (f)  O   freq(f+1) slip flag
!                               (bit8-7:rcv1 LLI, bit6-5:rcv2 LLI,
!                                bit2:parity unknown, bit1:slip)
!                .lock (f)  IO  freq(f+1) carrier lock count
!                .outc (f)  IO  freq(f+1) carrier outage count
!                .slipc(f)  IO  freq(f+1) cycle slip count
!                .rejc (f)  IO  freq(f+1) data reject count
!                .gf        IO  geometry-free phase (L1-L2) (m)
!                .gf2       IO  geometry-free phase (L1-L5) (m)
!            rtk%nfix       IO  number of continuous fixes of ambiguity
!            rtk%neb        IO  bytes of error message buffer
!            rtk%errbuf     IO  error message buffer
!            rtk%tstr       O   time string for debug
!            rtk%opt        I   processing options
!          obsd_t *obs      I   observation data for an epoch
!                               obs(i).rcv=1:rover,2:reference
!                               sorted by receiver and satellte
!          integer*4 n      I   number of observation data
!          nav_t  *nav      I   navigation messages
! return : status (0:no solution,1:valid solution)
! notes  : before calling function, base station position rtk%sol.rb(_) should
!          be properly set for relative mode except for moving-baseline
!-----------------------------------------------------------------------------
subroutine rtkpos(rtk, obs, n, nav, stat)
implicit none
integer*4, intent(in) :: n
type(rtk_t), intent(inout) :: rtk
type(obsd_t), intent(in) :: obs(n)
type(nav_t), intent(in) :: nav
integer*4, intent(out) :: stat
type(gtime_t) time
integer*4 :: nu,info
character(128) :: msg
real*8 azel(n,2)
nu=1; msg=''
! count rover/base station observations 
do while(nu<=n .and. obs(nu)%rcv==1)
    nu=nu+1
    if(nu>n) exit
enddo
nu=nu-1
time=rtk%sol%time  ! previous epoch 
! rover position by single point positioning 
call pntpos(obs,nu,nav,rtk%opt,rtk%sol,azel,rtk%ssat,msg,info)

if (info==0)then
    stat=0; return
endif
if (time%time/=0)then
    rtk%tt=timediff(rtk%sol%time,time)
endif
stat=1
end subroutine
end module rtkpos_f90_
