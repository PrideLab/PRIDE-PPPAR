!! reference and acknowledgement to RTKLIB
!! modified by Kunlun Zhang
! constants -----------------------------------------------------------------

module pntpos_f90_
use ephe_f90_
implicit none
integer*4, parameter :: NX     = (4+3)       ! # of estimated parameters 
integer*4, parameter :: MAXITR = 10          ! max number of iteration for point pos 
real*8, parameter :: ERR_ION   = 5.d0        ! ionospheric delay std (m) 
real*8, parameter :: ERR_TROP  = 3.d0        ! tropspheric delay std (m) 
real*8, parameter :: ERR_SAAS  = 0.3d0       ! saastamoinen model error std (m) 
real*8, parameter :: ERR_BRDCI = 0.5d0       ! broadcast iono model error factor 
real*8, parameter :: ERR_CBIAS = 0.3d0       ! code bias error std (m) 
real*8, parameter :: REL_HUMI  = 0.7d0       ! relative humidity for saastamoinen model 

contains
! pseudorange measurement error variance ------------------------------------
real*8 function varerr(opt, el, sys)
implicit none
type(prcopt_t), intent(in) :: opt
real*8, intent(in) :: el
integer*4, intent(in) :: sys
real*8 fact,varr
fact=rcond(sys==SYS_GLO,EFACT_GLO,EFACT_GPS)
varr=SQR(100.d0)*(SQR(0.003d0)+SQR(0.003d0)/dsin(el))
if (opt%ionoopt==IONOOPT_IFLC) varr=varr*SQR(3.d0); ! iono-free 
varerr=SQR(fact)*varr;
end function

! get tgd parameter (m) -----------------------------------------------------
real*8 function gettgd(sat, nav)
implicit none
integer*4, intent(in) :: sat
type(nav_t), intent(in) :: nav
integer*4 i
do i=1,nav%n
    if(nav%eph(i)%sat/=sat) cycle
    gettgd=CLIGHT*nav%eph(i)%tgd(1); return
enddo
gettgd=0.d0
end function

! psendorange with code bias correction -------------------------------------
subroutine prange(obs, nav, azel, iter, opt, var, range)  !жд
implicit none
type(obsd_t), intent(in) :: obs
type(nav_t), intent(in) :: nav
real*8, intent(in) :: azel(*)
integer*4, intent(in) :: iter
type(prcopt_t), intent(in) :: opt
real*8, intent(out) :: var, range
real*8 lam(NFREQ)
real*8 PC,P1,P2,P1_P2,P1_C1,P2_C2,gamma
integer*4 :: i,j,sys,prn
i=0; j=1
P1_P2=0
P1_C1=0
P2_C2=0
var=0.d0; range=0.d0
lam=nav%lam(obs%sat,:)  ! test values here
call satsys(obs%sat,prn,sys)
if(sys==0) return

! L1-L2 for GPS/GLO/QZS, L1-L5 for GAL/SBS 
if(NFREQ>=3 .and. and(sys,or(SYS_GAL,SYS_SBS))/=0) j=2
if(NFREQ<2 .or. dabs(lam(i+1))<=1d-20 .or. dabs(lam(j+1))<=1d-20) return

gamma=SQR(lam(j+1))/SQR(lam(i+1))  ! f1^2/f2^2 
P1=obs%P(i+1)
P2=obs%P(j+1)
P1_P2=nav%cbias(obs%sat,1)
P1_C1=nav%cbias(obs%sat,2)
P2_C2=nav%cbias(obs%sat,3)

! if no P1-P2 DCB, use TGD instead 
if(dabs(P1_P2)<=1d-20 .and. and(sys,or(SYS_GPS,or(SYS_GAL,SYS_QZS)))/=0)then
    P1_P2=(1.d0-gamma)*gettgd(obs%sat,nav)
endif
if(opt%ionoopt==IONOOPT_IFLC)then            ! dual-frequency 
    if(dabs(P1)<=1d-20 .or. dabs(P2)<=1d-20) return
    if(obs%code(i+1)==CODE_L1C) P1=P1+P1_C1  ! C1%P1 
    if(obs%code(j+1)==CODE_L2C) P2=P2+P2_C2  ! C2%P2 
    PC=(gamma*P1-P2)/(gamma-1.d0)            ! iono-free combination 
else  ! single-frequency 
    if (dabs(P1)<=1d-20) return
    if (obs%code(i+1)==CODE_L1C) P1=P1+P1_C1 ! C1%P1 
    PC=P1-P1_P2/(1.d0-gamma)
endif
!if (opt%sateph==EPHOPT_SBAS) PC=PC-P1_C1  ! sbas clock based C1 
var=SQR(ERR_CBIAS)
range=PC
end subroutine

! ionospheric correction ------------------------------------------------------
! compute ionospheric correction
! args   : gtime_t time      I   time
!          nav_t  *nav       I   navigation data
!          integer*4 sat     I   satellite number
!          real*8 *pos       I   receiver position {lat,lon,h} (rad|m)
!          real*8 *azel      I   azimuth/elevation angle {az,el} (rad)
!          integer*4 ionoopt I   ionospheric correction option (IONOOPT_???)
!          real*8 *ion       O   ionospheric delay (L1) (m)
!          real*8 *var       O   ionospheric delay (L1) variance (m^2)
! return : status(1:ok,0:error)
!-----------------------------------------------------------------------------
subroutine ionocorr(time, nav, sat, pos, azel, ionoopt, ion, var, stat)
implicit none
type(gtime_t), intent(in) :: time
type(nav_t), intent(in) :: nav
integer*4, intent(in) :: sat,ionoopt
real*8, intent(in) :: pos(3),azel(2)
real*8, intent(out) :: ion,var
integer*4, intent(out) :: stat
! broadcast model 
if(ionoopt==IONOOPT_BRDC)then
    ion=ionmodel(time,nav%ion_gps,pos,azel)
    var=SQR(ion*ERR_BRDCI)
    stat=1; return
endif
ion=0.d0
var=rcond(ionoopt==IONOOPT_OFF,SQR(ERR_ION),0.d0)
stat=1
end subroutine

! tropospheric correction -----------------------------------------------------
! compute tropospheric correction
! args   : gtime_t time     I   time
!          nav_t  *nav      I   navigation data
!          real*8 *pos      I   receiver position {lat,lon,h} (rad|m)
!          real*8 *azel     I   azimuth/elevation angle {az,el} (rad)
!          integer*4 tropopt I  tropospheric correction option (TROPOPT_???)
!          real*8 *trp      O   tropospheric delay (m)
!          real*8 *var      O   tropospheric delay variance (m^2)
! return : status(1:ok,0:error)
!-----------------------------------------------------------------------------
subroutine tropcorr(time, nav, pos, azel, tropopt, trp, var, stat)
implicit none
type(gtime_t), intent(in) :: time
type(nav_t), intent(in) :: nav
real*8, intent(in) :: pos(3),azel(2)
integer*4, intent(in) :: tropopt
real*8, intent(out) :: trp, var
integer*4, intent(out) :: stat
! saastamoinen model 
if(tropopt==TROPOPT_SAAS .or. tropopt==TROPOPT_EST .or. tropopt==TROPOPT_ESTG)then
    trp=tropmodel(time,pos,azel,REL_HUMI)
    var=SQR(ERR_SAAS/(dsin(azel(2))+0.1d0))
    stat=1; return
endif
! no correction 
trp=0.d0
var=rcond(tropopt==TROPOPT_OFF,SQR(ERR_TROP),0.d0)
stat=1
end subroutine

! pseudorange residuals -----------------------------------------------------
subroutine rescode(iter, obs, n, rs, dts, vare, svh, nav, x, opt, v, H, var, vara, valres, valrag, azel, vsat, resp, ns, stat)
implicit none
integer*4, intent(in) :: iter, n, svh(n), valres(:)
type(obsd_t), intent(in) :: obs(n)
real*8, intent(in) :: rs(n,6), dts(n,2), vare(n), x(NX)
type(nav_t), intent(in) :: nav
type(prcopt_t), intent(in) :: opt
real*8, intent(out) :: v(n+4), H(n+4,NX), var(n+4), vara(n), azel(n,2), resp(n)
integer*4, intent(out) :: valrag(:), vsat(n), ns, stat
real*8 r,dion,dtrp,vmeas,vion,vtrp,rr(3),pos(3),dtr,e(3),P,lam_L1,el
integer*4 :: i,j,nv,sys,mask(4),prn,info
logical*1 isopt
nv=1; mask=0
rr(1:3)=x(1:3)
dtr=x(4)
valrag=0
call ecef2pos(rr,pos)
ns=0; i=1
do while(i<=min(n,MAXOBS))
    vsat(i)=0; azel(i,:)=0.d0; resp(i)=0.d0; vara(i)=0.d0
    call satsys(obs(i)%sat,prn,sys)
    if (sys==0)then
        i=i+1; cycle
    endif
    
    ! mark the gross error (valpos)
    if(valres(obs(i)%sat)==-1)then
        i=i+1; cycle
    endif
    
    ! reject duplicated observation data 
    if(i<n .and. i<MAXOBS)then
        if(obs(i)%sat==obs(i+1)%sat)then
            i=i+2; cycle
        endif
    endif
    
    ! geometric distance/azimuth/elevation angle 
    call geodist(rs(i,:),rr,e,r)
    call satazel(pos,e,azel(i,:),el)
    if (r<=0.d0 .or. el<opt%elmin)then
        i=i+1; cycle
    endif
    
    ! psudorange with code bias correction 
    call prange(obs(i),nav,azel(i,:),iter,opt,vmeas,P)
    if (dabs(P)<=1d-20)then
        i=i+1; valrag(i)=-1; cycle
    endif
    
    ! excluded satellite? 
    isopt=.true.
    if (satexclude(obs(i)%sat,vare(i),svh(i),opt,isopt)/=0)then
        i=i+1; cycle
    endif
    
    ! ionospheric corrections 
    call ionocorr(obs(i)%time,nav,obs(i)%sat,pos,azel(i,:),icond(iter>0,opt%ionoopt,IONOOPT_BRDC),dion,vion,info)
    if (info==0)then
        i=i+1; cycle
    endif
    
    ! GPS-L1 -> L1/B1 
    lam_L1=nav%lam(obs(i)%sat,1)
    if (lam_L1>0.d0)then
        dion=dion*SQR(lam_L1/lam_carr(1))
    endif
    
    ! tropospheric corrections 
    call tropcorr(obs(i)%time,nav,pos,azel(i,:),icond(iter>0,opt%tropopt,TROPOPT_SAAS),dtrp,vtrp,info)
    if (info==0)then
        i=i+1; cycle
    endif
    
    ! pseudorange residual 
    v(nv)=P-(r+dtr-CLIGHT*dts(i,1)+dion+dtrp)
    
    ! design matrix 
    do j=1,NX
        !H(nv,j)=rcond(j<4,-e(j),rcond(j==4,1.d0,0.d0))
        if(j<4)then
            H(nv,j)=-e(j)
        else
            H(nv,j)=rcond(j==4,1.d0,0.d0)
        endif
    enddo
    
    ! time system and receiver bias offset correction 
    if(sys==SYS_GLO)then
        v(nv)=v(nv)-x(5); H(nv,5)=1.d0; mask(2)=1
    elseif(sys==SYS_GAL)then
        v(nv)=v(nv)-x(6); H(nv,6)=1.d0; mask(3)=1
    elseif(sys==SYS_CMP)then
        v(nv)=v(nv)-x(7); H(nv,7)=1.d0; mask(4)=1
    else
        mask(1)=1
    endif
    vsat(i)=1; resp(i)=v(nv); ns=ns+1
    
    ! error variance 
    var(nv)=varerr(opt,azel(i,2),sys)+vare(i)+vmeas+vion+vtrp  ! 
    vara(i)=var(nv) ! all
    nv=nv+1; i=i+1
enddo

! constraint to avoid rank-deficient 
do i=1,4
    if(mask(i)/=0) cycle
    v(nv)=0.d0
    do j=1,NX
        H(nv,j)=rcond(j==i+3,1.d0,0.d0)
    enddo
    var(nv)=0.01d0; nv=nv+1
enddo
stat=nv-1
end subroutine

! validate solution ---------------------------------------------------------
subroutine valsol(azel, vsat, n, opt, v, nv, nx, sol, msg, stat)
implicit none
integer*4, intent(in) :: n, vsat(n), nv, nx
real*8, intent(in) :: azel(n,2), v(nv)
type(prcopt_t), intent(in) :: opt
type(sol_t), intent(out) :: sol
character(*), intent(out) :: msg
integer*4, intent(out) :: stat
real*8 azels(MAXOBS,2),dop(4),vv
integer*4 :: i,ns

i=1; ns=1
! chi-square validation of residuals 
vv=dot(v,v,nv)
if(nv>nx)then
    if(vv>chisqr(nv-nx))then
        write(msg,"('chi-square error nv=',I4,' vv=',F14.4,' cs=',F14.4)") nv,vv,chisqr(nv-nx)
        stat=0; return
        !stat=1; return
    endif
endif
! large gdop check 
do i=1,n
    if (vsat(i)==0) cycle
    azels(ns,1)=azel(i,1)
    azels(ns,2)=azel(i,2)
    ns=ns+1
enddo
call dops(ns-1,azels(1:ns-1,:),opt%elmin,dop)
sol%dop=dop
if(dop(1)<=0.d0 .or. dop(1)>opt%maxgdop)then
    write(msg,"('gdop error nv=',I4' gdop=',F10.4)") nv,dop(1)
    stat=0; return
endif
stat=1
end subroutine

!integer*4 function valres()
subroutine valpos(v, R, nv, thres, vflg)
implicit none
integer*4, intent(in) :: nv
real*8, intent(in) :: v(nv), R(nv,nv), thres
integer*4, intent(out) :: vflg(nv)
real*8 fact
integer*4 i
fact=thres*thres
do i=1, nv
    if(v(i)*v(i)<=fact*R(i,i))then
        vflg(i)=0
    else
        vflg(i)=1
    endif
enddo
end subroutine

! estimate receiver position ------------------------------------------------
subroutine estpos(obs, n, rs, dts, vare, svh, nav, opt, sol, azel, vsat, resp, msg, stat)
implicit none
integer*4, intent(in) :: n, svh(n)
type(obsd_t), intent(in) :: obs(n)
real*8, intent(in) :: rs(n,6), dts(n,2), vare(n)
type(nav_t), intent(in) :: nav
type(prcopt_t), intent(in) :: opt
type(sol_t), intent(out) :: sol
real*8, intent(out) :: azel(n,2), resp(n)
integer*4, intent(out) :: vsat(n), stat
character(*), intent(out) :: msg
real*8 :: x(NX),dx(NX),Q(NX,NX),v(n+4,1),H(n+4,NX),var(n+4,1),sig,HT(NX,n+4)
integer*4 i,j,k,info,stat1,nv,ns
real*8 thres, fact, vara(n)  ! valpos
integer*4 valres(MAXSAT), nslim, nt, nf  ! number of true / false
integer*4 navsy2,ii,prn,sy2
integer*4 valrag(MAXSAT)
thres=3d0; fact=thres*thres
x=0.d0; valres=0; nslim=0; nt=0; nf=0; valrag=0
x(1:3)=sol%rr(1:3)
do i=1,MAXITR
    ! pseudorange residuals 
100 call rescode(i-1,obs,n,rs,dts,vare,svh,nav,x,opt,v,H,var,vara,valres,valrag,azel,vsat,resp,ns,nv)
    
    nslim=4; navsy2=SYS_NONE
    do ii=1,n
        if(valrag(ii)==0)then
            call satsys(obs(ii)%sat,prn,sy2)
            navsy2=or(navsy2,sy2)
        endif
    enddo
    
    if(and(navsy2,SYS_GPS)/=0) nslim=nslim+1
    if(and(navsy2,SYS_GLO)/=0) nslim=nslim+1
    if(and(navsy2,SYS_GAL)/=0) nslim=nslim+1
    if(and(navsy2,SYS_CMP)/=0) nslim=nslim+1
    
    if (ns<nslim)then
        write(msg,"('lack of valid sats ns=',I4)") ns; exit
    endif
    ! weight by variance 
    do j=1,nv
        sig=dsqrt(var(j,1))
        v(j,1)=v(j,1)/sig
        do k=1,NX
            H(j,k)=H(j,k)/sig
        enddo
    enddo
    ! least square estimation 
    call mattran(n+4,NX,H,HT)
    call lsq(HT,v,NX,nv,dx,Q,info)
    if(info/=0)then
        write(msg,"('lsq error info=',I4)") info; exit
    endif
    do j=1,NX
        x(j)=x(j)+dx(j)
    enddo
    
    nt=0; nf=0
    do j=1,n
        if(dabs(resp(j))>=1d-20 .and. resp(j)**2> fact*vara(j)) nf=nf+1
        if(dabs(resp(j))>=1d-20 .and. resp(j)**2<=fact*vara(j)) nt=nt+1
    enddo
    !if(norm(dx,NX)<1d-4)then  !valpos
    if(0<nf .and. nf<=4 .and. nf<nt)then  !valpos
        do j=1,n
            if(dabs(resp(j))>=1d-20 .and. resp(j)**2>fact*vara(j))then
                valres(obs(j)%sat)=1
            endif
        enddo
        do j=1,n
            if(valres(obs(j)%sat)==1)then
                valres(obs(j)%sat)=-1
            endif
        enddo
        goto 100
    endif
    if(norm(dx,NX)<1d-4)then
        sol%mytype=0
        sol%time=timeadd(obs(1)%time,-x(4)/CLIGHT)
        sol%dtr(1)=x(4)/CLIGHT  ! receiver clock bias (s) 
        sol%dtr(2)=x(5)/CLIGHT  ! glo-gps time offset (s) 
        sol%dtr(3)=x(6)/CLIGHT  ! gal-gps time offset (s) 
        sol%dtr(4)=x(7)/CLIGHT  ! bds-gps time offset (s) 
        do j=1,6
            sol%rr(j)=rcond(j<=3,x(j),0.d0)
        enddo
        do j=1,3
            sol%qr(j)=Q(j,j)
        enddo
        sol%qr(4)=Q(1,2)     ! cov xy 
        sol%qr(5)=Q(2,3)     ! cov yz 
        sol%qr(6)=Q(1,3)     ! cov zx 
        sol%ns=ns
        ! validate solution 
        call valsol(azel,vsat,n,opt,v,nv,NX,sol,msg,stat1)
        sol%stat=icond(opt%sateph==EPHOPT_SBAS,SOLQ_SBAS,SOLQ_SINGLE)
        stat=stat1; return
    endif
enddo
if (i>MAXITR) write(msg,"('iteration divergent i=',I4)") i-1
stat=0
end subroutine

subroutine init_sol(sol)
implicit none
type(sol_t), intent(out) :: sol
sol=sol_t(gtime_t(0,0.d0),gtime_t(0,0.d0),0.d0,0.d0,0.d0,0,0,0,0.d0)
end subroutine

! raim fde (failure detection and exclution) -------------------------------
subroutine raim_fde(obs, n, rs, dts, vare, svh, nav, opt, sol, azel, vsat, resp, msg, stat)
implicit none
integer*4, intent(in) :: n, svh(n)
type(obsd_t), intent(in) :: obs(n)
real*8, intent(in) :: rs(n,6), dts(n,2), vare(n)
type(nav_t), intent(in) :: nav
type(prcopt_t), intent(in) :: opt
type(sol_t), intent(out) :: sol
real*8, intent(out) :: azel(n,2), resp(n)
integer*4, intent(out) :: vsat(n), stat
character(*), intent(out) :: msg

type(obsd_t), pointer :: obs_e(:)
type(sol_t) sol_e
character(32) tstr
character(16) name
character(128) msg_e
real*8 :: rs_e(n,6),dts_e(n,2),vare_e(n),azel_e(n,2),resp_e(n),rms_e,rms
integer*4 :: i,j,k,nvsat,stat1,svh_e(n),vsat_e(n),sat,info

rms=100.d0; stat1=0; sat=0
allocate(obs_e(n))
if(size(obs_e,dim=1)<n)then
    stat=0; return
endif
call init_sol(sol_e)
do i=1,n
    ! satellite exclution 
    j=1; k=1
    do j=1,n
        if (j==i) cycle
        obs_e(k)=obs(j)
        rs_e(k,:)=rs(j,:)
        dts_e(k,:)=dts(j,:)
        vare_e(k)=vare(j)
        svh_e(k)=svh(j)
        k=k+1
    enddo
    ! estimate receiver position without a satellite 
    call estpos(obs_e,n-1,rs_e,dts_e,vare_e,svh_e,nav,opt,sol_e,azel_e,vsat_e,resp_e,msg_e,info)
    if (info==0) cycle
    
    j=1; nvsat=0; rms_e=0.d0
    do j=1,n-1
        if(vsat_e(j)==0) cycle
        rms_e=rms_e+SQR(resp_e(j))
        nvsat=nvsat+1
    enddo
    if (nvsat<5) cycle
    rms_e=dsqrt(rms_e/nvsat)
    if (rms_e>rms) cycle
    
    ! save result 
    j=1; k=1
    do j=1,n
        if (j==i) cycle
        azel(j,:)=azel_e(k,:)
        vsat(j)=vsat_e(k)
        resp(j)=resp_e(k)
        k=k+1
    enddo
    stat1=1
    sol=sol_e
    sat=obs(i)%sat
    rms=rms_e
    vsat(i)=0
    msg=msg_e
enddo
if(stat1/=0)then
    call time2str(obs(1)%time,tstr,2)
    call satno2id(sat,name)
endif
deallocate(obs_e)
stat=stat1
end subroutine

! single-point positioning ----------------------------------------------------
! compute receiver position, velocity, clock bias by single-point positioning
! with pseudorange and doppler observables
! args   : obsd_t *obs      I   observation data
!          integer*4 n      I   number of observation data
!          nav_t  *nav      I   navigation data
!          prcopt_t *opt    I   processing options
!          sol_t  *sol      IO  solution
!          real*8 *azel     IO  azimuth/elevation angle (rad) (NULL: no output)
!          ssat_t *ssat     IO  satellite status              (NULL: no output)
!          char   *msg      O   error message for error exit
! return : status(1:ok,0:error)
! notes  : assuming sbas-gps, galileo-gps, qzss-gps, compass-gps time offset and
!          receiver bias are negligible (only involving glonass-gps time offset
!          and receiver bias)
!-----------------------------------------------------------------------------
subroutine pntpos(obs, n, nav, opt, sol, azel, ssat, msg, stat)
implicit none
integer*4, intent(in) :: n
type(obsd_t), intent(in) :: obs(n)
type(nav_t), intent(in) :: nav
type(prcopt_t), intent(in) :: opt
type(sol_t), intent(out) :: sol
real*8, intent(out) :: azel(n,2)
type(ssat_t), intent(out) :: ssat(MAXSAT)
character(*), intent(out) :: msg
integer*4, intent(out) :: stat
type(prcopt_t) opt_
real*8 rs(n,6),dts(n,2),var(n),azel_(n,2),resp(n),clksod
integer*4 :: i,stat1,vsat(MAXOBS),svh(MAXOBS),clkmjd
character(3) clkprn
vsat=0; opt_=opt
sol%stat=SOLQ_NONE
if(n<=0)then
    msg='no observation data'; stat=0; return
endif
sol%time0=obs(1)%time  ! save timestamp of observation
sol%time=obs(1)%time; msg=''

if(opt_%mode/=PMODE_SINGLE)then  ! for precise positioning 
    opt_%ionoopt=IONOOPT_BRDC
    opt_%tropopt=TROPOPT_SAAS
endif
! satellite positons, velocities and clocks 
call satposs(sol%time,obs,n,nav,opt_%sateph,rs,dts,var,svh)

! output sat clock bias
if(clkfile_/="")then
    do i=1,n
        call time2mjd(sol%time, clkmjd, clksod)
        call satno2id(obs(i)%sat, clkprn)
        !write(FPCLK,"(F12.4\)") sol%time, clkprn, dts(i,1), dts(i,1) ! nowrap
        write(FPCLK,*) clkmjd, clksod, clkprn, dts(i,1), dts(i,2)
    enddo
endif

! estimate receiver position with pseudorange 
call estpos(obs,n,rs,dts,var,svh,nav,opt_,sol,azel_,vsat,resp,msg,stat1)

! raim fde 
if (stat1==0 .and. n>=6 .and. opt%posopt(5)/=0)then
    call raim_fde(obs,n,rs,dts,var,svh,nav,opt_,sol,azel_,vsat,resp,msg,stat1)
endif

azel=azel_
do i=1,MAXSAT
    ssat(i)%vs=0
    ssat(i)%azel(1:2)=0.d0
    ssat(i)%resp(1)=0.d0
    ssat(i)%resc(1)=0.d0
    ssat(i)%snr(1)=0
enddo
do i=1,n
    ssat(obs(i)%sat)%azel(1)=azel_(i,1)
    ssat(obs(i)%sat)%azel(2)=azel_(i,2)
    ssat(obs(i)%sat)%snr(1)=obs(i)%SNR(1)
    if (vsat(i)==0) cycle
    ssat(obs(i)%sat)%vs=1
    ssat(obs(i)%sat)%resp(1)=resp(i)
enddo
stat=stat1
end subroutine
end module pntpos_f90_
