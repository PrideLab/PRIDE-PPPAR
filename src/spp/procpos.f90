
! Process positioning -------------------------------------------------------
subroutine procpos(fp, popt, sopt, rtk, mode, ret)  ! output pos
implicit none
include 'file_para.h'
integer*4, intent(in) :: fp
type(prcopt_t), intent(in) :: popt
type(solopt_t), intent(in) :: sopt
type(rtk_t), intent(out) :: rtk
integer*4, intent(out) :: ret  ! 0-error, 1-right
integer*4, intent(in) :: mode
type(obsd_t) obs(MAXOBS*2)  ! for rover and base 
character(MAXOBS*2*3) satrec
character(3) satid
integer*4 i,nobs,n,prn,sys,info
type(sol_t), pointer :: sol_data(:)

integer*4 fpclk1  ! output sat clock
fpclk1=FPCLK
if(clkfile_/="") open(unit=fpclk1,file=clkfile_,status='replace',iostat=info)

call rtkinit(rtk,popt);
do while(.true.)
    call inputobs(obs,rtk%sol%stat,popt,nobs)
    if(nobs<0) exit
    ! exclude satellites 
    n=1; satrec=""
    do i=1,nobs
        call satsys(obs(i)%sat,prn,sys)
        if(and(sys,popt%navsys)/=0 .and. popt%exsats(obs(i)%sat)/=1)then
            call satno2id(obs(i)%sat,satid)
            if(index(satrec,satid)==0)then
                obs(n)=obs(i); n=n+1
                satrec=trim(satrec)//satid
            endif
        endif
    enddo
    if (n<=1) cycle
    call rtkpos(rtk,obs,n-1,navs,info)
    
    if (info==0) cycle
    if(mode==0 .and. sopt%issingle==0)then ! forward/backward 
        call outsol(fp,rtk%sol,rtk%rb,sopt)
    endif
    ! if(sopt%issingle==1)then  ! store solutions
        solindex_=solindex_+1
        if(solindex_<=size(allsol_))then
            allsol_(solindex_)=rtk%sol
        else
            allocate(sol_data(size(allsol_)+MAXSOLNUM))
            sol_data(1:size(allsol_))=allsol_(1:size(allsol_))
            deallocate(allsol_)
            allsol_=>sol_data
            nullify(sol_data)
        endif
    ! endif
enddo
call rtkfree(rtk);
if(clkfile_/="") close(unit=fpclk1,status='keep')  ! output sat clock

if(solindex_==0)then
    write(*,*) "Info : no solution are computed, procpos()"
    ret=0; return
endif
ret=1
end subroutine

! initialize rtk control ----------------------------------------------------
! initialize rtk control struct
! args   : rtk_t    *rtk    IO  rtk control/result struct
!          prcopt_t *opt    I   positioning options (see rtklib.h)
! return : none
!----------------------------------------------------------------------------
subroutine rtkinit(rtk, opt)
implicit none
include 'file_para.h'
type(rtk_t), intent(out) :: rtk
type(prcopt_t), intent(in) :: opt
type(sol_t) :: sol0=sol_t(gtime_t(0,0.d0),gtime_t(0,0.d0),0.d0,0.d0,0.d0,0,0,0,0.d0)
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

! free rtk control ----------------------------------------------------------
! free memory for rtk control struct
! args   : rtk_t    *rtk    IO  rtk control/result struct
! return : none
!----------------------------------------------------------------------------
subroutine rtkfree(rtk)
implicit none
include 'file_para.h'
type(rtk_t), intent(inout) :: rtk
rtk%nx=0; rtk%na=0
deallocate(rtk%x ); nullify(rtk%x)
deallocate(rtk%P ); nullify(rtk%P)
deallocate(rtk%xa); nullify(rtk%xa)
deallocate(rtk%Pa); nullify(rtk%Pa)
end subroutine

subroutine rtkpos(rtk, obs, n, nav, stat)
implicit none
include 'file_para.h'
integer*4, intent(in) :: n
type(rtk_t), intent(inout) :: rtk
type(obsd_t), intent(in) :: obs(n)
type(nav_t), intent(in) :: nav
integer*4, intent(out) :: stat
type(gtime_t) time
integer*4 :: nu,info
character(128) :: msg
real*8 azel(n,2)
real*8, external :: timediff
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

! single-point positioning --------------------------------------------------
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
!----------------------------------------------------------------------------
subroutine pntpos(obs, n, nav, opt, sol, azel, ssat, msg, stat)
implicit none
include 'file_para.h'
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

! raim fde (failure detection and exclution) --------------------------------
subroutine raim_fde(obs, n, rs, dts, vare, svh, nav, opt, sol, azel, vsat, resp, msg, stat)
implicit none
include 'file_para.h'
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
real*8, external :: SQR

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
