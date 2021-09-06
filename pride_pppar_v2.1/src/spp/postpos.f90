!! reference and acknowledgement to RTKLIB
!! modified by Kunlun Zhang
! constants -----------------------------------------------------------------

module postpos_f90_
use solution_f90_
use rtkpos_f90_
use rinex_f90_
implicit none

integer*4, parameter :: MAXPRCDAYS = 100          ! max days of continuous processing 
integer*4, parameter :: MAXINFILE  = 1000         ! max number of input files 

! constants/global variables ------------------------------------------------
private obss, navs, stas, nepoch, iobsu
private iobsr, isbs, ilex, revs, aborts
type(obs_t) :: obss=obs_t(0,0,0,null())           ! observation data 
type(nav_t) :: navs=nav_t(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
                        null(),null(),null(),0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,&
                        0.d0,0.d0,0.d0,0.d0,0,0.d0,0.d0,0.d0,0.d0,0.d0,'')  ! navigation data 
type(sta_t) stas(MAXRCV)        ! station infomation 
integer*4 :: nepoch = 0         ! number of observation epochs 
integer*4 :: iobsu  = 0         ! current rover observation data index 
integer*4 :: iobsr  = 0         ! current reference observation data index 
integer*4 :: isbs   = 0         ! current sbas message index 
integer*4 :: ilex   = 0         ! current lex message index 
integer*4 :: revs   = 0         ! analysis direction (0:forward,1:backward) 
integer*4 :: aborts = 0         ! abort status 

type(sol_t) solf                ! forward solutions 
type(sol_t) solb                ! backward solutions 
real*8 rbf(3)                   ! forward base positions 
real*8 rbb(3)                   ! backward base positions 
integer*4 :: isolf=0            ! current forward solutions index 
integer*4 :: isolb=0            ! current backward solutions index 
character :: proc_rov (64)=''   ! rover for current processing 
character :: proc_base(64)=''   ! base station for current processing 
character :: rtcm_file(1024)='' ! rtcm data file 
character :: rtcm_path(1024)='' ! rtcm data path 
!rtcm_t rtcm                    ! rtcm control struct 
!FILE *fp_rtcm=NULL             ! rtcm data file pointer 

contains
! output header -------------------------------------------------------------
subroutine outheader(fp, file, n, popt, sopt)
implicit none
integer*4, intent(in) :: fp, n
character(*), intent(in) :: file(n)
type(prcopt_t), intent(in) :: popt
type(solopt_t), intent(in) :: sopt
character(4) :: s1(3)
type(gtime_t) ts,te
real*8 t1,t2
integer*4 i,j,w1,w2
character(32) s2,s3

s1=(/'GPST','UTC ','JST '/)
if (sopt%posf==SOLF_NMEA .or. sopt%posf==SOLF_STAT) return
if(sopt%outhead/=0)then
     
    if (len_trim(sopt%prog)<=0)then
        write(fp,"(A,' program   : RTKLIB ver.',A)") COMMENTH,VER_RTKLIB
    else
        write(fp,"(A,' program   : ',A)") COMMENTH,sopt%prog
    endif
    
    do i=1,n
        write(fp,"(A,' inp file  : ',A)") COMMENTH,trim(file(i))
    enddo
    do i=1,obss%n
        if(obss%mydata(i)%rcv==1) exit
    enddo
    do j=obss%n,1,-1
        if(obss%mydata(j)%rcv==1) exit
    enddo
    if(j<i)then
        write(fp,"(/A,' no rover obs data')") COMMENTH; return
    endif
    ts=obss%mydata(i)%time
    te=obss%mydata(j)%time
    call time2gpst(ts,w1,t1)
    call time2gpst(te,w2,t2)
    if (sopt%times>=1) ts=gpst2utc(ts)
    if (sopt%times>=1) te=gpst2utc(te)
    if (sopt%times==2) ts=timeadd(ts,9*3600.d0)
    if (sopt%times==2) te=timeadd(te,9*3600.d0)
    call time2str(ts,s2,1)
    call time2str(te,s3,1)
    write(fp,"(A,' obs start : ',A,' ',A,' (week',I4.4,' ',F8.1,'s)')") COMMENTH,trim(s2),trim(s1(sopt%times+1)),w1,t1
    write(fp,"(A,' obs end   : ',A,' ',A,' (week',I4.4,' ',F8.1,'s)')") COMMENTH,trim(s3),trim(s1(sopt%times+1)),w2,t2
    if(sopt%outopt/=0) call outprcopt(fp,popt)
    if(sopt%outhead/=0 .or. sopt%outopt/=0) write(fp,"(A)") COMMENTH
    call outsolhead(fp,sopt)
endif
end subroutine

! output header (custom format 1) -------------------------------------------
subroutine outheader2(fp, sta, obs)
implicit none
integer*4, intent(in) :: fp
type(sta_t), intent(in) :: sta
type(obs_t), intent(in) :: obs
if(fp/=6)then  ! disallow output to screen
    write(fp,"(A20,A10,A4,A26,A7)") "Kinematic Trajectory",'',trim(sta%name(1:4)),'',"COMMENT"
    write(fp,"(F9.2,A51,A8)") obs%tint,'',"INTERVAL"
    write(fp,"(A60,A13)") '',"END OF HEADER"
endif
end subroutine

! search next observation data index ----------------------------------------
subroutine nextobsf(obs, i, rcv, stat)
implicit none
type(obs_t), intent(in) :: obs
integer*4, intent(out) :: i, stat
integer*4, intent(in) :: rcv
real*8 tt
integer*4 n
n=0
do while(i<obs%n)
    if(obs%mydata(i+1)%rcv==rcv) exit
    i=i+1
enddo
do while(i+n<obs%n)
    tt=timediff(obs%mydata(i+1+n)%time,obs%mydata(i+1)%time)
    if(obs%mydata(i+1+n)%rcv/=rcv .or. tt>DTTOL) exit
    n=n+1
enddo
stat=n
end subroutine

subroutine nextobsb(obs, i, rcv, stat)
implicit none
type(obs_t), intent(in) :: obs
integer*4, intent(out) :: i, stat
integer*4, intent(in) :: rcv
real*8 tt
integer*4 n
n=0
do while(i>=0)
    if(obs%mydata(i+1)%rcv==rcv) exit
    i=i-1
enddo
do while(i-n>=0)
    tt=timediff(obs%mydata(i+1-n)%time,obs%mydata(i+1)%time)
    if(obs%mydata(i+1-n)%rcv/=rcv .or. tt>DTTOL) exit
    n=n+1
enddo
stat=n
end subroutine

! input obs data, navigation messages and sbas correction -------------------
subroutine inputobs(obs, solq, popt, stat)
implicit none
type(obsd_t), intent(out) :: obs(:)
integer*4, intent(in) :: solq
type(prcopt_t), intent(in) :: popt
integer*4, intent(out) :: stat
type(gtime_t) :: time
integer*4 :: i,nu,n
time=gtime_t(0,0.d0)
n=0
if(revs==0)then  ! input forward data 
    call nextobsf(obss,iobsu,1,nu)
    if (nu<=0)then
        stat=-1; return
    endif
    i=0
    do while(i<nu .and. n<MAXOBS*2)
        !if(obss%mydata(iobsu+1+i)%sat<33)then
        if(obss%mydata(iobsu+1+i)%sat<=MAXSAT)then
            obs(n+1)=obss%mydata(iobsu+1+i); n=n+1
        endif
        i=i+1
    enddo
    iobsu=iobsu+nu
endif
stat=n
end subroutine

! process positioning -------------------------------------------------------
subroutine procpos(fp, popt, sopt, rtk, mode, ret)  ! output pos
implicit none
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

!integer*4 fpres  ! output res
!fpres=FPADD
!open(unit=fpres,file='D:\Desktop\data\out.res',status='replace',iostat=info)

integer*4 fpclk1  ! output sat clock
fpclk1=FPCLK
if(clkfile_/="") open(unit=fpclk1,file=clkfile_,status='replace',iostat=info)

!integer*4 fpclk2  ! output rcv clock
!fpclk2=FPCLK
!open(unit=fpclk2,file='D:\Desktop\out.clk',status='replace',iostat=info)

!integer*4 fppos1  ! output sat position
!fppos1=FPPOS
!open(unit=fppos1,file='D:\Desktop\out.pos',status='replace',iostat=info)

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
    
    !do i=1,MAXSAT  !add res
    !    write(fpres,"(F12.4\)") rtk%ssat(i)%resp(1)
    !enddo
    !write(fpres,"(A\)") char(10)  ! \n]
    !write(fpclk2,*) rtk%sol%dtr(:)  ! rcv clock
enddo
!close(unit=fpres, status='keep')  ! output res
if(clkfile_/="") close(unit=fpclk1,status='keep')  ! output sat clock
!close(unit=fpclk2,status='keep')  ! output rcv clock
!close(unit=fppos1,status='keep')  ! output sat position

if(solindex_==0)then
    ret=0; return
endif
ret=1
end subroutine

! print the final result
subroutine printresult(sopt, ret)
implicit none
type(solopt_t), intent(in) :: sopt
integer*4, intent(out) :: ret  ! 0-error, 1-right
real*8, allocatable :: rr(:,:)
real*8 avexyz(3)
type(gtime_t) tss, tee
real*8 epss(6), timespan
avexyz=0
! if(sopt%issingle==1)then  ! output a single solution
    if(solindex_==0)then
        ret=0; return
    endif
    
    allocate(rr(solindex_,3))
    call errorfilter(allsol_(1:solindex_)%rr(1),rr(:,1))
    call errorfilter(allsol_(1:solindex_)%rr(2),rr(:,2))
    call errorfilter(allsol_(1:solindex_)%rr(3),rr(:,3))
    avexyz(1)=raver(rr(:,1))
    avexyz(2)=raver(rr(:,2))
    avexyz(3)=raver(rr(:,3))
    deallocate(rr)
    
    if(avexyz(1)==0 .and. avexyz(2)==0 .and. avexyz(3)==0)then
        ret=0; return
    endif
    !if(sopt%issingle==1)then
        write(unit=6,fmt="(A11,3F16.4)") "Position : ", avexyz(1), avexyz(2), avexyz(3)
        !write(unit=6,fmt="(A5,3F16.4)") trim(stas(1)%name(1:4)), raver(rr(:,1)), raver(rr(:,2)), raver(rr(:,3))
    !endif
    ! output time stamp
    tss=allsol_(1)%time0
    tee=allsol_(solindex_)%time0
    call time2epoch(tss, epss)
    timespan=timediff(tee, tss)
    write(unit=6,fmt="(A11,I4.4,' ',I2.2,' ',I2.2,' ',I2.2,' ',I2.2,' ',F5.2,' 'F11.2)") &
          "Duration : ", int(epss(1)), int(epss(2)), int(epss(3)), &
          int(epss(4)), int(epss(5)), epss(6), timespan  ! , round(timespan)
! endif
ret=1
end subroutine

! read obs and nav data -----------------------------------------------------
subroutine readobsnav(ts_in, te_in, ti, infile, myindex, n, prcopt, obs, nav, sta, stat)
implicit none
type(gtime_t), intent(in) :: ts_in, te_in
real*8, intent(in) :: ti
character(*), intent(in) :: infile(:)
integer*4, intent(in) :: myindex(:), n
type(prcopt_t), intent(in) :: prcopt
type(obs_t), intent(out) :: obs
type(nav_t), intent(out) :: nav
type(sta_t), intent(out) :: sta(:)
integer*4, intent(out) :: stat
integer*4 :: i,j,ind,nobs,rcv,info
type(sta_t) statmp
type(gtime_t) ts, te
ts=ts_in; te=te_in
ind=0; nobs=0; rcv=1
nullify(obs%mydata); obs%n =0; obs%nmax =0
nullify(nav%eph);    nav%n =0; nav%nmax =0
nullify(nav%geph);   nav%ng=0; nav%ngmax=0
nepoch=0
do i=1,n
    if(myindex(i)/=ind)then
        if(obs%n>nobs) rcv=rcv+1
        ind=myindex(i); nobs=obs%n
    endif
    ! read rinex obs and nav file 
    if(rcv<=2)then
        call readrnxt(infile(i),rcv,ts,te,ti,prcopt%rnxopt(icond(rcv<=1,1,2)),obs,nav,sta(rcv),info)
    else
        call readrnxt(infile(i),rcv,ts,te,ti,prcopt%rnxopt(icond(rcv<=1,1,2)),obs,nav,statmp,info)
    endif
    if (info<0)then
        stat=0; return
    endif
enddo
if(obs%n<=0)then
    stat=0; return
endif
if(nav%n<=0 .and. nav%ng<=0 .and. nav%ns<=0)then
    stat=0; return
endif
! sort observation data 
nepoch=obs%n

! delete duplicated ephemeris 
call uniqnav(nav)

! set time span for progress display 
if(ts%time==0 .or. te%time==0)then 
    do i=1,obs%n
        if (obs%mydata(i)%rcv==1) exit
    enddo
    do j=obs%n,1,-1
        if (obs%mydata(j)%rcv==1) exit
    enddo
    if(i<j)then
        if (ts%time==0) ts=obs%mydata(i)%time
        if (te%time==0) te=obs%mydata(j)%time
        !settspan(ts,te)
    endif
endif
stat=1
end subroutine

! free obs and nav data -----------------------------------------------------
subroutine freeobsnav(obs, nav)
implicit none
type(obs_t), intent(out) :: obs
type(nav_t), intent(out) :: nav
deallocate(obs%mydata); nullify(obs%mydata); obs%n =0; obs%nmax =0
if(associated(nav%eph))then
    deallocate(nav%eph)
    nullify(nav%eph)
    nav%n =0; nav%nmax =0
endif
if(associated(nav%geph))then
    deallocate(nav%geph)
    nullify(nav%geph)
    nav%ng=0; nav%ngmax=0
endif
end subroutine

! write header to output file -----------------------------------------------
integer*4 function outhead(outfile, infile, n, popt, sopt, sta, obs)    ! change the output format here_1
implicit none
character(*), intent(in) :: outfile, infile(:)
integer*4, intent(in) :: n
type(prcopt_t), intent(in) :: popt
type(solopt_t), intent(in) :: sopt
type(sta_t), intent(in) :: sta
type(obs_t), intent(in) :: obs
integer*4 :: fp=6, info
if(len_trim(outfile)/=0)then
    fp=FPOUT
    open(unit=fp,file=trim(outfile),status='replace',iostat=info)
    if(info/=0)then
        !showmsg('error : open output file %s',outfile)
        write(unit=6,fmt="('error : open output file ',A)") trim(outfile)
        outhead=0; return
    endif
endif
! output header 
!call outheader(fp,infile,n,popt,sopt)
call outheader2(fp,sta,obs)
if(len_trim(outfile)/=0) close(unit=fp,status='keep')
outhead=1
end function

! open output file for append -----------------------------------------------
integer*4 function openfile(outfile)
implicit none
character(*), intent(in) :: outfile
integer*4 :: fp=FPOUT, info
if(len_trim(outfile)/=0)then
    open(unit=fp,file=outfile,status='old',iostat=info,position='append')
    if(info/=0)then
        openfile=6; return
    else
        openfile=fp; return
    endif
else
    openfile=6; return
endif
end function

! open output file for additional -------------------------------------------
integer*4 function addfile(outfile)
implicit none
character(*), intent(in) :: outfile
integer*4 :: fp=FPADD, info
if(len_trim(outfile)/=0)then
    open(unit=fp,file=outfile,status='old',iostat=info,position='append')
    if(info/=0)then
        addfile=6; return
    else
        addfile=fp; return
    endif
else
    addfile=6; return
endif
end function

! 0-error, 1-right
! execute processing session ------------------------------------------------
integer*4 function execses(ts, te, ti, popt, sopt, flag, infile, myindex, n, outfile)  ! output additional files
implicit none
type(gtime_t), intent(in) :: ts, te
real*8, intent(in) :: ti
type(prcopt_t), intent(in) :: popt
type(solopt_t), intent(in) :: sopt
integer*4, intent(in) :: flag, myindex(:), n
character(*), intent(in) :: infile(:), outfile
integer*4 fp, info
type(rtk_t) rtk
type(prcopt_t) popt_
type(gtime_t) stamp_
popt_=popt

! read obs and nav data 
!stamp_=timeget(); write(*,*) "Before readobsnav : ", stamp_%time*1d0+stamp_%sec
call readobsnav(ts,te,ti,infile,myindex,n,popt_,obss,navs,stas,info)
!stamp_=timeget(); write(*,*) "After readobsnav  : " , stamp_%time*1d0+stamp_%sec

if (info==0)then
    execses=0; return
endif
! write header to output file 
if(flag/=0 .and. sopt%issingle==0 .and. headwritten_==0)then
    if(outhead(outfile,infile,n,popt_,sopt,stas(1),obss)==0)then
        call freeobsnav(obss,navs)
        execses=0; return
    endif
    headwritten_=1
endif
iobsu=0
if(popt_%mode==PMODE_SINGLE .or. popt_%soltype==0)then
    fp=openfile(outfile)
    call procpos(fp,popt_,sopt,rtk,0,info)
    if(fp/=6)then
        if(sopt%issingle==0) close(unit=fp,status='keep')
        !if(sopt%issingle==1) close(unit=fp,status='delete')
    endif
    if(info==0)then
        call freeobsnav(obss,navs)
        execses=0; return
    endif
endif
!stamp_=timeget(); write(*,*) "After procpos :     " , stamp_%time*1d0+stamp_%sec

! free obs and nav data 
call freeobsnav(obss,navs)
execses=1
end function

! execute processing session for each rover ---------------------------------
integer*4 function execses_r(ts, te, ti, popt, sopt, flag, infile, myindex, n, outfile, rov)
implicit none
type(gtime_t), intent(in) :: ts, te
real*8, intent(in) :: ti
type(prcopt_t), intent(in) :: popt
type(solopt_t), intent(in) :: sopt
integer*4, intent(in) :: flag, myindex(:), n
character(*), intent(in) :: infile(:), outfile, rov
integer*4 :: stat
stat=0
! execute processing session 
stat=execses(ts,te,ti,popt,sopt,flag,infile,myindex,n,outfile)
execses_r=stat
end function

! execute processing session for each base station --------------------------
integer*4 function execses_b(ts, te, ti, popt, sopt, flag, infile, myindex, n, outfile, rov, base)
implicit none
type(gtime_t), intent(in) :: ts, te
real*8, intent(in) :: ti
type(prcopt_t), intent(in) :: popt
type(solopt_t), intent(in) :: sopt
integer*4, intent(in) :: flag, myindex(:), n
character(*), intent(in) :: infile(:), outfile, rov, base
integer*4 :: stat
stat=0
stat=execses_r(ts,te,ti,popt,sopt,flag,infile,myindex,n,outfile,rov)
execses_b=stat
end function

! post-processing positioning -------------------------------------------------
! post-processing positioning
! args   : gtime_t ts       I   processing start time (ts.time==0: no limit)
!        : gtime_t te       I   processing end time   (te.time==0: no limit)
!          real*8 ti        I   processing interval  (s) (0:all)
!          real*8 tu        I   processing unit time (s) (0:all)
!          prcopt_t *popt   I   processing options
!          solopt_t *sopt   I   solution options
!          filopt_t *fopt   I   file options
!          char   **infile  I   input files (see below)
!          integer*4    n   I   number of input files
!          char   *outfile  I   output file ('':stdout, see below)
!          char   *rov      I   rover id list        (separated by ' ')
!          char   *base     I   base station id list (separated by ' ')
! return : status (0:ok,0>:error,1:aborted)
! notes  : input files should contain observation data, navigation data, precise 
!          ephemeris/clock (optional), sbas log file (optional), ssr message
!          log file (optional) and tec grid file (optional). only the first 
!          observation data file in the input files is recognized as the rover
!          data.
!
!          the type of an input file is recognized by the file extention as )
!          follows:
!              .sp3,.SP3,.eph*,.EPH*: precise ephemeris (sp3c)
!              .sbs,.SBS,.ems,.EMS  : sbas message log files (rtklib or ems)
!              .lex,.LEX            : qzss lex message log files
!              .rtcm3,.RTCM3        : ssr message log files (rtcm3)
!              .*i,.*I              : tec grid files (ionex)
!              .fcb,.FCB            : satellite fcb
!              others               : rinex obs, nav, gnav, hnav, qnav or clock
!
!          inputs files can include wild-cards (*). if an file includes
!          wild-cards, the wild-card expanded multiple files are used.
!
!          inputs files can include keywords. if an file includes keywords,
!          the keywords are replaced by date, time, rover id and base station
!          id and multiple session analyses run. refer reppath() for the
!          keywords.
!
!          the output file can also include keywords. if the output file does
!          not include keywords. the results of all multiple session analyses
!          are output to a single output file.
!
!          ssr corrections are valid only for forward estimation.
!-----------------------------------------------------------------------------
integer*4 function postpos(ts, te, ti, tu, popt, sopt, infile, n, outfile, rov, base)
implicit none
type(gtime_t), intent(in) :: ts, te
real*8, intent(in) :: ti, tu
type(prcopt_t), intent(in) :: popt
type(solopt_t), intent(in) :: sopt
character(*), intent(in) :: infile(:), outfile, rov, base
integer*4, intent(in) :: n
integer*4 :: i,stat,myindex(MAXINFILE)
stat=0; myindex=0
do i=1,n
    myindex(i)=i-1
enddo
! execute processing session 
stat=execses_b(ts,te,ti,popt,sopt,1,infile,myindex,n,outfile,rov,base)
postpos=stat
end function
end module postpos_f90_
