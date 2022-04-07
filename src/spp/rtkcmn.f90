
! Common function -----------------------------------------------------------
! conditional operation (integer) -------------------------------------------
integer*4 function icond(key,a,b)
implicit none
include 'file_para.h'
logical*4, intent(in) :: key
integer*4, intent(in) :: a,b
if(key)then
    icond=a; return
else
    icond=b; return
endif
end function

! conditional operation (real) ----------------------------------------------
real*8 function rcond(key,a,b)
implicit none
include 'file_para.h'
logical*4, intent(in) :: key
real*8, intent(in) :: a,b
if(key)then
    rcond=a; return
else
    rcond=b; return
endif
end function

! screen by time ------------------------------------------------------------
! screening by time start, time end, and time interval
! args   : gtime_t time  I      time
!          gtime_t ts    I      time start (ts.time==0:no screening by ts)
!          gtime_t te    I      time end   (te.time==0:no screening by te)
!          double  tint  I      time interval (s) (0.0:no screen by tint)
! return : 1:on condition, 0:not on condition
! ---------------------------------------------------------------------------
integer*4 function screent(time, ts, te, tint)
implicit none
include 'file_para.h'
type(gtime_t), intent(in) :: time, ts, te
real*8, intent(in) :: tint
integer*4 week
real*8 gpst, dwnd, timediff
external :: timediff
logical(1) l1, l2, l3
call time2gpst(time,week,gpst)
dwnd=min(tint/10.d0,0.3d0)
! l1=(tint<=0.d0 .or. dmod(gpst+DTTOL,tint)<=DTTOL*2.d0)
! l1=(tint<=0.d0 .or. dmod(timediff(time,ts),tint)>=tint-DTTOL .or. dmod(timediff(time,ts),tint)<=DTTOL)
l1=(tint<=0.d0 .or. dmod(timediff(time,ts),tint)>=tint-dwnd .or. dmod(timediff(time,ts),tint)<=dwnd)
l2=(ts%time==0 .or. timediff(time,ts)>=-dwnd)
l3=(te%time==0 .or. timediff(time,te)<= dwnd)  ! < dwnd
if(l1 .and. l2 .and. l3)then
    screent=1
else
    screent=0
endif
end function

! compute dops --------------------------------------------------------------
! * * compute DOP (dilution of precision)
! * * args   : integer*4    ns        I   number of satellites
! * *          real*8 *azel     I   satellite azimuth/elevation angle (rad)
! * *          real*8 elmin     I   elevation cutoff angle (rad)
! * *          real*8 *dop      O   DOPs {GDOP,PDOP,HDOP,VDOP}
! * * return : none
! * * notes  : dop(0)-(3) return 0 in case of dop computation error
! * *------------------------------------------------------------------------
subroutine dops(ns, azel, elmin, dop)
implicit none
include 'file_para.h'
integer*4, intent(in) :: ns
real*8, intent(in) :: azel(ns,2),elmin
real*8, intent(out) :: dop(4)
real*8 H(MAXSAT,4),Q(4,4),cosel,sinel,SQRT2
integer*4 :: i,n,info
external :: SQRT2
i=1; n=1
dop=0.d0
do i=1,min(ns,MAXSAT)
    if (azel(i,2)<elmin .or. azel(i,2)<=0.0) cycle
    cosel=dcos(azel(i,2))
    sinel=dsin(azel(i,2))
    H(n,1)=cosel*dsin(azel(i,1))
    H(n,2)=cosel*dcos(azel(i,1))
    H(n,3)=sinel
    H(n,4)=1.d0
    n=n+1
enddo
if (n-1<4) return

call matmul2("TN",4,4,n-1,H(1:n-1,:),H(1:n-1,:),Q)
call matinv(Q,4,info)
if (info==0)then
    dop(1)=SQRT2(Q(1,1)+Q(2,2)+Q(3,3)+Q(4,4))  ! GDOP 
    dop(2)=SQRT2(Q(1,1)+Q(2,2)+Q(3,3))         ! PDOP 
    dop(3)=SQRT2(Q(1,1)+Q(2,2))                ! HDOP 
    dop(4)=SQRT2(Q(3,3))                       ! VDOP 
endif
end subroutine

! satellite azimuth/elevation angle -----------------------------------------
! compute satellite azimuth/elevation angle
! args   : real*8 *pos      I   geodetic position {lat,lon,h} (rad,m)
!          real*8 *e        I   receiver-to-satellilte unit vevtor (ecef)
!          real*8 *azel     IO  azimuth/elevation {az,el} (rad) (NULL: no output)
!                               (0.d0<=azel(0)<2*pi,-pi/2<=azel(1)<=pi/2)
! return : elevation angle (rad)
!----------------------------------------------------------------------------
subroutine satazel(pos, e, azel, el)
implicit none
include 'file_para.h'
real*8, intent(in) :: pos(3),e(3,1)
real*8, intent(out) :: azel(2),el
real*8 :: az, enu(3), dot, rcond
external :: dot, rcond
az=0; el=pi/2.d0
if(pos(3)>-RE_WGS84)then
    call ecef2enu(pos,e,enu)
    az=rcond(dot(enu,enu,2)<1d-12,0.d0,datan2(enu(1),enu(2)))
    if (az<0.d0) az=az+2.d0*PI
    el=dasin(enu(3))
endif
azel(1)=az; azel(2)=el
end subroutine

! save slips ----------------------------------------------------------------
subroutine saveslips(slips, mydata)
implicit none
include 'file_para.h'
integer*4, intent(out) :: slips(MAXSAT,NFREQ)  !slips(:,:)  !(:,NFREQ)
type(obsd_t), intent(in):: mydata
integer*4 i
do i=1, NFREQ
    if(and(mydata%LLI(i),1)/=0)then
        slips(mydata%sat,i)=or(slips(mydata%sat,i),LLI_SLIP)
    endif
enddo
end subroutine

! restore slips -------------------------------------------------------------
subroutine restslips(slips, mydata)
implicit none
include 'file_para.h'
integer*4, intent(out) :: slips(MAXSAT,NFREQ)  !slips(:,:)
type(obsd_t), intent(out) :: mydata
integer*4 i
do i=1, NFREQ
    if (and(slips(mydata%sat,i),1)/=0) mydata%LLI(i)=or(mydata%LLI(i),LLI_SLIP)
    slips(mydata%sat,i)=0
enddo
end subroutine

! set signal index ----------------------------------------------------------
subroutine set_index(ver, sys, opt, tobs, ind)
implicit none
include 'file_para.h'
real*8, intent(in) :: ver
integer*4, intent(in) :: sys
character(*), intent(in) :: opt
character(3), intent(in) :: tobs(MAXOBSTYPE)
type(sigind_t), intent(out) :: ind
character(1) p
real*8 shift
integer*4 :: i,j,k,n,info,icond,getcodepri
logical*1 ltmp
external :: icond,getcodepri
ltmp=.false.
i=1; n=0
do while(tobs(i)/='end')
    if(tobs(i)=='')then  !? to be verified
        ind%frq(i)=0; ind%code(i)=0
        ind%mytype(i)=0;
        if(sys==SYS_GPS)then
            ind%pri(i)=6
        elseif(sys==SYS_SBS.or.sys==SYS_GLO.or.sys==SYS_GAL.or.sys==SYS_QZS.or.sys==SYS_CMP.or.sys==SYS_IRN)then
            ind%pri(i)=14
        elseif(sys==SYS_ALL.or.sys==SYS_LEO.or.sys==SYS_NONE)then
            ind%pri(i)=0
        endif
        ind%pos(i)=-1
        i=i+1; if(i>size(tobs,1)) exit
        n=n+1; cycle
    endif
    call obs2code(tobs(i)(2:3),ind%frq(i),ind%code(i))
    info=index(obstyps,tobs(i)(1:1))
    ind%mytype(i)=icond(info/=0,info-1,-1)
    ind%pri(i)=getcodepri(sys,ind%code(i),opt)
    ind%pos(i)=-1
    if (sys==SYS_CMP)then
        if (ind%frq(i)==5) ind%frq(i)=2; ! B2 
        if (ind%frq(i)==4) ind%frq(i)=3; ! B3 
    endif
    i=i+1; n=n+1
    if(i>size(tobs,1)) exit
enddo
! assign index for highest priority code 
do i=1,NFREQ
    k=0
    do j=1,n
        if(k<1)then
            ltmp=.true.
        elseif(ind%pri(j)>ind%pri(k))then
            ltmp=.true.
        else
            ltmp=.false.
        endif
        if(ind%frq(j)==i.and.ind%pri(j)/=0.and.ltmp) k=j
    enddo
    if(k<1) cycle
    do j=1,n
        if(ind%code(j)==ind%code(k)) ind%pos(j)=i-1
        if(ind%frq(j)==ind%frq(k))   ind%pos(j)=i-1  !add
    enddo
enddo
!! assign index of extended obs data 
!do i=1,NEXOBS
!    do j=1,n
!        if(ind%code(j)/=0 .and. ind%pri(j)/=0 .and. ind%pos(j)<0) exit
!    enddo
!    if(j>n) exit
!    do k=1,n
!        if(ind%code(k)==ind%code(j)) ind%pos(k)=i-1  !NFREQ+i-1
!    enddo
!enddo
do i=1,n
    if (ind%code(i)==0 .or. ind%pri(i)==0 .or. ind%pos(i)>=0) cycle
enddo
ind%n=n
end subroutine

! galileo sisa value (m) to sisa index --------------------------------------
integer*4 function sisa_index(value)
implicit none
include 'file_para.h'
real*8, intent(in) :: value
if (value<0.d0  .or.  value>6.d0)then  ! unknown or NAPA 
    sisa_index=255
elseif (value<=0.5d0)then
    sisa_index=int(value/0.01d0)
elseif (value<=1.d0)then
    sisa_index=int((value-0.5d0)/0.02d0)+50
elseif (value<=2.d0)then
    sisa_index=int((value-1.d0)/0.04d0)+75
else
    sisa_index=int(int(value-2.d0)/0.16d0)+100
endif
end function

! ura value (m) to ura index ------------------------------------------------
integer*4 function uraindex(value)
implicit none
include 'file_para.h'
real*8, intent(in) :: value
integer*4 i
do i=1,15
    if(ura_eph(i)>=value) exit
enddo
uraindex=i-1
end function

!! variance by ura ephemeris (ref (1) 20.3.3.3.1.1) -------------------------
real*8 function var_uraeph(sys, ura)
implicit none
include 'file_para.h'
integer*4, intent(in) :: sys, ura
real*8, external :: SQR,rcond
if(sys==SYS_GAL)then    ! galileo sisa (ref (7) 5.1.11) 
    if (ura<= 49) var_uraeph=SQR(ura*0.01d0)
    if (ura<= 74) var_uraeph=SQR(0.5+(ura-50)*0.02d0)
    if (ura<= 99) var_uraeph=SQR(1.d0+(ura-75)*0.04d0)
    if (ura<=125) var_uraeph=SQR(2.d0+(ura-100)*0.16d0)
    if (ura> 125) var_uraeph=SQR(STD_GAL_NAPA)
else    ! gps ura (ref (1) 20.3.3.3.1.1) 
    var_uraeph=rcond(ura<0.or.14<ura,SQR(6144.d0),SQR(ura_value(ura+1)))
endif
end function

! pseudorange measurement error variance ------------------------------------
real*8 function varerr(opt, el, sys)
implicit none
include 'file_para.h'
type(prcopt_t), intent(in) :: opt
real*8, intent(in) :: el
integer*4, intent(in) :: sys
real*8 fact,varr,rcond,SQR
external :: rcond,SQR
fact=rcond(sys==SYS_GLO,EFACT_GLO,EFACT_GPS)
varr=SQR(100.d0)*(SQR(0.003d0)+SQR(0.003d0)/dsin(el))
if (opt%ionoopt==IONOOPT_IFLC) varr=varr*SQR(3.d0); ! iono-free 
varerr=SQR(fact)*varr;
end function

! set system mask -----------------------------------------------------------
integer*4 function set_sysmask(opt)  ! -SYS=GR
implicit none
include 'file_para.h'
character(*), intent(in) :: opt
character(1) p
integer*4 :: mask, i
mask=SYS_NONE
i=index(opt,'-SYS=')
if (i==0)then
    set_sysmask=SYS_ALL; return
endif
i=i+5
do while(i<len_trim(opt))
    p=opt(i:i)
    select case(p)
    case('G')
        mask=or(mask,SYS_GPS)
    case('R')
        mask=or(mask,SYS_GLO)
    case('E')
        mask=or(mask,SYS_GAL)
    case('J')
        mask=or(mask,SYS_QZS)
    case('C')
        mask=or(mask,SYS_CMP)
    case('I')
        mask=or(mask,SYS_IRN)
    case('S')
        mask=or(mask,SYS_SBS)
    end select
    i=i+1
    if(p==' ') exit
enddo
set_sysmask=mask
end function
