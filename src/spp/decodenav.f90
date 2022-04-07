
! Decode navigation data ----------------------------------------------------
! decode ephemeris ----------------------------------------------------------
subroutine decode_eph(ver, sat, toc, mydata, eph, stat)
implicit none
include 'file_para.h'
real*8, intent(in) :: ver, mydata(64)
integer*4, intent(in) :: sat
type(gtime_t), intent(in) :: toc
type(eph_t), intent(out) :: eph
integer*4, intent(out) :: stat  ! 0-error, 1-normal
integer*4 sys, prn
integer*4, external :: uraindex,sisa_index
real*8, external :: SQR,rcond
type(gtime_t), external :: gpst2time,adjweek

call satsys(sat,prn,sys)
if (and(sys,or(SYS_GPS,or(SYS_GAL,or(SYS_QZS,or(SYS_CMP,SYS_IRN)))))==0)then
    stat=0; return
endif
call init_eph(eph)

eph%sat=sat
eph%toc=toc
eph%f0=mydata(1)
eph%f1=mydata(2)
eph%f2=mydata(3)

eph%A=SQR(mydata(11)); eph%e=mydata( 9); eph%i0  =mydata(16); eph%OMG0=mydata(14)
eph%omg =mydata(18); eph%M0 =mydata( 7); eph%deln=mydata( 6); eph%OMGd=mydata(19)
eph%idot=mydata(20); eph%crc=mydata(17); eph%crs =mydata( 5); eph%cuc =mydata( 8)
eph%cus =mydata(10); eph%cic=mydata(13); eph%cis =mydata(15)
if(sys==SYS_GPS .or. sys==SYS_QZS)then
    eph%iode=int(mydata( 4))      ! IODE 
    eph%iodc=int(mydata(27))      ! IODC 
    eph%toes=   (mydata(12))      ! toe (s) in gps week 
    eph%week=int(mydata(22))      ! gps week 
    eph%toe=adjweek(gpst2time(eph%week,mydata(12)),toc)
    eph%ttr=adjweek(gpst2time(eph%week,mydata(28)),toc)
    
    eph%code=int(mydata(21))      ! GPS: codes on L2 ch 
    eph%svh =int(mydata(25))      ! sv health 
    eph%sva =uraindex(mydata(24)) ! ura (m%index) 
    eph%flag=int(mydata(23))      ! GPS: L2 P data flag 
    
    eph%tgd(1)=   mydata(26)      ! TGD 
    if(sys==SYS_GPS)then
        eph%fit=mydata(29)        ! fit interval (h) 
    else
        eph%fit=rcond(dabs(mydata(29))<=1d-20,1.d0,2.d0) ! fit interval (0:1h,1:>2h) 
    endif
elseif(sys==SYS_GAL)then
    eph%iode=int(mydata( 4))      ! IODnav 
    eph%toes=   (mydata(12))      ! toe (s) in galileo week 
    eph%week=int(mydata(22))      ! gal week = gps week 
    eph%toe=adjweek(gpst2time(eph%week,mydata(12)),toc)
    eph%ttr=adjweek(gpst2time(eph%week,mydata(28)),toc)
    
    eph%code=int(mydata(21))      ! data sources 
                                  ! bit 0 set: I/NAV E1-B 
                                  ! bit 1 set: F/NAV E5a-I 
                                  ! bit 2 set: F/NAV E5b-I 
                                  ! bit 8 set: af0-af2 toc are for E5a.E1 
                                  ! bit 9 set: af0-af2 toc are for E5b.E1 
    eph%svh =int(mydata(25))      ! sv health 
                                  ! bit     0: E1B DVS 
                                  ! bit   1-2: E1B HS 
                                  ! bit     3: E5a DVS 
                                  ! bit   4-5: E5a HS 
                                  ! bit     6: E5b DVS 
                                  ! bit   7-8: E5b HS 
    eph%sva =sisa_index(mydata(24)) ! sisa (m->index) 
    eph%tgd(1)=mydata(26)         ! BGD E5a/E1 
    eph%tgd(2)=mydata(27)         ! BGD E5b/E1 
elseif(sys==SYS_CMP)then
elseif(sys==SYS_IRN)then
endif
stat=1
end subroutine

! decode nav header ---------------------------------------------------------
subroutine decode_navh(buff, nav)
implicit none
include 'file_para.h'
character(*), intent(in) :: buff
type(nav_t), intent(out) :: nav
integer*4 i,j,info
character(20) :: label
label=buff(61:)

if(index(label,'ION ALPHA')/=0)then
    read(buff(3:50),"(4D12.4)",iostat=info) nav%ion_gps(1:4)
elseif(index(label,'ION BETA')/=0)then
    read(buff(3:50),"(4D12.4)",iostat=info) nav%ion_gps(5:8)
elseif(index(label,'DELTA-UTC: A0,A1,T,W')/=0)then  ! opt ver.2 
    read(buff(4:59),"(2D19.12,2D9.0)",iostat=info) nav%utc_gps(1:4)
elseif(index(label,'IONOSPHERIC CORR')/=0)then
    if(buff(1:4)=='GPSA') read(buff(6:53),"(4D12.4)",iostat=info) nav%ion_gps(1:4)
    if(buff(1:4)=='GPSB') read(buff(6:53),"(4D12.4)",iostat=info) nav%ion_gps(5:8)
    if(buff(1:4)=='GAL ') read(buff(6:53),"(4D12.4)",iostat=info) nav%ion_gal(1:4)
    if(buff(1:4)=='BDSA') read(buff(6:53),"(4D12.4)",iostat=info) nav%ion_cmp(1:4)
    if(buff(1:4)=='BDSB') read(buff(6:53),"(4D12.4)",iostat=info) nav%ion_cmp(5:8)
elseif(index(label,'TIME SYSTEM CORR')/=0)then      ! opt ver.3 
    !if(buff(1:4)=='GPUT') read(buff(5:50),"(D18.10,D16.9,I7,I5)",iostat=info) nav%utc_gps(1:4)
    !if(buff(1:4)=='GLUT') read(buff(5:50),"(D18.10,D16.9,I7,I5)",iostat=info) nav%utc_glo(1:4)
    if(buff(1:4)=='GPUT') read(buff(5:50),*,iostat=info) nav%utc_gps(1:4)
    if(buff(1:4)=='GLUT') read(buff(5:50),*,iostat=info) nav%utc_glo(1:4)
    if(buff(1:4)=='GAUT') read(buff(5:50),*,iostat=info) nav%utc_gal(1:4)  ! v.3.02 
    if(buff(1:4)=='QZUT') read(buff(5:50),*,iostat=info) nav%utc_qzs(1:4)
    if(buff(1:4)=='BDUT') read(buff(5:50),*,iostat=info) nav%utc_cmp(1:4)
    if(buff(1:4)=='SBUT') read(buff(5:50),*,iostat=info) nav%utc_sbs(1:4)
    if(buff(1:4)=='IRUT') read(buff(5:50),*,iostat=info) nav%utc_irn(1:4)
elseif(index(label,'LEAP SECONDS')/=0)then
    read(buff(1:6),"(I6)",iostat=info) nav%leaps
endif
end subroutine

subroutine decode_geph(ver, sat, toc_in, mydata, geph, stat)
implicit none
include 'file_para.h'
real*8, intent(in) :: ver, mydata(64)
integer*4, intent(in) :: sat
type(gtime_t), intent(in) :: toc_in
type(geph_t), intent(out) :: geph
integer*4, intent(out) :: stat  ! 0-error, 1-normal
type(gtime_t) toc, tof
real*8 tow,tod
integer*4 week,dow,prn,sys
real*8, external :: rcond
type(gtime_t), external :: gpst2time,adjday,utc2gpst
toc=toc_in
call satsys(sat,prn,sys)
if(sys/=SYS_GLO)then
    stat=0; return
endif
call init_geph(geph)
geph%sat=sat

! toc rounded by 15 min in utc
call time2gpst(toc,week,tow)
toc=gpst2time(week,floor((tow+450.d0)/900.d0)*900.d0)
dow=int(floor(tow/86400.d0))

! time of frame in utc 
tod=rcond(ver<=2.99,mydata(3),dmod(mydata(3),86400.d0))  ! tod (v.2), tow (v.3) in utc 
tof=gpst2time(week,tod+dow*86400.d0)
tof=adjday(tof,toc)
geph%toe=utc2gpst(toc)    ! toc (gpst)
geph%tof=utc2gpst(tof)    ! tof (gpst)

! iode = tb (7bit), tb =index of UTC+3H within current day 
geph%iode=int(dmod(tow+10800.d0,86400.d0)/900.d0+0.5d0)
geph%taun=-mydata(1)       ! -taun 
geph%gamn= mydata(2)       ! +gamman 

geph%pos(1)=mydata(4)*1d3; geph%pos(2)=mydata(8)*1d3;  geph%pos(3)=mydata(12)*1d3
geph%vel(1)=mydata(5)*1d3; geph%vel(2)=mydata(9)*1d3;  geph%vel(3)=mydata(13)*1d3
geph%acc(1)=mydata(6)*1d3; geph%acc(2)=mydata(10)*1d3; geph%acc(3)=mydata(14)*1d3

geph%svh=int(mydata( 7))
geph%frq=int(mydata(11))
geph%age=int(mydata(15))

! some receiver output >128 for minus frequency number
if (geph%frq>50 .or. geph%frq<-50)then
    stat=0; return
endif
!if (geph%frq>128) geph%frq=geph%frq-256
!if (geph%frq<MINFREQ_GLO .or. MAXFREQ_GLO<geph%frq) {}
stat=1
end subroutine

! decode gnav header --------------------------------------------------------
subroutine decode_gnavh(buff, nav)
implicit none
include 'file_para.h'
character(*), intent(in) :: buff
type(nav_t), intent(out) :: nav
character(20) :: label
integer*4 info
label=buff(61:)
!if(index(label,'CORR TO SYTEM TIME')/=0) ;
if(index(label,'LEAP SECONDS')/=0)then
    read(buff(1:6),"(I6)",iostat=info) nav%leaps
endif
end subroutine
