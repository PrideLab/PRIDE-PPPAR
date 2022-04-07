
! Time function -------------------------------------------------------------
! adjust time considering week handover -------------------------------------
type(gtime_t) function adjday(t, t0)
implicit none
include 'file_para.h'
type(gtime_t), intent(in) :: t,t0
type(gtime_t) timeadd
real*8 tt, timediff
external :: timeadd, timediff
tt=timediff(t,t0)
if (tt<-43200.d0)then
    adjday=timeadd(t, 86400.d0); return
endif
if (tt> 43200.d0)then
    adjday=timeadd(t,-86400.d0); return
endif
adjday=t
end function

! adjust time considering week handover -------------------------------------
type(gtime_t) function adjweek(t, t0)
implicit none
include 'file_para.h'
type(gtime_t), intent(in) :: t,t0
type(gtime_t) timeadd
real*8 tt,timediff
external :: timeadd,timediff
tt=timediff(t,t0)
if (tt<-302400.d0)then
    adjweek=timeadd(t, 604800.d0); return
endif
if (tt> 302400.d0)then
    adjweek=timeadd(t,-604800.d0); return
endif
adjweek=t
end function

! convert calendar day/time to time -----------------------------------------
! convert calendar day/time to gtime_t struct
! args   : real*8 *ep       I   day/time {year,month,day,hour,min,sec}
! return : gtime_t struct
! notes  : proper in 1970-2037 or 1970-2099 (64bit time_t)
!----------------------------------------------------------------------------
type(gtime_t) function epoch2time(ep)
implicit none
include 'file_para.h'
real*8, intent(in) :: ep(6)
integer*4 :: doy(12)=(/1,32,60,91,121,152,182,213,244,274,305,335/)
type(gtime_t) :: mytime
integer*4 days,sec,year,mon,day,icond
external :: icond
mytime=gtime_t(0,0.d0)
year=int(ep(1)); mon=int(ep(2)); day=int(ep(3))
if (year<1970 .or. 2099<year .or. mon<1 .or. 12<mon)then
    epoch2time=mytime; return
endif
days=(year-1970)*365+(year-1969)/4+doy(mon)+day-2+icond(mod(year,4)==0.and.mon>=3,1,0)
sec=floor(ep(6))
mytime%time=days*86400+int(ep(4))*3600+int(ep(5))*60+sec
mytime%sec=ep(6)-sec
epoch2time=mytime
end function

! gps time to time ----------------------------------------------------------
! convert week and tow in gps time to gtime_t struct
! args   : integer*4 week   I   week number in gps time
!          real*8 sec       I   time of week in gps time (s)
! return : gtime_t struct
!----------------------------------------------------------------------------
type(gtime_t) function gpst2time(week, sec)
implicit none
include 'file_para.h'
integer*4, intent(in) :: week
real*8, intent(in) :: sec
type(gtime_t) t, epoch2time
real*8 sec2
external :: epoch2time
sec2=sec; t=epoch2time(gpst0)
if (sec2<-1d9 .or. 1d9<sec2) sec2=0.d0
t%time=t%time+86400*7*week+int(sec2)
t%sec=sec2-int(sec2)
gpst2time=t
end function

! gpstime to utc ------------------------------------------------------------
! convert gpstime to utc considering leap seconds
! args   : gtime_t t        I   time expressed in gpstime
! return : time expressed in utc
! notes  : ignore slight time offset under 100 ns
!----------------------------------------------------------------------------
type(gtime_t) function gpst2utc(t)
implicit none
include 'file_para.h'
type(gtime_t), intent(in) :: t
type(gtime_t) tu, timeadd, epoch2time
real*8 timediff
integer*4 i
external :: timeadd, epoch2time, timediff
tu=gtime_t(0,0.d0)
i=1
do while(.true.)
    if(leaps(i,1)<=0) exit
    tu=timeadd(t,leaps(i,7))
    if(timediff(tu,epoch2time(leaps(i,1:6)))>=0.d0)then
        gpst2utc=tu; return
    endif
    i=i+1
enddo
gpst2utc=t
end function

! string to time ------------------------------------------------------------
! convert substring in string to gtime_t struct
! args   : char   *s        I   string ('... yyyy mm dd hh mm ss ...')
!          integer*4 i,n    I   substring position and width
!          gtime_t *t       O   gtime_t struct
! return : status (0:ok,0>:error)
!----------------------------------------------------------------------------
subroutine str2time(s, i, n, t, stat)
implicit none
include 'file_para.h'
character(*), intent(in) :: s
integer*4, intent(in) :: i,n
type(gtime_t), intent(out) :: t
integer*4, intent(out) :: stat  ! -1: error, 0-normal
type(gtime_t) epoch2time
real*8 ep(6),rcond
integer*4 err
external :: rcond, epoch2time
t=gtime_t(0,0.d0); stat=0
if (i<1 .or. len_trim(s)<i+n-1 .or. n<1)then
    stat=-1; return
endif
read(s(i:i+n-1),*,iostat=err) ep(1:6)
if(err/=0)then
    stat=-1; return
endif
if (ep(1)<100.d0) ep(1)=ep(1)+rcond(ep(1)<80.d0,2000.d0,1900.d0)
t=epoch2time(ep)
end subroutine

! time to calendar day/time -------------------------------------------------
! convert gtime_t struct to calendar day/time
! args   : gtime_t t        I   gtime_t struct
!          real*8 *ep       O   day/time {year,month,day,hour,min,sec}
! return : none
! notes  : proper in 1970-2037 or 1970-2099 (64bit time_t)
!----------------------------------------------------------------------------
subroutine time2epoch(t, ep)
implicit none
include 'file_para.h'
type(gtime_t), intent(in) :: t
real*8, intent(out) :: ep(6)
integer*4 days,sec,mon,day
integer*4, parameter :: mday(48)=(/&
    31,28,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31,&
    31,29,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31/)
days=t%time/86400
sec=int(t%time-days*86400)
day=mod(days,1461)
do mon=1,48,1
    if (day>=mday(mon))then
        day=day-mday(mon)
    else
        exit
    endif
enddo
ep(1)=1970+days/1461*4+(mon-1)/12; ep(2)=mod(mon-1,12)+1; ep(3)=day+1
ep(4)=sec/3600; ep(5)=mod(sec,3600)/60; ep(6)=mod(sec,60)+t%sec
end subroutine

! time to gps time ----------------------------------------------------------
! convert gtime_t struct to week and tow in gps time
! args   : gtime_t t        I   gtime_t struct
!          integer*4 *week  IO  week number in gps time (NULL: no output)
! return : time of week in gps time (s)
!----------------------------------------------------------------------------
subroutine time2gpst(t, week, gpst)
implicit none
include 'file_para.h'
type(gtime_t), intent(in) :: t
integer*4, intent(out) :: week
real*8, intent(out) :: gpst
type(gtime_t) t0,epoch2time
integer*4 sec
external :: epoch2time
t0=epoch2time(gpst0); sec=t%time-t0%time
week=sec/(86400*7)
gpst=sec-week*86400*7.d0+t%sec
end subroutine

! convert time to mjd and sod -----------------------------------------------
subroutine time2mjd(t, mjd, sod)
implicit none
include 'file_para.h'
type(gtime_t), intent(in) :: t
integer*4, intent(out) :: mjd
real*8, intent(out) :: sod
integer*4 week,ymd2mjd
real*8 gpst,ep(6)
external :: ymd2mjd
call time2epoch(t,ep)
mjd=ymd2mjd(int(ep(1:3)))
call time2gpst(t,week,gpst)
sod=dmod(gpst,86400.d0)
end subroutine

! time to day and sec -------------------------------------------------------
subroutine time2sec(time, day, sec)
implicit none
include 'file_para.h'
type(gtime_t), intent(in) :: time
type(gtime_t), intent(out) :: day
real*8, intent(out) :: sec
type(gtime_t) epoch2time
external :: epoch2time
real*8 ep(6)
call time2epoch(time,ep)
sec=ep(4)*3600.d0+ep(5)*60.d0+ep(6)
ep(4:6)=0.d0
day=epoch2time(ep)
end subroutine

! time to string ------------------------------------------------------------
! convert gtime_t struct to string
! args   : gtime_t t        I   gtime_t struct
!          char   *s        O   string ('yyyy/mm/dd hh:mm:ss.ssss')
!          integer*4 n      I   number of decimals
! return : none
!----------------------------------------------------------------------------
subroutine time2str(t,s,n)
implicit none
include 'file_para.h'
type(gtime_t), intent(in) :: t
integer*4, intent(in) :: n
character(*), intent(out) :: s
real*8 ep(6)
type(gtime_t) t1
integer*4 n1,f1,f2,icond
character(60) fmtstr,f6
external :: icond
t1=t; n1=n
if(n1<0) n1=0; if(n1>12) n1=12
if(1.d0-t1%sec<0.5d0/(10.d0**n1))then
    t1%time=t1%time+1; t1%sec=0.d0
endif
call time2epoch(t1,ep)
write(f6,"(I2'.'I2)") icond(n1<=0,2,n1+3),icond(n1<=0,0,n1)
fmtstr="(I4.4'/'I2.2'/'I2.2' 'I2.2':'I2.2':'f"//trim(f6)//")"
write(s,fmtstr) int(ep(1:5)),ep(6)
end subroutine

! add time ------------------------------------------------------------------
! add time to gtime_t struct
! args   : gtime_t t        I   gtime_t struct
!          real*8 sec       I   time to add (s)
! return : gtime_t struct (t+sec)
!----------------------------------------------------------------------------
type(gtime_t) function timeadd(t, sec)
implicit none
include 'file_para.h'
type(gtime_t), intent(in) :: t
real*8, intent(in) :: sec
type(gtime_t) t1
real*8 tt
t1=t
t1%sec=t1%sec+sec; tt=floor(t1%sec); t1%time=t1%time+int(tt); t1%sec=t1%sec-tt
timeadd=t1
end function

! time difference -----------------------------------------------------------
! difference between gtime_t structs
! args   : gtime_t t1,t2    I   gtime_t structs
! return : time difference (t1-t2) (s)
!----------------------------------------------------------------------------
real*8 function timediff(t1,t2)
implicit none
include 'file_para.h'
type(gtime_t), intent(in) :: t1,t2
timediff=t1%time-t2%time+t1%sec-t2%sec
end function

! get current time in utc ---------------------------------------------------
! get current time in utc
! args   : none
! return : current time in utc
!----------------------------------------------------------------------------
type(gtime_t) function timeget()
implicit none
include 'file_para.h'
character(8) date
character(10) time_s
character(5) zone
integer*4 values(8)
real*8 :: ep(6)
type(gtime_t) time,epoch2time,timeadd
external :: epoch2time,timeadd
ep=0.d0
call date_and_time(date,time_s,zone,values)
ep(1:3)=values(1:3)  ! year, month, day
ep(4:6)=values(5:7)  ! h, m, s
ep(6)=ep(6)+values(8)*1d-3  ! ms
time=epoch2time(ep)
time=timeadd(time,-values(4)*60.d0)
time=timeadd(time,timeoffset_)
timeget=time
end function

! set current time in utc ---------------------------------------------------
! set current time in utc
! args   : gtime_t          I   current time in utc
! return : none
! notes  : just set time offset between cpu time and current time
!          the time offset is reflected to only timeget()
!          not reentrant
!----------------------------------------------------------------------------
subroutine timeset(t)
implicit none
include 'file_para.h'
type(gtime_t), intent(in) :: t
type(gtime_t), external :: timeget
real*8, external :: timediff
timeoffset_=timeoffset_+timediff(t,timeget())
end subroutine

! utc to gmst ---------------------------------------------------------------
! convert utc to gmst (Greenwich mean sidereal time)
! args   : gtime_t t        I   time expressed in utc
!          real*8 ut1_utc   I   UT1-UTC (s)
! return : gmst (rad)
!----------------------------------------------------------------------------
real*8 function utc2gmst(t, ut1_utc)
implicit none
include 'file_para.h'
type(gtime_t), intent(in) :: t
real*8, intent(in) :: ut1_utc
real*8 :: ep2000(6),timediff
type(gtime_t) tut,tut0,timeadd,epoch2time
real*8 ut,t1,t2,t3,gmst0,gmst
external :: timeadd,epoch2time,timediff
ep2000=(/2000,1,1,12,0,0/)
tut=timeadd(t,ut1_utc)
call time2sec(tut,tut0,ut)
t1=timediff(tut0,epoch2time(ep2000))/86400.d0/36525.d0
t2=t1*t1; t3=t2*t1
gmst0=24110.54841d0+8640184.812866d0*t1+0.093104d0*t2-6.2d-6*t3
gmst=gmst0+1.002737909350795d0*ut
utc2gmst=dmod(gmst,86400.d0)*PI/43200.d0
end function

! utc to gpstime ------------------------------------------------------------
! convert utc to gpstime considering leap seconds
! args   : gtime_t t        I   time expressed in utc
! return : time expressed in gpstime
! notes  : ignore slight time offset under 100 ns
!----------------------------------------------------------------------------
type(gtime_t) function utc2gpst(t)
implicit none
include 'file_para.h'
type(gtime_t), intent(in) :: t
type(gtime_t) epoch2time,timeadd
real*8 timediff,ep(6)
external :: epoch2time,timeadd,timediff
integer*4 i
i=1
do while(.true.)
    if(leaps(i,1)<=0) exit
    ep=leaps(i,1:6)
    if(timediff(t,epoch2time(ep))>=0.d0)then
        utc2gpst=timeadd(t,-leaps(i,7)); return
    endif
    i=i+1
end do
utc2gpst=t
end function
