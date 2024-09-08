
! Read obs and nav data -----------------------------------------------------
subroutine readobsnav(ts, te, ti, infile, myindex, n, prcopt, obs, nav, sta, stat)
implicit none
include 'file_para.h'
type(gtime_t), intent(in) :: ts, te
real*8, intent(in) :: ti
character(*), intent(in) :: infile(n)
integer*4, intent(in) :: myindex(MAXINFILE), n
type(prcopt_t), intent(in) :: prcopt
type(obs_t), intent(out) :: obs
type(nav_t), intent(out) :: nav
type(sta_t), intent(out) :: sta(MAXRCV)
integer*4, intent(out) :: stat  ! 0-error, 1-normal
integer*4 :: i,j,ind,nobs,rcv,info
integer*4, external :: icond
type(sta_t) statmp
nobs=0; rcv=1
nullify(obs%mydata); obs%n =0; obs%nmax =0
nepoch=0

! read rinex obs and nav file 
call readrnxt(infile(1),1,ts,te,ti,'',obs,nav,site,info)  ! -1,0: error, 1-normal
if(info<=0) write(*,*) "File reading error : ", trim(infile(1))

if (info<0)then
    stat=0; return
endif
if(obs%n<=0)then
    write(*,*) "Error : no observation data!"
    stat=0; return
endif

if(nav%n<=0 .and. nav%ng<=0 .and. nav%ns<=0)then
    write(*,*) "Error : no navigation data!"
    stat=0; return
endif
! sort observation data 
nepoch=obs%n

! delete duplicated ephemeris 
call uniqnav(nav)
stat=1
end subroutine

! read rinex obs and nav files ----------------------------------------------
subroutine readrnxt(filepath, rcv, ts, te, tint, opt, obs, nav, sta, stat)
implicit none
include 'file_para.h'
character(*), intent(in) :: filepath, opt
integer*4, intent(in) :: rcv
type(gtime_t), intent(in) :: ts, te
real*8, intent(in) :: tint
type(obs_t), intent(out) :: obs
type(nav_t), intent(out) :: nav
type(sta_t), intent(out) :: sta
integer*4, intent(out) :: stat  ! -1,0: error, 1-normal
integer*4 :: i,n,stat1
character(1) :: mytype
character(150) :: files(2)
stat1=0
mytype=''; files=''
files(1)=filepath
! read rinex files 
write(*,*) files(1),ts,te,tint,opt,0,rcv,mytype
call readrnxfile(files(1),ts,te,tint,'',0,1,mytype,obss,navs,site,stat1)  ! -1,0: error, 1-normal

stat=stat1
end subroutine

! uncompress and read rinex file --------------------------------------------
subroutine readrnxfile(filepath, ts, te, tint, opt, flag, index, mytype, obs, nav, sta, stat)
implicit none
include 'file_para.h'
character(*), intent(in) :: filepath, opt
type(gtime_t), intent(in) :: ts, te
real*8, intent(in) :: tint
integer*4, intent(in) :: flag, index
character(*), intent(out) :: mytype
type(obs_t), intent(out) :: obs
type(nav_t), intent(out) :: nav
type(sta_t), intent(out) :: sta
integer*4, intent(out) :: stat  ! -1,0: error, 1-normal
integer*4 :: fp=FPNAV, stat1, info

call init_sta(sta)
open(unit=fp,file=filepath,status='old',iostat=info)
if(info/=0)then
    write(*,*) "File opening error : ", trim(filepath)
    stat=-1; return
endif
! read rinex file 
call readrnxfp(fp,ts,te,tint,opt,flag,index,mytype,obs,nav,sta,stat1)  ! 0-error, 1-normal
if(sta%name=='')then
    call decodemarker(filepath, sta%name)
endif
close(unit=fp,status='keep')
stat=stat1
end subroutine

! read rinex file -----------------------------------------------------------
subroutine readrnxfp(fp, ts, te, tint, opt, flag, index, mytype, obs, nav, sta, stat)
implicit none
include 'file_para.h'
integer*4, intent(in) :: fp, flag, index
type(gtime_t), intent(in) :: ts, te
real*8, intent(in) :: tint
character(*), intent(in) :: opt
character(*), intent(in) :: mytype
type(obs_t), intent(out) :: obs
type(nav_t), intent(out) :: nav
type(sta_t), intent(out) :: sta
integer*4, intent(out) :: stat  ! 0-error, 1-normal
real*8 ver
integer*4 :: sys,tsys,info
character(3) :: tobs(NUMSYS,MAXOBSTYPE)

tsys=TSYS_GPS
tobs=''
! read rinex header 
call readrnxh(fp,ver,mytype,sys,tsys,tobs,obs,nav,sta,info)
if (info==0)then
    write(*,*) "Error reading rinex header ..."
    stat=0; return
endif
! read rinex body 

call readrnxobs(fp,ts,te,tint,opt,index,ver,tsys,tobs,obs,sta,stat); return  ! 0-error, 1-normal

stat=0
end subroutine

! free obs and nav data -----------------------------------------------------
subroutine freeobsnav(obs, nav)
implicit none
include 'file_para.h'
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

! free nav data ----------------------------------------------------------
subroutine freenav(nav)
implicit none
include 'file_para.h'
type(nav_t), intent(out) :: nav
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