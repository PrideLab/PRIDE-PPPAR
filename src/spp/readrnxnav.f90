! read rinex nav files ------------------------------------------------------
subroutine readrnxnavfiles(rnxnavs, nav, stat)
implicit none
include 'file_para.h'
character(*), intent(in) :: rnxnavs(FLNUMLIM)
type(nav_t), intent(out) :: nav
integer*4, intent(out) :: stat  ! -1,0: error, 1-normal

! local 
integer*4 :: i
character(1024) :: filetmp

! read rinex navs
nullify(nav%eph);    nav%n =0; nav%nmax =0
nullify(nav%geph);   nav%ng=0; nav%ngmax=0
do i=1,FLNUMLIM
    ! check if nav file is blank
    call getfname(rnxnavs(i),filetmp)
    if(len_trim(filetmp)==0) cycle

    ! read one rinex nav file
    call readrnxnavfile(rnxnavs(i),navs,stat)
    if(stat<=0) write(*,*) "File reading error : ", trim(rnxnavs(i))
enddo
end subroutine

! read rinex nav file -------------------------------------------------------
subroutine readrnxnavfile(filepath, nav, stat)
implicit none
include 'file_para.h'
character(*), intent(in) :: filepath
type(nav_t), intent(out) :: nav
integer*4, intent(out) :: stat  ! -1,0: error, 1-normal

! local
integer*4 :: fp=FPNAV, info, sys, tsys, i
real*8 :: ver
character(MAXRNXLEN) :: buff
character(20) :: label
character(1) :: mytype

! open rinex file 
open(unit=fp,file=filepath,status='old',iostat=info)
if(info/=0)then
    write(*,*) "File opening error : ", trim(filepath)
    stat=-1; return
endif

! read rinex nav header 
stat=0
do while(.true.)
    read(fp,"(A)",iostat=info) buff
    label=buff(61:)
    if(info/=0) exit
    if (len_trim(buff)<=60)then
        cycle
    elseif(index(label,'RINEX VERSION / TYPE')/=0)then
        read(buff(1:9),*,iostat=info) ver
        read(buff(21:40),"(A)",iostat=info) mytype
        ! satellite system 
        select case(buff(41:41))
        case(' ','G')
            sys=SYS_GPS; tsys=TSYS_GPS
        case('R')
            sys=SYS_GLO; tsys=TSYS_UTC
        case('E')
            sys=SYS_GAL; tsys=TSYS_GAL  ! v.2.12 
        case('S')
            sys=SYS_SBS; tsys=TSYS_GPS
        case('J')
            sys=SYS_QZS; tsys=TSYS_QZS  ! v.3.02 
        case('C')
            sys=SYS_CMP; tsys=TSYS_CMP  ! v.2.12 
        case('I')
            sys=SYS_IRN; tsys=TSYS_IRN  ! v.3.03 
        case('M')
            sys=SYS_NONE; tsys=TSYS_GPS  ! mixed 
        case default
        end select
        cycle
    elseif(index(label,'PGM / RUN BY / DATE')/=0)then
        cycle
    elseif(index(label,'COMMENT')/=0)then
        cycle
    endif
    ! file type 
    select case(mytype)
    case('N')
        call decode_navh(buff,nav)
    case('G')
        call decode_gnavh(buff,nav)
    end select
    if (index(label,'END OF HEADER')/=0)then
        stat=1; exit
    endif
    i=i+1
    if (i>=MAXPOSHEAD .and. mytype==' ') exit  ! no rinex file 
enddo
if (stat==0)then
    write(*,*) "Error reading rinex header ..."
    stat=0; return
endif

select case(mytype)
case('N')
    call readrnxnav(fp,'',ver,sys,nav,stat); call uniqnav(navs)
    return  ! 0-error, 1-normal
case('G')
    call readrnxnav(fp,'',ver,SYS_GLO,nav,stat); call uniqnav(navs)
    return
case('L')
    call readrnxnav(fp,'',ver,SYS_GAL,nav,stat); call uniqnav(navs)
    return  ! extension
end select

close(unit=fp,status='keep')
stat=0
end subroutine

! Read rinex nav/gnav/geo nav -----------------------------------------------
subroutine readrnxnav(fp, opt, ver, sys, nav, stat)
implicit none
include 'file_para.h'
integer*4, intent(in) :: fp, sys
character(*), intent(in) :: opt
real*8, intent(in) :: ver
type(nav_t), intent(out) :: nav
integer*4, intent(out) :: stat  ! 0-error, 1-normal

! local
type(eph_t) eph
type(geph_t) geph
integer*4 stat1, mytype

! function
real*8, external :: timediff
integer*4, external :: screent

! read rinex navigation data body 
do while(.true.)
    call readrnxnavb(fp,opt,ver,sys,mytype,eph,geph,stat1)  ! 0,-1: error, 1-normal
    if(stat1<0) exit  ! skip"-1"
    ! add ephemeris to navigation data 
    if(stat1/=0)then  ! stat=1
        select case(mytype)
        case(1)
            call add_geph(nav,geph,stat1)
        case default
            call add_eph (nav,eph, stat1)
        end select
    endif
enddo
if(nav%n>0 .or. nav%ng>0)then
    stat=1
else
    stat=0
endif
end subroutine

! read rinex navigation data body -------------------------------------------
subroutine readrnxnavb(fp, opt, ver, sys_in, mytype, eph, geph, stat)
implicit none
include 'file_para.h'
integer*4, intent(in) :: fp, sys_in
character(*), intent(in) :: opt
real*8, intent(in) :: ver
integer*4, intent(out) :: mytype, stat  ! 0,-1: error, 1-normal
type(eph_t), intent(out) :: eph
type(geph_t), intent(out) :: geph
type(gtime_t) toc
real*8 mydata(64)
integer*4 :: sys,i,j,prn,sat,sp,mask,info,prntmp
character(MAXRNXLEN) buff, p
character(8) :: id
integer*4, external :: set_sysmask,satid2no,satno
real*8, external :: str2num

i=1; sat=0; sp=3
id=''
! set system mask 
sys=sys_in
mask=set_sysmask(opt)
do while(.true.)
    read(fp,"(A)",iostat=info) buff
    if(info/=0) exit
    if(i==1)then
        ! decode satellite field 
        if(ver>=3.d0 .or. sys==SYS_GAL .or. sys==SYS_QZS)then  ! ver.3 or GAL/QZS 
            id=buff(1:3)
            sat=satid2no(id); sp=4
            if(ver>=3.d0) call satsys(sat,prntmp,sys)
        else
            prn=int(str2num(buff,1,2))
            if (sys==SYS_SBS) sat=satno(SYS_SBS,prn+100)
            if(sys==SYS_GLO)then
                sat=satno(SYS_GLO,prn)
            elseif(93<=prn .and. prn<=97)then
                sat=satno(SYS_QZS,prn+100)
            else
                sat=satno(SYS_GPS,prn)
            endif
        endif
        ! decode toc field 
        call str2time(buff(1+sp:),1,19,toc,info)
        if (info/=0)then
            stat=0; return  ! error
        endif
        ! decode data fields 
        p=buff(1+sp+19:)
        do j=1,3
            mydata(i)=str2num(p,1,19)
            i=i+1; p=p(20:)
        enddo
    else
        ! decode data fields 
        p=buff(1+sp:)
        do j=1,4
            mydata(i)=str2num(p,1,19)
            i=i+1; p=p(20:)
        enddo
        ! decode ephemeris 
        if(sys==SYS_GLO .and. i>15)then
            if (and(mask,sys)==0)then
                stat=0; return
            endif
            mytype=1
            call decode_geph(ver,sat,toc,mydata,geph,stat); return  ! 0-error, 1-normal
        elseif(i>31)then
            if (and(mask,sys)==0)then
                stat=0; return
            endif
            mytype=0
            call decode_eph(ver,sat,toc,mydata,eph,stat); return  ! 0-error, 1-normal
        endif
    endif
enddo
stat=-1
end subroutine

! unique ephemerides --------------------------------------------------------
! unique ephemerides in navigation data and update carrier wave length
! args   : nav_t *nav    IO     navigation data
! return : number of epochs
!----------------------------------------------------------------------------
subroutine uniqnav(nav)
implicit none
include 'file_para.h'
type(nav_t), intent(out) :: nav
integer*4 i,j
real*8, external :: satwavelen

! update carrier wave length 
do i=1,MAXSAT
    do j=1,NFREQ
        nav%lam(i,j)=satwavelen(i,j-1,nav)
    enddo
enddo
end subroutine

!! select ephememeris -------------------------------------------------------
function seleph(time, sat, iode, nav)
implicit none
include 'file_para.h'
type(eph_t), pointer :: seleph
type(gtime_t), intent(in) :: time
integer*4, intent(in) :: sat,iode
type(nav_t), intent(in) :: nav
real*8 t,tmax,tmin,timediff
integer*4 :: i,j
integer*4 prn, sys, sel
external :: timediff
j=-1; sel=0
call satsys(sat,prn,sys)
select case(sys)
case(SYS_GPS)
    tmax=MAXDTOE+1.d0    ; sel=eph_sel(1)
case(SYS_GAL)
    tmax=MAXDTOE_GAL     ; sel=eph_sel(3)
case(SYS_QZS)
    tmax=MAXDTOE_QZS+1.d0; sel=eph_sel(4)
case(SYS_CMP)
    tmax=MAXDTOE_CMP+1.d0; sel=eph_sel(5)
case default
    tmax=MAXDTOE+1.d0
end select
tmin=tmax+1.d0

! write(*,*) char(10),"sat_t=",t
do i=1,nav%n
    if (nav%eph(i)%sat/=sat) cycle
    if (iode>=0 .and. nav%eph(i)%iode/=iode) cycle
    if (sys==SYS_GAL .and. sel/=0)then
        if (sel==1 .and. and(nav%eph(i)%code,lshift(1,9))==0) cycle  ! I/NAV 
        if (sel==2 .and. and(nav%eph(i)%code,lshift(1,8))==0) cycle  ! F/NAV 
    endif
    t=dabs(timediff(nav%eph(i)%toe,time))
    ! write(*,*) "nav_t=",t  ! for test
    if (t>tmax) cycle
    if (iode>=0)then
        seleph=>nav%eph(i); return
    endif
    if (t<=tmin)then  ! toe closest to time 
        j=i; tmin=t
    endif
enddo
if (iode>=0 .or. j<0)then
    nullify(seleph); return
endif
seleph=>nav%eph(j)
end function

! select glonass ephememeris ------------------------------------------------
function selgeph(time, sat, iode, nav)
implicit none
include 'file_para.h'
type(geph_t), pointer :: selgeph
type(gtime_t), intent(in) :: time
integer*4, intent(in) :: sat,iode
type(nav_t), intent(in) :: nav
real*8 :: t,tmax,tmin,timediff
integer*4 :: i,j
external :: timediff
j=-1
tmax=MAXDTOE_GLO; tmin=tmax+1.d0
do i=1,nav%ng
    if (nav%geph(i)%sat/=sat) cycle
    if (iode>=0 .and. nav%geph(i)%iode/=iode) cycle
    t=dabs(timediff(nav%geph(i)%toe,time))
    if (t>tmax) cycle
    if (iode>=0)then
        selgeph=>nav%geph(i); return
    endif
    if (t<=tmin)then  ! toe closest to time 
        j=i; tmin=t
    endif
enddo
if (iode>=0 .or. j<0)then
    nullify(selgeph); return
endif
selgeph=>nav%geph(j)
end function
