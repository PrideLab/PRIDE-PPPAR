!! reference and acknowledgement to RTKLIB
!! modified by Kunlun Zhang
! constants and macros ------------------------------------------------------

module solution_f90_
use rtkcmn_f90_
implicit none

integer*4, parameter :: MAXFIELD = 64           ! max number of fields in a record 
integer*4, parameter :: MAXNMEA  = 256          ! max length of nmea sentence 
real*8, parameter    :: KNOT2M = 0.514444444d0  ! m/knot 

integer*4, parameter :: solq_nmea(10)=(/&       ! nmea quality flags to rtklib sol quality 
    SOLQ_NONE ,SOLQ_SINGLE, SOLQ_DGPS, SOLQ_PPP , SOLQ_FIX,&
    SOLQ_FLOAT,SOLQ_DR    , SOLQ_NONE, SOLQ_NONE, SOLQ_NONE/)

contains
subroutine outecef(buff, s, sol, opt, stat1)
implicit none
character(*), intent(out) :: buff
character(*), intent(in) :: s
type(sol_t), intent(in) :: sol
type(solopt_t), intent(in) :: opt
integer*4, intent(out) :: stat1
character(1) :: sep=char(9)  !'\t'
character(9) solqr(6)
character(10) soltmp
integer*4 i
write(buff,"(A,3(A,F14.4),2(A,I3))") trim(s), sep,sol%rr(1), sep,sol%rr(2), sep,sol%rr(3), sep,sol%stat, sep,sol%ns
do i=1,3
    if(SQRT2(sol%qr(i))<=-1000)then
        write(solqr(i),"(F10.4)") SQRT2(sol%qr(i))
        buff=trim(buff)//sep//solqr(i)
    elseif(SQRT2(sol%qr(i))<=-100 .or. SQRT2(sol%qr(i))>=1000)then
        write(solqr(i),"(F9.4)") SQRT2(sol%qr(i))
        buff=trim(buff)//sep//solqr(i)
    else
        write(solqr(i),"(F8.4)") SQRT2(sol%qr(i))
        buff=trim(buff)//sep//solqr(i)
    endif
enddo

do i=4,6
    if(sqvar(sol%qr(i))<=-1000)then
        write(soltmp,"(F10.4)") sqvar(sol%qr(i))
        buff=trim(buff)//sep//soltmp
    elseif(sqvar(sol%qr(i))<=-100 .or. sqvar(sol%qr(i))>=1000)then
        write(solqr(i),"(F9.4)") sqvar(sol%qr(i))
        buff=trim(buff)//sep//solqr(i)
    else
        write(solqr(i),"(F8.4)") sqvar(sol%qr(i))
        buff=trim(buff)//sep//solqr(i)
    endif
enddo
!if(opt%outvel)   ! output velocity 
stat1=len_trim(buff)
end subroutine

! output ecef (custom format 3) -----------------------------------------------
subroutine outecef2(buff, s, sol, opt, stat1)
implicit none
character(*), intent(out) :: buff
character(*), intent(in) :: s
type(sol_t), intent(in) :: sol
type(solopt_t), intent(in) :: opt
integer*4, intent(out) :: stat1
character(1) :: sep=char(9)  !'\t'
character(9) solqr(6)
character(10) soltmp
integer*4 i
write(buff,"(A,A2,3(F13.3))") trim(s),'',sol%rr(1),sol%rr(2),sol%rr(3)
stat1=len_trim(buff)
end subroutine

! output processing options ---------------------------------------------------
! output processing options to buffer
! args   : unsigned char *buff IO output buffer
!          prcopt_t *opt   I   processign options
! return : number of output bytes
!------------------------------------------------------------------------------
subroutine outprcopts(buff, opt, stat1)
implicit none
character(*), intent(out) :: buff
type(prcopt_t), intent(in) :: opt
integer*4, intent(out) :: stat1

integer*4 :: sys(8)=(/SYS_GPS,SYS_GLO,SYS_GAL,SYS_QZS,SYS_CMP,SYS_IRN,SYS_SBS,0/), i
character(100) :: s1(10)=(/"single       ","dgps         ","kinematic    ",&
                           "static       ","moving-base  ","fixed        ",&
                           "ppp-kinematic","ppp-static   ","ppp-fixed    ",&
                           "             "/)
character(100) :: s4(12)=(/"off          ","broadcast    ","sbas         ",&
                           "iono-free    ","estimation   ","ionex tec    ",&
                           "qzs          ","lex          ","vtec_sf      ",&
                           "vtec_ef      ","gtec         ","             "/)
character(100) :: s5(6) =(/"off          ","saastamoinen ","sbas         ",&
                           "est ztd      ","est ztd+grad ","             "/)
character(100) :: s6(7) =(/"broadcast        ","precise          ","broadcast+sbas   ",&
                           "broadcast+ssr com","broadcast+ssr com","qzss lex         ",&
                           "                 "/)
character(100) :: s7(8) =(/"gps          ","glonass      ","galileo      ",&
                           "qzss         ","beidou       ","irnss        ",&
                           "sbas         ","             "/)

write(buff                   , "(A,' pos mode  : ',2A)") COMMENTH, trim(s1(opt%mode+1)), char(10)  !'\n'
write(buff(len_trim(buff)+1:), "(A,' elev mask : ',F4.1,' deg',A)") COMMENTH, opt%elmin*R2D,char(10)
write(buff(len_trim(buff)+1:), "(A,' ionos opt : ',2A)") COMMENTH, trim(s4(opt%ionoopt+1)), char(10)
write(buff(len_trim(buff)+1:), "(A,' tropo opt : ',2A)") COMMENTH, trim(s5(opt%tropopt+1)), char(10)
write(buff(len_trim(buff)+1:), "(A,' ephemeris : ',2A)") COMMENTH, trim(s6(opt%sateph +1))

if(opt%navsys/=SYS_GPS)then
    write(buff(len_trim(buff)+1:), "(A,A,' navi sys  :')") char(10), COMMENTH
    do i=1,size(sys,dim=1)
        if(sys(i)==0) exit
        if(and(opt%navsys,sys(i))/=0) write(buff(len_trim(buff)+1:), "(' ',A)") trim(s7(i))
    enddo
    buff=trim(buff)  !//char(10)
endif
!write(buff(len_trim(buff)+1:),"(A)") COMMENTH
stat1=len_trim(buff)
end subroutine

! output solution header ------------------------------------------------------
! output solution header to buffer
! args   : character(1) *buff IO output buffer
!          solopt_t *opt  I   solution options
! return : number of output bytes
!------------------------------------------------------------------------------
subroutine outsolheads(buff, opt, stat1)
implicit none
character(*), intent(out) :: buff
type(solopt_t), intent(in) :: opt
integer*4, intent(out) :: stat1
character(1) :: sep=char(9)
if(opt%outhead/=0)then
    write(buff,"(A1,' (')") COMMENTH
    buff=trim(buff)//"x/y/z-ecef=WGS84"
    buff=trim(buff)//",Q=1:fix,2:float,3:sbas,4:dgps,5:single,6:ppp,ns=# of satellites)"
    buff=trim(buff)//char(10)
endif
write(buff(len_trim(buff)+1:),"(A,A6,2A,3(A14A),2(A3A),5(A8A),A8)")&
    COMMENTH,"GPST",sep,sep,'x-ecef(m)',sep,'y-ecef(m)',sep,'z-ecef(m)',&
    sep,'Q',sep,'ns',sep,'sdx(m)',sep,'sdy(m)',sep,'sdz(m)',sep,'sdxy(m)',&
    sep,'sdyz(m)',sep,'sdzx(m)'  !,sep,'age(s)',sep,'ratio'
!buff=trim(buff)//char(10)
stat1=len_trim(buff)
end subroutine

! std-dev of soltuion -------------------------------------------------------
real*8 function sol_std(sol)
implicit none
type(sol_t), intent(in) :: sol

! approximate as max std-dev of 3-axis std-devs 
if (sol%qr(1)>sol%qr(2) .and. sol%qr(1)>sol%qr(3))then
    sol_std=SQRT2(sol%qr(1)); return
endif
if (sol%qr(2)>sol%qr(3))then
    sol_std=SQRT2(sol%qr(2)); return
endif
sol_std=SQRT2(sol%qr(3))
end function

! output solution body --------------------------------------------------------
! output solution body to buffer
! args   : character(1) *buff IO output buffer
!          sol_t  *sol      I   solution
!          real*8 *rb       I   base station position {x,y,z} (ecef) (m)
!          solopt_t *opt    I   solution options
! return : number of output bytes
!------------------------------------------------------------------------------
subroutine outsols(buff, sol, rb, opt, stat1)
implicit none
character(*), intent(out) :: buff
type(sol_t), intent(in) :: sol
real*8, intent(in) :: rb(*)
type(solopt_t), intent(in) :: opt
integer*4, intent(out) :: stat1
type(gtime_t) time
real*8 gpst
integer*4 :: week,timeu,info
character(1) :: sep=char(9)
character(64) s,fmtstr
timeu=0
! suppress output if std is over opt%maxsolstd 
if(opt%maxsolstd>0.d0 .and. sol_std(sol)>opt%maxsolstd)then
    stat1=0; return
endif
time=sol%time
if (opt%times>=TIMES_UTC) time=gpst2utc(time)
if (opt%times==TIMES_JST) time=timeadd(time,9*3600.d0)
call time2gpst(time,week,gpst)
if(86400.d0*7-gpst<0.5/(10.d0**timeu))then
    week=week+1; gpst=0.d0
endif
if(timeu<=0)then
    fmtstr="(I6,A,I6)"
    write(s,trim(fmtstr)) week,sep,ROUND(gpst)
else
    write(fmtstr,"('(I4,A,F',I1,'.',I1,')')") 6+icond(timeu<=0,0,timeu+1),timeu
    write(s,trim(fmtstr)) week,sep,gpst
endif
call outecef(buff,s,sol,opt,info)
stat1=len_trim(buff)  ! p-(char*)buff
end subroutine

! output solution body (custom format 2) -------------------------------------
subroutine outsols2(buff, sol, rb, opt, stat1)  ! ignore 3
implicit none
character(*), intent(out) :: buff
type(sol_t), intent(in) :: sol
real*8, intent(in) :: rb(*)
type(solopt_t), intent(in) :: opt
integer*4, intent(out) :: stat1
type(gtime_t) time
real*8 sod
integer*4 :: mjd,timeu,info
character(1) :: sep=char(9)
character(64) s,fmtstr
timeu=0
! suppress output if std is over opt%maxsolstd 
if(opt%maxsolstd>0.d0 .and. sol_std(sol)>opt%maxsolstd)then
    stat1=0; return
endif
time=sol%time
if (opt%times>=TIMES_UTC) time=gpst2utc(time)
if (opt%times==TIMES_JST) time=timeadd(time,9*3600.d0)
!call time2gpst(time,week,gpst)
call time2mjd(time,mjd,sod)
if(86400.d0-sod<0.5/(10.d0**timeu))then
    mjd=mjd+1; sod=0.d0
endif
fmtstr="(I5,F9.2)"
write(s,trim(fmtstr)) mjd,ROUNDF(sod,2)
call outecef2(buff,s,sol,opt,info)
stat1=len_trim(buff)  ! p-(char*)buff
end subroutine

! output processing option ---------------------------------------------------
! output processing option to file
! args   : FILE   *fp       I   output file pointer
!          prcopt_t *opt    I   processing options
! return : none
!-----------------------------------------------------------------------------
subroutine outprcopt(fp, opt)
implicit none
integer*4, intent(in) :: fp
type(prcopt_t), intent(in) :: opt
character(MAXSOLMSG+1) buff
integer*4 n
call outprcopts(buff,opt,n)
if(n>0)then
    write(fp,"(A)") trim(buff)
endif
end subroutine
 
! output solution header ------------------------------------------------------
! output solution heade to file
! args   : FILE   *fp       I   output file pointer
!          solopt_t *opt    I   solution options
! return : none
!-----------------------------------------------------------------------------
subroutine outsolhead(fp, opt)
implicit none
integer*4, intent(in) :: fp
type(solopt_t), intent(in) :: opt
character(255) buff
integer*4 n
call outsolheads(buff,opt,n)
if(n>0)then
    write(fp,"(A)") trim(buff)
endif
end subroutine

! output solution body --------------------------------------------------------
! output solution body to file
! args   : FILE   *fp       I   output file pointer
!          sol_t  *sol      I   solution
!          real*8 *rb       I   base station position {x,y,z} (ecef) (m)
!          solopt_t *opt    I   solution options
! return : none
!-----------------------------------------------------------------------------
subroutine outsol(fp, sol, rb, opt)    ! change the output format here_2
implicit none
integer*4, intent(in) :: fp
type(sol_t), intent(in) :: sol
real*8, intent(in) :: rb(*)
type(solopt_t), intent(in) :: opt

character(256) buff
integer*4 n
!call outsols(buff,sol,rb,opt,n)
call outsols2(buff,sol,rb,opt,n)

if(n>0)then
    write(fp,"(A)") trim(buff)
endif
end subroutine
end module solution_f90_