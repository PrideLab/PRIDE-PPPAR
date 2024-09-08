
! Solution ------------------------------------------------------------------
! write header to output file -----------------------------------------------
integer*4 function outhead(outfile, sta, interval)
implicit none
include 'file_para.h'
character(*), intent(in) :: outfile
type(sta_t), intent(in) :: sta
real*8, intent(in) :: interval

! local
integer*4 :: fp=6, info

if(len_trim(outfile)/=0)then
    fp=FPOUT
    open(unit=fp,file=trim(outfile),status='replace',iostat=info)
    if(info/=0)then
        write(unit=6,fmt="('Error : open output file ',A)") trim(outfile)
        outhead=0; return
    endif
endif
! output header 
call outheader(fp,sta,interval)
if(len_trim(outfile)/=0) close(unit=fp,status='keep')
outhead=1
end function

subroutine outheader(fp, sta, interval)
implicit none
include 'file_para.h'
integer*4, intent(in) :: fp
type(sta_t), intent(in) :: sta
real*8, intent(in) :: interval

if(fp/=6)then  ! disallow output to screen
    write(fp,"(A20,A10,A4,A26,A7)") "Kinematic Trajectory",'',trim(sta%name(1:4)),'',"COMMENT"
    write(fp,"(F9.2,A51,A8)") interval,'',"INTERVAL"
    write(fp,"(A60,A13)") '',"END OF HEADER"
endif
end subroutine

! output solution body ------------------------------------------------------
! output solution body to file
! args   : FILE   *fp       I   output file pointer
!          sol_t  *sol      I   solution
!          real*8 *rb       I   base station position {x,y,z} (ecef) (m)
!          solopt_t *opt    I   solution options
! return : none
!----------------------------------------------------------------------------
subroutine outsol(fp, sol, opt)
implicit none
include 'file_para.h'
integer*4, intent(in) :: fp
type(sol_t), intent(in) :: sol
type(solopt_t), intent(in) :: opt

character(256) buff
integer*4 n
call outsols(buff,sol,opt,n)

if(n>0)then
    if(fp/=6) write(fp,"(A)") trim(buff)  ! disallow output to screen
endif
end subroutine

! output solution body ------------------------------------------------------
! output solution body to buffer
! args   : character(1) *buff IO output buffer
!          sol_t  *sol      I   solution
!          real*8 *rb       I   base station position {x,y,z} (ecef) (m)
!          solopt_t *opt    I   solution options
! return : number of output bytes
!----------------------------------------------------------------------------
subroutine outsols(buff, sol, opt, stat1)
implicit none
include 'file_para.h'
character(*), intent(out) :: buff
type(sol_t), intent(in) :: sol
type(solopt_t), intent(in) :: opt
integer*4, intent(out) :: stat1
type(gtime_t) time,gpst2utc,timeadd
real*8 sod,sol_std,ROUNDF
integer*4 :: mjd,timeu,info
character(64) s,fmtstr
external :: sol_std,gpst2utc,timeadd,ROUNDF
timeu=0
! suppress output if std is over opt%maxsolstd 
if(opt%maxsolstd>0.d0 .and. sol_std(sol)>opt%maxsolstd)then
    stat1=0; return
endif
! align the timestamp with observations insead of estimates (2022, Jihang Lin)
!   reason: avoid mismatching timestamps with float seconds
time=sol%time0
call time2mjd(time,mjd,sod)
if (opt%times>=TIMES_UTC) time=gpst2utc(time)
if (opt%times==TIMES_JST) time=timeadd(time,9*3600.d0)
call time2mjd(time,mjd,sod)
if(dabs(86400.d0-sod)<1d-3)then
    mjd=mjd+1; sod=0.d0
endif
fmtstr="(I5,F10.2)"
write(s,trim(fmtstr)) mjd,sod
call outecef(buff,s,sol,opt,info)
stat1=len_trim(buff)  ! p-(char*)buff
end subroutine

! output ecef (custom format 3) ---------------------------------------------
subroutine outecef(buff, s, sol, opt, stat1)
implicit none
include 'file_para.h'
character(*), intent(out) :: buff
character(*), intent(in) :: s
type(sol_t), intent(in) :: sol
type(solopt_t), intent(in) :: opt
integer*4, intent(out) :: stat1
character(1) :: sep=char(9)  !'\t'
character(9) solqr(6)
character(10) soltmp
integer*4 i
write(buff,"(A,A2,3(F14.4))") trim(s),'',sol%rr(1),sol%rr(2),sol%rr(3)
stat1=len_trim(buff)
end subroutine

! std-dev of soltuion -------------------------------------------------------
real*8 function sol_std(sol)
implicit none
include 'file_para.h'
type(sol_t), intent(in) :: sol
real*8, external :: SQRT2
! approximate as max std-dev of 3-axis std-devs 
if (sol%qr(1)>sol%qr(2) .and. sol%qr(1)>sol%qr(3))then
    sol_std=SQRT2(sol%qr(1)); return
endif
if (sol%qr(2)>sol%qr(3))then
    sol_std=SQRT2(sol%qr(2)); return
endif
sol_std=SQRT2(sol%qr(3))
end function
