
! Post-processing positioning -----------------------------------------------
integer*4 function postpos(ts, te, ti, tu, popt, sopt, infile, n, outfile, rov, base)
implicit none
include 'file_para.h'
type(gtime_t), intent(in) :: ts, te
real*8, intent(in) :: ti, tu
type(prcopt_t), intent(in) :: popt
type(solopt_t), intent(in) :: sopt
character(*), intent(in) :: infile(n), outfile, rov, base
integer*4, intent(in) :: n
integer*4 :: i,stat,myindex(MAXINFILE)
integer*4, external :: execses
stat=0; myindex=0
do i=1,n
    myindex(i)=i-1
enddo
! execute processing session 
stat=execses(ts,te,ti,popt,sopt,1,infile,myindex,n,outfile)
postpos=stat
end function

! execute processing session ------------------------------------------------
integer*4 function execses(ts, te, ti, popt, sopt, flag, infile, myindex, n, outfile)  ! 0-error, 1-right
implicit none
include 'file_para.h'
type(gtime_t), intent(in) :: ts, te
real*8, intent(in) :: ti
type(prcopt_t), intent(in) :: popt
type(solopt_t), intent(in) :: sopt
integer*4, intent(in) :: flag, myindex(MAXINFILE), n
character(*), intent(in) :: infile(n), outfile
integer*4 fp, info, outhead, openfile
external :: outhead, openfile
type(rtk_t) rtk
type(prcopt_t) popt_
type(gtime_t) stamp_
popt_=popt

! read obs and nav data 
!stamp_=timeget(); write(*,*) "Before readobsnav : ", stamp_%time*1d0+stamp_%sec
call readobsnav(ts,te,ti,infile,myindex,n,popt_,obss,navs,stas,info)
!stamp_=timeget(); write(*,*) "After readobsnav  : " , stamp_%time*1d0+stamp_%sec

if (info==0)then
    write(*,*) "Function execution failure, readobsnav()"
    execses=0; return
endif
! write header to output file 
if(flag/=0 .and. sopt%issingle==0 .and. headwritten_==0)then
    if(outhead(outfile,infile,n,popt_,sopt,stas(1),obss)==0)then
        call freeobsnav(obss,navs)
        write(*,*) "Function execution failure, outhead()"
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
        write(*,*) "Function execution failure, procpos()"
        execses=0; return
    endif
endif
!stamp_=timeget(); write(*,*) "After procpos :     " , stamp_%time*1d0+stamp_%sec

! free obs and nav data 
call freeobsnav(obss,navs)
execses=1
end function
