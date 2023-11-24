
! Input obs data, navigation messages and sbas correction -------------------
subroutine inputobs(obs, solq, popt, stat)
implicit none
include 'file_para.h'
type(obsd_t), intent(out) :: obs(MAXOBS*2)
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
        if(obss%mydata(iobsu+1+i)%sat<=MAXSAT)then
            obs(n+1)=obss%mydata(iobsu+1+i); n=n+1
        endif
        i=i+1
    enddo
    iobsu=iobsu+nu
endif
stat=n
end subroutine

! search next observation data index ----------------------------------------
subroutine nextobsf(obs, i, rcv, stat)
implicit none
include 'file_para.h'
type(obs_t), intent(in) :: obs
integer*4, intent(out) :: i, stat
integer*4, intent(in) :: rcv
real*8 tt,timediff
integer*4 n
external :: timediff
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
include 'file_para.h'
type(obs_t), intent(in) :: obs
integer*4, intent(out) :: i, stat
integer*4, intent(in) :: rcv
real*8 tt,timediff
integer*4 n
external :: timediff
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
