
! Information display -------------------------------------------------------
! print the final result ----------------------------------------------------
subroutine printresult(sopt, ret)
implicit none
include 'file_para.h'
type(solopt_t), intent(in) :: sopt
integer*4, intent(out) :: ret  ! 0-error, 1-right
real*8, allocatable :: rr(:,:)
real*8 avexyz(3)
type(gtime_t) tss, tee
real*8 epss(6), timespan
real*8, external :: raver, timediff
avexyz=0

if(solindex_==0)then
    ret=0; return
endif

allocate(rr(solindex_,3))
call errorfilter(allsol_(1:solindex_)%rr(1),rr(:,1),solindex_)
call errorfilter(allsol_(1:solindex_)%rr(2),rr(:,2),solindex_)
call errorfilter(allsol_(1:solindex_)%rr(3),rr(:,3),solindex_)
avexyz(1)=raver(rr(:,1),solindex_)
avexyz(2)=raver(rr(:,2),solindex_)
avexyz(3)=raver(rr(:,3),solindex_)
deallocate(rr)

if(avexyz(1)==0 .and. avexyz(2)==0 .and. avexyz(3)==0)then
    ret=0; return
endif
write(unit=6,fmt="(A11,3F16.4)") "Position : ", avexyz(1), avexyz(2), avexyz(3)

! output time stamp
tss=allsol_(1)%time0
tee=allsol_(solindex_)%time0
call time2epoch(tss, epss)
timespan=timediff(tee, tss)
write(unit=6,fmt="(A11,I4.4,' ',I2.2,' ',I2.2,' ',I2.2,' ',I2.2,' ',F5.2,' 'F11.2)") &
        "Duration : ", int(epss(1)), int(epss(2)), int(epss(3)), &
                        int(epss(4)), int(epss(5)), epss(6), timespan  !, round(timespan)
! endif
ret=1
end subroutine

subroutine printhelp()
implicit none
include 'file_para.h'
integer*4 i
character(63) :: helptext(11)=(/&
    "                                                               ",&
    " usage: spp [option]... file file [...]                        ",&
    "                                                               ",&
    " -?/-h      print help                                         ",&
    " -o file    set output file [stdout]                           ",&
    " -ts ds ts  start day/time (ds=y/m/d ts=h:m:s) [obs start time]",&
    " -te de te  end day/time   (de=y/m/d te=h:m:s) [obs end time]  ",&
    " -ti tint   time interval  (sec) [all]                         ",&
    " -trop      trop option    (non saas) [saas]                   ",&
    " -elev      minimum elevation (deg) [10]                       ",&
    "                                                               "/)
do i=1,size(helptext,dim=1)
    write(*,*) trim(helptext(i))
enddo
end subroutine
