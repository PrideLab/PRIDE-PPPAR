
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
character(63) :: helptext(42)=(/&
    "                                                               ",&
    "spp version 3.0,  based on open-source software rtklib.        ",&
    "  Wuhan University, Oct. 2024                                  ",&
    "                                                               ",&
    "Usage: spp [options] <arguments> -o outname rinexobs rinexnav  ",&
    "                                                               ",&
    "Description:                                                   ",&
    "  spp is a module of PRIDE PPP-AR, used for calculating the    ",&
    "  initial coordinates of the station using pesudo-range.       ",&
    "                                                               ",&
    "Required arguments:                                            ",&
    "  -o outname                                                   ",&
    "    output file.                                               ",&
    "  rinexobs                                                     ",&
    "    rinex obs file.                                            ",&
    "  rinexnav                                                     ",&
    "    rinex nav file.                                            ",&
    "                                                               ",&
    "Optional arguments:                                            ",&
    "  -elev [mask]                                                 ",&
    "    elevation mask. Default is 7. (Unit: degree)               ",&
    "  -trop [model]                                                ",&
    "    tropsphere correction model. Default is SAAS.              ",&
    "      NON = not correct                                        ",&
    "      SAAS = saastamoinen model                                ",&
    "  -ts [year/month/day hour:minute:second]                      ",&
    "    start time. Default is the start time in rinex file.       ",&
    "  -te [year/month/day hour:minute:second]                      ",&
    "    end time. Default is the end time in rinex file.           ",&
    "  -ti [interval]                                               ",&
    "    sampleing rate. Default is all records. (Unit: second)     ",&
    "  -twnd [window]                                               ",&
    "    processing time window in seconds, for not standard or     ",&
    "    super high rate rinex obs. Default is 0.01. (Unit: second) ",&
    "                                                               ",&
    "Examples:                                                      ",&
    "  spp -o kin_brux brux0010.24o brdm0010.24p                    ",&
    "  spp -elev 7 -trop SAAS -o kin_brux brux0010.24o brdm0010.24p ",&
    "                                                               ",&
    "More details refer to PRIDE PPP-AR manual and repository       ",&
    "  https://github.com/PrideLab/PRIDE-PPPAR/                     ",&
    "                                                               "/)
do i=1,size(helptext,dim=1)
    write(*,*) trim(helptext(i))
enddo
end subroutine
