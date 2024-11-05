
! Main function entry -------------------------------------------------------
program spp
implicit none
include 'file_para.h'
type(gtime_t) :: ts, te
real*8 :: tint, twnd
type(prcopt_t) :: prcopt
type(solopt_t) :: solopt
integer*4 :: nrnxo, nrnxn
character(1024) :: outfile, rnxobslist(FLNUMLIM), rnxnavlist(FLNUMLIM)

! Initialize variables
call InitGlobal()
prcopt=prcopt_default
solopt=solopt_default
ts=gtime_t(0,0.d0)
te=gtime_t(0,0.d0)
tint=0.d0
rnxobslist=""
rnxnavlist=""
outfile=""
nrnxn=0
nrnxo=0

! Get command line arguments
call get_spp_args(ts, te, tint, twnd, prcopt, rnxobslist, nrnxo, rnxnavlist, nrnxn, outfile)

call rnx2rtkp(ts, te, tint, twnd, prcopt, solopt, rnxobslist, nrnxo, rnxnavlist, nrnxn, outfile)

end program spp

! rnx2rtkp ------------------------------------------------------------------
subroutine rnx2rtkp(ts, te, tint, twnd, prcopt, solopt, rnxobslist, nrnxo, rnxnavlist, nrnxn, outfile)
implicit none
include 'file_para.h'
type(gtime_t),intent(in) :: ts, te
real*8,intent(in) :: tint, twnd
type(prcopt_t),intent(in) :: prcopt
type(solopt_t),intent(in) :: solopt
character(1024),intent(in) :: rnxobslist(FLNUMLIM)
integer*4,intent(in) :: nrnxo
character(1024),intent(in) :: rnxnavlist(FLNUMLIM)
integer*4,intent(in) :: nrnxn
character(1024),intent(in) :: outfile

! local
integer*4 :: i, ret, info
character(1024) :: filetmp
logical*1 :: isexceed

! function
integer*4 , external :: postpos
real*8, external :: timediff

! open output file
open(unit=FPOUT,file=outfile,status='REPLACE',iostat=info)
if(info/=0)then
    write(*,*) "Error : open output file ", trim(outfile)
    return
endif

isexceed=.false.
do i=1,FLNUMLIM
    ! cancel loop
    if(isexceed) exit

    ! check if obs file is exist
    call getfname(rnxobslist(i),filetmp)
    if(len_trim(filetmp)==0) cycle

    ret=postpos(ts,te,tint,twnd,prcopt,solopt,rnxobslist(i),rnxnavlist,FPOUT)  ! 0-error, 1-right
    if(ret==0) exit

    if (timediff(te,gtime_t(0,0.d0))>0.0 .and. timediff(allsol_(solindex_)%time,te) >= 0)then
        isexceed=.true.
    ! don't set te, only read one rinex obs
    elseif (te%time==0) then 
        isexceed=.true.
    endif
enddo

close(FPOUT)

! print the final result
call printresult(solopt,ret)
if(ret==0)then
    write(unit=6,fmt="(A40)") ''; call exit(1)
endif
if(ret==1) call exit(0)

end subroutine
