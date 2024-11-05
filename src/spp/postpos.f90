
! Post-processing positioning -----------------------------------------------
integer*4 function postpos(ts, te, ti, twnd, popt, sopt, rnxobsfile, rnxnavs , fout)
implicit none
include 'file_para.h'
type(gtime_t), intent(in) :: ts, te
real*8, intent(inout) :: ti, twnd
type(prcopt_t), intent(in) :: popt
type(solopt_t), intent(in) :: sopt
character(*), intent(in) :: rnxobsfile
character(*), intent(in) :: rnxnavs(FLNUMLIM)
integer*4, intent(in) :: fout

! local 
real*8 :: ver, obs_int
integer*4 :: nsat, i, nobs, prn, sys, stat, dobs, dlst
character(3) :: tobs(NUMSYS,MAXOBSTYPE), satid
type(obsd_t) obs(MAXOBS)
character(MAXOBS*2*3) satrec
character*128 :: msg=''
real*8 azel(MAXOBS,2)
type(sol_t) sol0

type(ssat_t) :: ssat0(MAXSAT)

! function
integer*4, external :: outhead
real*8, external :: timediff

! initialize
stat=0
call init_sta(site)
tobs=''
ssat0=ssat_t(0,0,0.d0,0.d0,0.d0,0.d0,0,0,0)
sol0=sol_t(gtime_t(0,0.d0),gtime_t(0,0.d0),0.d0,0.d0,0.d0,0,0,0,0.d0)

! open rinex obs file
open(unit=FPREAD,file=rnxobsfile,status='old',iostat=stat)
if(stat/=0)then
    write(*,*) "File opening error : ", trim(rnxobsfile)
    postpos=-1; return
endif

! read rinex obs header 
call readrnxobsh(FPREAD, ver, tobs, site, obs_int, stat)
if (stat==0)then
    write(*,*) "Error reading rinex header ...", trim(rnxobsfile)
    postpos=-1; return
endif
! no site name in obs header
if(site%name=='') call decodemarker(rnxobsfile, site%name)
! don't set compute inertval
if(ti==0) ti=obs_int

! write header to output file 
if(sopt%issingle==0 .and. headwritten_==0)then
    call outheader(fout,site,ti)
    headwritten_=1
endif

do while(.true.)
    ! read obs epoch data
    call get_obsepoch(FPREAD,ts,te,ti,twnd,ver,tobs,obs,nobs)
    if(nobs<0) exit

    ! check result time stamp
    if(solindex_>0 .and. timediff(allsol_(solindex_)%time0,obs(1)%time)>=0) cycle

    ! first time read rinex nav
    if(brdm_idx==0) then
        brdm_idx=brdm_idx+1
        call readrnxnavfile(rnxnavs(brdm_idx),navs,stat)
        if(stat<=0)then
            write(*,*) "File reading error : ", trim(rnxnavs(i))
            exit
        endif
    ! shift rinex nav file according day diff
    elseif(timediff(t_prev,gtime_t(0,0.d0))>0.0)then
        dobs=int(obs(1)%time%time/86400) 
        dlst=int(t_prev%time/86400)
        if((dobs-dlst)>FLNUMLIM) then
            exit
        elseif((dobs-dlst)>0) then
            brdm_idx=brdm_idx+dobs-dlst
            call freenav(navs)
            call readrnxnavfile(rnxnavs(brdm_idx),navs,stat)
            if(stat<=0)then
                write(*,*) "File reading error : ", trim(rnxnavs(i))
                exit
            endif
        endif
    endif
    t_prev=obs(1)%time
    
    ! exclude satellites 
    nsat=0; satrec=""
    do i=1,nobs
        call satsys(obs(i)%sat,prn,sys)
        if(and(sys,popt%navsys)/=0 .and. popt%exsats(obs(i)%sat)/=1)then
            call satno2id(obs(i)%sat,satid)
            if(index(satrec,satid)==0)then
                obs(nsat+1)=obs(i); nsat=nsat+1
                satrec=trim(satrec)//satid
            endif
        endif
    enddo
    if (nsat<=0) cycle
    
    ! single point positioning
    call pntpos(obs,nsat,navs,popt,sol0,azel,ssat0,msg,stat)
    if (stat==0) cycle

    ! write positioning result to outfile
    if(sopt%issingle==0) call outsol(fout,sol0,sopt)

    ! add current sloution to all solutions
    call add_sol(sol0)
enddo

close(FPREAD)

if(solindex_==0)then
    write(*,*) "Info : no solution are computed, postpos()"
    postpos=0; return
endif

postpos=1
end function

! add solution to all solution -------------------------------------------------------
subroutine add_sol(sol)
implicit none
include 'file_para.h'
type(sol_t), intent(in) :: sol
type(sol_t), pointer :: sol_data(:)

solindex_=solindex_+1
if(solindex_>size(allsol_))then
    allocate(sol_data(size(allsol_)+MAXSOLNUM))
    sol_data(1:size(allsol_))=allsol_(1:size(allsol_))
    deallocate(allsol_)
    allsol_=>sol_data
    nullify(sol_data)
endif
allsol_(solindex_)=sol

end subroutine
