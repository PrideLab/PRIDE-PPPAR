subroutine readrnxobsh(fp, ver, tobs, station, tint, stat)
implicit none
include 'file_para.h'
integer*4, intent(in) :: fp
real*8, intent(out) :: ver
character(3), intent(out) :: tobs(NUMSYS,MAXOBSTYPE)
type(sta_t), intent(out) :: station
real*8, intent(out) :: tint
integer*4, intent(out) :: stat

! local
character(1024) :: buff
character(20) :: label
integer*4 :: info, tsys

ver=2.10
stat=-1
do while(.true.)
    read(fp,"(A)",iostat=info) buff
    if(info/=0) exit
    label=buff(61:)
    if (len_trim(buff)<=60)then
        cycle
    elseif(index(label,'RINEX VERSION / TYPE')/=0)then
        read(buff(1:9),*,iostat=info) ver
        cycle
    elseif(index(label,'PGM / RUN BY / DATE')/=0)then
        cycle
    elseif(index(label,'COMMENT')/=0)then
        cycle
    endif
    call decode_obsh(fp,buff,ver,tsys,tobs,tint,station)
    if (index(label,'END OF HEADER')/=0) then
        stat=1
        return
    endif
enddo
stat=-1
end subroutine

! Read one epoch rinex obs ----------------------------------------------------
subroutine get_obsepoch(fp, ts, te, tint, ver, tobs, obs, stat)
implicit none
include 'file_para.h'
integer*4, intent(in) :: fp
type(gtime_t), intent(in) :: ts, te
real*8, intent(in) :: tint, ver
character(3), intent(in) :: tobs(NUMSYS,MAXOBSTYPE)
type(obsd_t), intent(out) :: obs(MAXOBS)
integer*4, intent(out) :: stat  ! 0-error, 1-normal

! local
integer*4 :: line,nsat,flag,sats(MAXOBS),info,n
real*8 ti
character(MAXRNXLEN) buff
type(sigind_t) :: myindex(7)
type(gtime_t) :: time

! function
type(gtime_t), external :: utc2gpst
integer*4, external :: screent

! initialize variables
flag=0; ti=0
ti=dabs(tint)
time=gtime_t(0,0.d0)
myindex=sigind_t(0,0,0,0,0,0,0.d0)
line=0; n=1; nsat=0; sats(MAXOBS)=0

! set signal index 
call set_index(0.d0,SYS_GPS,'',tobs(1,:),myindex(1))
call set_index(0.d0,SYS_GLO,'',tobs(2,:),myindex(2))
call set_index(0.d0,SYS_GAL,'',tobs(3,:),myindex(3))
call set_index(0.d0,SYS_QZS,'',tobs(4,:),myindex(4))
call set_index(0.d0,SYS_SBS,'',tobs(5,:),myindex(5))
call set_index(0.d0,SYS_CMP,'',tobs(6,:),myindex(6))
call set_index(0.d0,SYS_IRN,'',tobs(7,:),myindex(7))

! read record 
do while(.true.)
    read(fp,"(A)",iostat=info) buff
    if(info/=0) exit
    if(index(buff,"RINEX VERSION / TYPE")/=0)then
        do while(.true.)
            read(fp,"(A)",iostat=info) buff
            if(info/=0) exit
            if(index(buff,'END OF HEADER')/=0)exit
        enddo
    endif
    ! decode obs epoch 
    if(line==0)then
        call decode_obsepoch(fp,buff,ver,time,flag,sats,nsat)
        if(nsat<=0) cycle
    elseif(flag<=2 .or. flag==6)then
        obs(n)%time=time
        obs(n)%sat=sats(line)
        obs(n)%rcv=1
        ! decode obs data 
        if(nsat>0 .and. screent(time,ts,te,ti)==0)then
            call cycle_obsdata(fp,ver,myindex)
        else
            call decode_obsdata(fp,buff,ver,SYS_ALL,myindex,obs(n),info)
        endif
        if(info/=0 .and. n<=MAXOBS) n=n+1
    elseif(flag==3 .or. flag==4)then
        backspace fp
    endif
    line=line+1
    if(line>nsat .and. screent(time,ts,te,ti)==0)then
        line=0
    elseif(line>nsat .and. screent(time,ts,te,ti)==1)then
        stat=n-1
        return
    endif
enddo
stat=-1

end subroutine

! Read rinex obs ------------------------------------------------------------
subroutine readrnxobs(fp, ts, te, tint, opt, rcv, ver, tsys, tobs, obs, sta, stat)
implicit none
include 'file_para.h'
integer*4, intent(in) :: fp, rcv, tsys
type(gtime_t), intent(in) :: ts, te
real*8, intent(in) :: tint, ver
character(*), intent(in) :: opt
character(3), intent(in) :: tobs(NUMSYS,MAXOBSTYPE)  !tobs(:,:)
type(obs_t), intent(out) :: obs
type(sta_t), intent(in) :: sta
integer*4, intent(out) :: stat  ! 0-error, 1-normal
type(obsd_t), pointer :: mydata(:)
integer*4 :: i,n,flag,stat1,slips(MAXSAT,NFREQ)
real*8 ti
type(gtime_t), external :: utc2gpst
integer*4, external :: screent
flag=0; stat1=0; slips=0; ti=0
if(rcv>MAXRCV)then
    stat=0; return
endif
allocate(mydata(MAXOBS))

if(tint/=0) obs%tint=dabs(tint)
ti=dabs(tint)

! read rinex obs data body 
do while(.true.)
    call readrnxobsb(fp,opt,ver,tsys,tobs,flag,mydata,sta,n)  ! n=-1: error
    if(.not.(n>=0 .and. stat1>=0)) exit
    do i=1,n
        ! utc % gpst 
        if(tsys==TSYS_UTC) mydata(i)%time=utc2gpst(mydata(i)%time)
        ! save cycle-slip 
        call saveslips(slips,mydata(i))
    enddo
    if(n>0 .and. screent(mydata(1)%time,ts,te,ti)==0) cycle
    do i=1,n
        ! restore cycle-slip 
        call restslips(slips,mydata(i))
        ! save obs data 
        mydata(i)%rcv=1
        call addobsdata(obs,mydata(i),stat1)  ! 1-normal
        if (stat1<0) exit
    enddo
enddo
deallocate(mydata)
stat=stat1
end subroutine

! read rinex obs data body --------------------------------------------------
subroutine readrnxobsb(fp, opt, ver, tsys, tobs, flag, mydata, sta, stat)
implicit none
include 'file_para.h'
integer*4, intent(in) :: fp, tsys
character(*), intent(in) :: opt
real*8, intent(in) :: ver
character(3), intent(in) :: tobs(NUMSYS,MAXOBSTYPE)  !tobs(:,:)
integer*4, intent(out) :: flag, stat  ! -1: error
type(obsd_t), intent(out) :: mydata(MAXOBS)  !mydata(:)
type(sta_t), intent(in) :: sta
type(gtime_t) :: time
type(sigind_t) :: myindex(7)
type(obs_t) obstmp
type(nav_t) navtmp
character(MAXRNXLEN) buff
integer*4 :: i,n,nsat,sats(MAXOBS),mask,info,set_sysmask
external :: set_sysmask
time=gtime_t(0,0.d0)
myindex=sigind_t(0,0,0,0,0,0,0.d0)
i=0; n=1; nsat=0; sats(MAXOBS)=0

! set system mask 
mask=set_sysmask(opt)
! set signal index 

call set_index(ver,SYS_GPS,opt,tobs(1,:),myindex(1))
call set_index(ver,SYS_GLO,opt,tobs(2,:),myindex(2))
call set_index(ver,SYS_GAL,opt,tobs(3,:),myindex(3))
call set_index(ver,SYS_QZS,opt,tobs(4,:),myindex(4))
call set_index(ver,SYS_SBS,opt,tobs(5,:),myindex(5))
call set_index(ver,SYS_CMP,opt,tobs(6,:),myindex(6))
call set_index(ver,SYS_IRN,opt,tobs(7,:),myindex(7))

! read record 
do while(.true.)
    read(fp,"(A)",iostat=info) buff
    !if(info/=0 .or. (buff=='' .and. ver>3d0)) exit  ! to be verified
    if(info/=0) exit
    if(index(buff,"RINEX VERSION / TYPE")/=0)then
        do while(.true.)
            ! decode obs header 
            ! call decode_obsh(fp,buff,ver,tsys,tobs,obstmp,navtmp,sta)
            read(fp,"(A)",iostat=info) buff
            if(info/=0)exit
            if(index(buff,'END OF HEADER')/=0)exit
        enddo
    endif
    ! decode obs epoch 
    if(i==0)then
        call decode_obsepoch(fp,buff,ver,time,flag,sats,nsat)
        if(nsat<=0)cycle
    elseif(flag<=2 .or. flag==6)then
        mydata(n)%time=time
        mydata(n)%sat=sats(i)
        ! decode obs data 
        call decode_obsdata(fp,buff,ver,mask,myindex,mydata(n),info)
        if(info/=0 .and. n<=MAXOBS) n=n+1
    elseif(flag==3 .or. flag==4)then
        backspace fp
    endif
    i=i+1
    if(i>nsat)then
        stat=n-1; return
    endif
enddo
stat=-1
end subroutine
