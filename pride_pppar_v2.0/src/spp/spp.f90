!! reference and acknowledgement to RTKLIB
!! modified by Kunlun Zhang
!! This is a module for standard point positioning which is modified and rewritten in Fortran by Kunlun Zhang with reference to RTKLIB. 
!! We thank RTKLIB and the copyright still belongs to RTKLIB.
! constants -----------------------------------------------------------------

module spp_all_
use postpos_f90_
implicit none

character(8), parameter  :: PROGNAME = 'rnx2rtkp'      ! program name 
integer*4, parameter     :: MAXFILE  = 16              ! max number of input files 

character(10), parameter :: SWTOPT = '0:off,1:on'
character(108),parameter :: MODOPT = '0:single,1:dgps,2:kinematic,3:static,4:static-start,5:movingbase,6:fixed,7:ppp-kine'  !,8:ppp-static,9:ppp-fixed'
character(31), parameter :: FRQOPT = '1:l1,2:l1+l2,3:l1+l2+l5,4:l1+l5'
character(31), parameter :: TYPOPT = '0:forward,1:backward,2:combined'
character(82), parameter :: IONOPT = '0:off,1:brdc,2:sbas,3:dual-freq,4:est-stec,5:ionex-tec,6:qzs-brdc,7:qzs-lex,8:stec'
character(49), parameter :: TRPOPT = '0:off,1:saas,2:sbas,3:est-ztd,4:est-ztdgrad,5:ztd'
character(56), parameter :: EPHOPT = '0:brdc,1:precise,2:brdc+sbas,3:brdc+ssrapc,4:brdc+ssrcom'
character(39), parameter :: NAVOPT = '1:gps+2:sbas+4:glo+8:gal+16:qzs+32:comp'
character(20), parameter :: GAROPT = '0:off,1:on,2:autocal'
character(24), parameter :: SOLOPT = '0:llh,1:xyz,2:enu,3:nmea'
character(18), parameter :: TSYOPT = '0:gpst,1:utc,2:jst'
character(11), parameter :: TFTOPT = '0:tow,1:hms'
character(11), parameter :: DFTOPT = '0:deg,1:dms'
character(24), parameter :: HGTOPT = '0:ellipsoidal,1:geodetic'
character(50), parameter :: GEOOPT = '0:internal,1:egm96,2:egm08_2.5,3:egm08_1,4:gsi2000'
character(14), parameter :: STAOPT = '0:all,1:single'
character(24), parameter :: STSOPT = '0:off,1:state,2:residual'
character(49), parameter :: ARMOPT = '0:off,1:continuous,2:instantaneous,3:fix-and-hold'
character(55), parameter :: POSOPT = '0:llh,1:xyz,2:single,3:posfile,4:rinexhead,5:rtcm,6:raw'
character(16), parameter :: TIDEOPT= '0:off,1:on,2:otl'
character(20), parameter :: PHWOPT = '0:off,1:on,2:precise'

type(prcopt_t), save :: prcopt_
type(solopt_t), save :: solopt_
character(1024),save :: exsats_
real*8, save :: elmask_

character(62), save :: helptext(13)=(/&
    "                                                              ",&
    " usage: spp [option]... file file [...]                       ",&
    "                                                              ",&
    " -?        print help                                         ",&
    " -s        output single solution                             ",&
    " -k file   input options from configuration file [off]        ",&
    " -o file   set output file [stdout]                           ",&
    " -ts ds ts start day/time (ds=y/m/d ts=h:m:s) [obs start time]",&
    " -te de te end day/time   (de=y/m/d te=h:m:s) [obs end time]  ",&
    " -ti tint  time interval  (sec) [all]                         ",&
    " -trop     trop option    (non saas) [saas]                   ",&
    " -elev     minimum elevation (deg) [10]                       ",&
    "                                                              "/)

contains
subroutine InitGlobal2()
implicit none
call InitGlobal()

! global variables in spp.f90
prcopt_=prcopt_default
solopt_=solopt_default
exsats_=''
elmask_=15.0

! global variables in rtklib.f90
solindex_=0
allsol_=sol_t(gtime_t(0,0.d0),0.d0,0.d0,0.d0,0,0,0)
sysopts=(/&
    opt_t('pos1-posmode',    3,  ivar(prcopt_%mode),       MODOPT ),&
    opt_t('pos1-frequency',  3,  ivar(prcopt_%nf),         FRQOPT ),&
    opt_t('pos1-soltype',    3,  ivar(prcopt_%soltype),    TYPOPT ),&
    opt_t('pos1-elmask',     1,  rvar(elmask_),            'deg'  ),&
    opt_t('pos1-ionoopt',    3,  ivar(prcopt_%ionoopt),    IONOPT ),&
    opt_t('pos1-tropopt',    3,  ivar(prcopt_%tropopt),    TRPOPT ),&
    opt_t('pos1-sateph',     3,  ivar(prcopt_%sateph),     EPHOPT ),&
    
    opt_t('pos1-posopt1',    3,  ivar(prcopt_%posopt(1)),  SWTOPT ),&
    opt_t('pos1-posopt2',    3,  ivar(prcopt_%posopt(2)),  SWTOPT ),&
    opt_t('pos1-posopt3',    3,  ivar(prcopt_%posopt(3)),  PHWOPT ),&
    opt_t('pos1-posopt4',    3,  ivar(prcopt_%posopt(4)),  SWTOPT ),&
    opt_t('pos1-posopt5',    3,  ivar(prcopt_%posopt(5)),  SWTOPT ),&
    opt_t('pos1-posopt6',    3,  ivar(prcopt_%posopt(6)),  SWTOPT ),&
    
    opt_t('pos1-exclsats',   2,  cvar(exsats_),            'prn ...'),&
    opt_t('pos1-navsys',     0,  ivar(prcopt_%navsys),     NAVOPT ),&
    opt_t('pos2-rejionno',   1,  rvar(prcopt_%maxinno),    'm'    ),&
    opt_t('pos2-rejgdop',    1,  rvar(prcopt_%maxgdop),    ''     ),&
    opt_t('pos2-niter',      0,  ivar(prcopt_%niter),      ''     ),&
    
    opt_t('out-solformat',   3,  ivar(solopt_%posf),       SOLOPT ),&
    opt_t('out-outhead',     3,  ivar(solopt_%outhead),    SWTOPT ),&
    opt_t('out-outopt',      3,  ivar(solopt_%outopt),     SWTOPT ),&
    opt_t('out-outvel',      3,  ivar(solopt_%outvel),     SWTOPT ),&
    opt_t('out-timesys',     3,  ivar(solopt_%times),      TSYOPT ),&
    
    opt_t('out-timeform',    3,  ivar(solopt_%timef),      TFTOPT ),&
    opt_t('out-timendec',    0,  ivar(solopt_%timeu),      ''     ),&
    opt_t('out-degform',     3,  ivar(solopt_%degf),       DFTOPT ),&
    opt_t('out-fieldsep',    2,  cvar(solopt_%sep),        ''     ),&
    opt_t('out-maxsolstd',   1,  rvar(solopt_%maxsolstd),  'm'    ),&
    
    opt_t('out-height',      3,  ivar(solopt_%height),     HGTOPT ),&
    opt_t('out-geoid',       3,  ivar(solopt_%geoid),      GEOOPT ),&
    opt_t('out-solstatic',   3,  ivar(solopt_%solstatic),  STAOPT ),&
    opt_t('out-nmeaintv1',   1,  rvar(solopt_%nmeaintv(1)),'s'    ),&
    opt_t('out-nmeaintv2',   1,  rvar(solopt_%nmeaintv(2)),'s'    ),&
    
    opt_t('out-outstat',     3,  ivar(solopt_%sstat),      STSOPT ),&
    opt_t('misc-rnxopt1',    2,  cvar(prcopt_%rnxopt(1)),  ''     ),&
    opt_t('misc-rnxopt2',    2,  cvar(prcopt_%rnxopt(2)),  ''     ),&
    opt_t('out-issingle',    3,  ivar(solopt_%issingle),   SWTOPT ),&
    opt_t('',0,vvar(),'') &
/)
end subroutine

subroutine printhelp()
implicit none
integer*4 i
do i=1,size(helptext,dim=1)
    write(*,*) trim(helptext(i))
enddo
end subroutine

subroutine ApplySysopts()
implicit none
prcopt_%mode =        sysopts( 1)%var%ivalue
prcopt_%nf =          sysopts( 2)%var%ivalue
prcopt_%soltype =     sysopts( 3)%var%ivalue
elmask_ =             sysopts( 4)%var%rvalue
prcopt_%ionoopt =     sysopts( 5)%var%ivalue
prcopt_%tropopt =     sysopts( 6)%var%ivalue
prcopt_%sateph =      sysopts( 7)%var%ivalue
prcopt_%posopt(1) =   sysopts( 8)%var%ivalue
prcopt_%posopt(2) =   sysopts( 9)%var%ivalue
prcopt_%posopt(3) =   sysopts(10)%var%ivalue
prcopt_%posopt(4) =   sysopts(11)%var%ivalue
prcopt_%posopt(5) =   sysopts(12)%var%ivalue
prcopt_%posopt(6) =   sysopts(13)%var%ivalue
exsats_ =             sysopts(14)%var%cvalue
prcopt_%navsys =      sysopts(15)%var%ivalue
prcopt_%maxinno =     sysopts(16)%var%rvalue
prcopt_%maxgdop =     sysopts(17)%var%rvalue
prcopt_%niter =       sysopts(18)%var%ivalue
solopt_%posf =        sysopts(19)%var%ivalue
solopt_%outhead =     sysopts(20)%var%ivalue
solopt_%outopt =      sysopts(21)%var%ivalue
solopt_%outvel =      sysopts(22)%var%ivalue
solopt_%times =       sysopts(23)%var%ivalue
solopt_%timef =       sysopts(24)%var%ivalue
solopt_%timeu =       sysopts(25)%var%ivalue
solopt_%degf =        sysopts(26)%var%ivalue
solopt_%sep =         sysopts(27)%var%cvalue
solopt_%maxsolstd =   sysopts(28)%var%rvalue
solopt_%height =      sysopts(29)%var%ivalue
solopt_%geoid =       sysopts(30)%var%ivalue
solopt_%solstatic =   sysopts(31)%var%ivalue
solopt_%nmeaintv(1) = sysopts(32)%var%rvalue
solopt_%nmeaintv(2) = sysopts(33)%var%rvalue
solopt_%sstat =       sysopts(34)%var%ivalue
prcopt_%rnxopt(1) =   sysopts(35)%var%cvalue
prcopt_%rnxopt(2) =   sysopts(36)%var%cvalue
solopt_%issingle =    sysopts(37)%var%ivalue
end subroutine

! system options buffer to options ------------------------------------------
subroutine SimplifyStr(str, splitter)
implicit none
character(*), intent(inout) :: str
character(*), intent(in) :: splitter
character(1024) split1, split2
integer*4 lensplit
split1=''; split2=''
lensplit=0
if(splitter=='')then
    str=trim(adjustl(str))
    do while(index(trim(str),'  ')/=0)
        str(index(trim(str),'  '):)=str(index(trim(str),'  ')+1:)
    enddo
else
    str=trim(adjustl(str))
    split1=trim(adjustl(splitter))
    split2=trim(split1)//trim(split1)
    lensplit=len_trim(split1)
    do while(index(str,trim(split2))/=0)
        str(index(str,trim(split2)):)=str(index(str,trim(split2))+lensplit:)
    enddo
endif
end subroutine

subroutine StringSplit(source, splitter, nsize, buff)
character(*), intent(in) :: source, splitter
integer*4, intent(out) :: nsize
character(*), intent(out) :: buff(:)
character(1024) str, split1
buff=''; nsize=0

if(splitter=='')then
    str=trim(source)
    call SimplifyStr(str,'')
    if(len_trim(str)==0) return
    nsize=1
    buff(1)=trim(str)
    do while(index(trim(str),' ')/=0)
        nsize=nsize+1
        buff(nsize-1)=str(1:index(str,' ')-1)
        buff(nsize)  =trim(str(index(str,' ')+1:))
        str=trim(buff(nsize))
    enddo
else
    str=trim(source)
    call SimplifyStr(str,splitter)
    split1=trim(adjustl(splitter))
    
    if(len_trim(str)<=len_trim(split1))then
        if(str==split1) return
        nsize=1; buff(1)=trim(str); return
    endif
    
    do while(index(str(len_trim(str)-len_trim(split1)+1:),trim(split1))==1)
        str=str(1:len_trim(str)-len_trim(splitter))
    enddo
    do while(index(str,trim(split1))==1)
        str=str(len_trim(split1)+1:)
        str=adjustl(str)
    enddo
    nsize=1; buff(1)=trim(str)
    do while(index(str,trim(split1))/=0)
        nsize=nsize+1
        buff(nsize-1)=str(1:index(str,trim(split1))-1)
        buff(nsize)  =trim(str(index(str,trim(split1))+len_trim(split1):))
        str=trim(buff(nsize))
    enddo
endif
end subroutine

subroutine buff2sysopts()
implicit none
character(1024) buff, p, id
character(100) buff2(100)
integer*4 i,sat,nsize
prcopt_%elmin=elmask_*D2R

! excluded satellites 
prcopt_%exsats=0
if(len_trim(exsats_)/=0)then
    buff=exsats_
    call StringSplit(buff,' ',nsize,buff2)
    do i=1,nsize
        p=buff2(i)
        if(p(1:1)=='+')then
            id=p(2:)
        else
            id=p(1:)
        endif
        sat=satid2no(id)
        if(sat==0) cycle
        prcopt_%exsats(sat)=icond(p(1:1)=='+',2,1)
    enddo
endif
! number of frequency (4:L1+L5) 
if (prcopt_%nf==4) prcopt_%nf=3
end subroutine

subroutine getsysopts(popt, sopt)
implicit none
type(prcopt_t), intent(out) :: popt
type(solopt_t), intent(out) :: sopt

call buff2sysopts()
popt=prcopt_
sopt=solopt_
end subroutine

type(opt_t) function searchopt(name, opts, optindex)
implicit none
character(*), intent(in) :: name
type(opt_t), intent(in) :: opts(:)
integer*4, intent(out) :: optindex  ! value changed
integer*4 i
i=1
do while(len_trim(opts(i)%name)/=0)
    if(index(opts(i)%name,trim(name))/=0)then
        searchopt=opts(i); optindex=i; return
    endif
    i=i+1
enddo
searchopt=opt_t('',0,vvar(),'')
end function

subroutine str2enum(str, comment, val, stat)
implicit none
character(*), intent(in) :: str, comment
integer*4, intent(out) :: val, stat
character(1024) p,p1
character(32) s
integer*4 len1, i

p=trim(comment)
110 if(index(p,trim(str))==0) goto 100
    if(index(p,trim(str))-1<1) goto 100
if(p(index(p,trim(str))-1:index(p,trim(str))-1)/=':')then
    p=p(index(p,trim(str))+1:); goto 110
endif
p1='_'//p(1:index(p,trim(str))-1)
len1=len_trim(p1)
do i=len1-1,1,-1
    if(('0'<=p1(i:i).and.p1(i:i)<='9').or.p1(i:i)==' ') cycle
    read(p1(i+1:len1-1),*,iostat=stat) val
    if(stat==0)then
        stat=1
    else
        stat=0
    endif
    return
enddo

100 write(s,"(A30':')") str
if(index(comment,trim(s))/=0)then  ! number 
    p=comment(index(comment,trim(s)):)
    read(p,*,iostat=stat) val
    if(stat==0)then
        stat=1
    else
        stat=0
    endif
    return
endif
stat=0
end subroutine

subroutine str2opt(opt, str, stat)
implicit none
type(opt_t), intent(out) :: opt
character(*), intent(in) :: str
integer*4, intent(out) :: stat
integer*4 itmp, info
real*8 rtmp
select case (opt%myformat)
case(0)
    read(str,*,iostat=info) itmp
    opt%var=ivar(itmp)
case(1)
    read(str,*,iostat=info) rtmp
    opt%var=rvar(rtmp)
case(2)
    opt%var=cvar(str)
case(3)
    call str2enum(str,opt%comment,itmp,stat)
    opt%var=ivar(itmp); return
case default
    stat=0; return
end select
stat=1
end subroutine

! discard space characters at tail ------------------------------------------
subroutine chop(str)
implicit none
character(*), intent(out) :: str

if(index(str,'#')/=0)then
    str=str(1:index(str,'#')-1)
else
    str=str
endif
end subroutine

integer*4 function loadopts(filepath, opts)
implicit none
character(*), intent(in) :: filepath
type(opt_t), intent(in) :: opts(:)
integer*4 :: fp=FPOPT, n, info, optindex
type(opt_t) opt
character(2048) buff, p
n=0
open(unit=fp,file=filepath,status='old',iostat=info)
if(info/=0)then
    loadopts=0; return
endif
do while(.true.)
    read(fp,"(A)",iostat=info) buff
    if(info/=0) exit
    n=n+1
    call chop(buff)
    if (len_trim(buff)==0) cycle

    if(index(buff,'=')==0)then
        write(unit=6,fmt="('invalid option ',A,' (',A,':',I4,')',/)") trim(buff),trim(filepath),n
        cycle
    endif
    p=buff(index(buff,'=')+1:)
    buff=buff(1:index(buff,'=')-1)
    call chop(buff)
    opt=searchopt(buff,opts,optindex)
    if (len_trim(opt%name)==0) cycle
    call str2opt(opt,p,info)
    if(info==0)then
        write(unit=6,fmt="('invalid option value ',A,' (',A,':',I4,')',/)") trim(buff),trim(filepath),n
        cycle
    else
        sysopts(optindex)=opt
    endif
enddo
close(unit=fp,status='keep')
loadopts=1
end function

! rnx2rtkp ---------------------------------------------------
subroutine rnx2rtkp(argcIn, argvIn)
implicit none
integer*4, intent(in) :: argcIn
character(*), intent(in) :: argvIn(:)
type(prcopt_t) prcopt
type(solopt_t) solopt
type(gtime_t) ts, te
real*8 :: tint, es(6),ee(6)
integer*4 :: i, n, ret
character(1024) :: infile(MAXFILE), outfile, buff(6)
integer*4 argc, info, nsize
character(1024), pointer :: argv(:)

call InitGlobal2()
prcopt=prcopt_default; solopt=solopt_default
ts=gtime_t(0,0.d0); te=gtime_t(0,0.d0)
tint=0.d0; es=(/2000,1,1,0,0,0/); ee=(/2000,12,31,23,59,59/)
i=0; n=0; ret=0
infile(MAXFILE)=''; outfile=''
exsats_=''; elmask_=0.d0

argc=argcIn
allocate(argv(argc))
argv(1:argc)=argvIn(1:argc)
!write(unit=6,fmt="('argc =',I2)") argc
!do i=1,argc
!    write(unit=6,fmt="(I2,' ',A)") i, argv(i)
!enddo
!print*

if(argc==0)then
    call printhelp(); call exit(0)
endif

do i=1,argc
    if(argv(i)=='-?'.or.argv(i)=='-h')then
        call printhelp(); call exit(0)
    endif
enddo

do i=1,argc
    if(argv(i)=='-k'.and.(i+1)<=argc)then
        if (loadopts(argv(i+1),sysopts)==0)then
            deallocate(argv); call exit(3)
        else
            call ApplySysopts()
        endif
        call getsysopts(prcopt,solopt)
    endif
enddo

n=0; i=1
do while(i<=argc)
    if(argv(i)=='-o'.and.(i+1)<=argc)then
        outfile=argv(i+1); i=i+1
    elseif(argv(i)=='-ts'.and.(i+2)<=argc)then
        call StringSplit(argv(i+1),'/',nsize,buff(1:3))
        call StringSplit(argv(i+2),':',nsize,buff(4:6))
        read(buff(1:3),*,iostat=info) es(1),es(2),es(3)
        read(buff(4:6),*,iostat=info) es(4),es(5),es(6); i=i+2
        ts=epoch2time(es)
    elseif(argv(i)=='-te'.and.(i+2)<=argc)then
        call StringSplit(argv(i+1),'/',nsize,buff(1:3))
        call StringSplit(argv(i+2),':',nsize,buff(4:6))
        read(buff(1:3),*,iostat=info) ee(1),ee(2),ee(3)
        read(buff(4:6),*,iostat=info) ee(4),ee(5),ee(6); i=i+2
        te=epoch2time(ee)
    elseif(argv(i)=='-ti'.and.(i+1)<=argc)then
        read(argv(i+1),*,iostat=info) tint; i=i+1
    elseif(argv(i)=='-k'.and.(i+1)<=argc)then
        i=i+2; cycle
    elseif(argv(i)=='-s')then
        solopt%issingle=1
    elseif(argv(i)=='-trop'.and.(i+1)<=argc)then
        if(argv(i+1)=='non')  prcopt%tropopt=0
        if(argv(i+1)=='saas') prcopt%tropopt=1
        i=i+1
    elseif(argv(i)=='-elev'.and.(i+1)<=argc)then
        read(argv(i+1),*,iostat=info) prcopt%elmin
        prcopt%elmin=prcopt%elmin*D2R; i=i+1
    elseif(n<MAXFILE)then
        n=n+1; infile(n)=argv(i)
    endif
    i=i+1
enddo

if(n<=0)then
    write(*,*) 'error : no input file'
    deallocate(argv); call exit(2)
endif

ret=postpos(ts,te,tint,0.d0,prcopt,solopt,infile,n,outfile,'','')  ! 0-error, 1-right
deallocate(argv);
if(ret==0)then
    write(unit=6,fmt="(A40)") ''; call exit(1)
endif
if(ret==1) call exit(0)
end subroutine
end module spp_all_

! main ---------------------------------------------------
program spp
use spp_all_
implicit none
integer*4 argc, i
character(1024), pointer :: argv(:)
argc=iargc()
allocate(argv(argc))
do i=1,argc
    call getarg(i,argv(i))
enddo
call rnx2rtkp(argc,argv)
deallocate(argv)
end program spp

!! main ---------------------------------------------------
!program spp
!use spp_all_
!implicit none
!integer*4 argc, calnum, info, i, nsize
!character(1024), pointer :: argv(:)
!character(1024) argv34(1024), buff, buffsplit(100), filepath
!
!!argc=iargc()
!!allocate(argv(argc))
!!do i=1,argc
!!    call getarg(i,argv(i))
!!enddo
!filepath='d:/desktop/data/'
!!filepath='d:\\desktop\\data\\'
!
!argc=4  !8
!allocate(argv(argc))
!!argv(5)='-k'
!!argv(6)=trim(filepath)//'rtkoption.conf'
!argv(1)='-o'
!!argv(7)='-ti'
!!argv(8)='300'
!
!! -------- 
!calnum=0
!open(unit=17,file=trim(filepath)//'path.txt',status='old',iostat=info)
!do while(.true.)
!    read(17,"(a)",iostat=info) buff
!    if(info/=0) exit
!    call simplifystr(buff,'')
!    if(buff(1:1)=='#'.or.buff=='')then
!        cycle
!    else
!        call stringsplit(buff,'',nsize,buffsplit)
!        if(nsize>=2)then
!            calnum=calnum+2
!            argv34(calnum-1)=trim(filepath)//trim(buffsplit(1))
!            argv34(calnum)  =trim(filepath)//trim(buffsplit(2))
!        endif
!    endif
!enddo
!close(unit=17,status='keep')
!
!do i=1,calnum/2
!    argv(3)=argv34(i*2-1)
!    argv(4)=argv34(i*2)
!    argv(2)=trim(filepath)//'outpos_'//argv(4)(len_trim(argv(4))-11:len_trim(argv(4))-4)//'_f.pos'
!    call rnx2rtkp(argc,argv)
!enddo
!deallocate(argv)
!pause
!end program spp
