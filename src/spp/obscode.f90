
! Obs code to obs code string -----------------------------------------------
!* * convert obs code to obs code string
!* * args   : character(1) code I obs code (CODE_???)
!* *          integer*4 *freq   IO     frequency (NULL: no output)
!* * (1:L1/E1, 2:L2/B1, 3:L5/E5a/L3, 4:L6/LEX/B3,
!* *  5:E5b/B2, 6:E5(a+b), 7:S)
!* * return : obs code string ('1C','1P','1P',...)
!* * notes  : obs codes are based on reference (6) and qzss extension
!* *-------------------------------------------------------------------------
subroutine code2obs(code, freq, obs)
implicit none
include 'file_para.h'
integer*4,intent(in) :: code
integer*4,intent(out) :: freq
character(*),intent(out) :: obs
freq=0
if (code<=CODE_NONE .or. MAXCODE<code)then
    obs=''; return
endif
freq=obsfreqs(code+1)
obs=obscodes(code+1)
end subroutine

! obs type string to obs code -----------------------------------------------
!* * convert obs code type string to obs code
!* * args   : char   *str      I      obs code string ('1C','1P','1Y',...)
!* *          integer*4 *freq  IO     frequency (1:L1,2:L2,3:L5,4:L6,5:L7,6:L8,0:err)
!* *                               (NULL: no output)
!* * return : obs code (CODE_???)
!* * notes  : obs codes are based on reference (6) and qzss extension
!* *-------------------------------------------------------------------------
subroutine obs2code(obs, freq, code)
implicit none
include 'file_para.h'
character(*),intent(in) :: obs
integer*4,intent(out) :: freq, code
integer*4 i
code=CODE_NONE; freq=0
do i=2,size(obscodes,1)
    if(len_trim(obscodes(i))==0) exit
    if(obs/=obscodes(i)) cycle
    freq=obsfreqs(i)
    code=i-1; return
enddo
code=CODE_NONE
end subroutine

! convert rinex obs type ver.2 % ver.3 --------------------------------------
subroutine convcode(ver, sys, str, mytype)
implicit none
include 'file_para.h'
real*8, intent(in) :: ver
integer*4, intent(in) :: sys
character(*), intent(in) :: str
character(*), intent(out) :: mytype
mytype=""
if(str=="P1")then
    if(sys==SYS_GPS) mytype="C1W"
    if(sys==SYS_GLO) mytype="C1P"
elseif(str=="P2")then
    if(sys==SYS_GPS) mytype="C2W"
    if(sys==SYS_GLO) mytype="C2P"
elseif(str=="C1")then
    if (ver<2.12d0.and.sys==SYS_GPS) mytype="C1C"
    if (ver<2.12d0.and.sys==SYS_GLO) mytype="C1C"
    if (ver<2.12d0.and.sys==SYS_GAL) mytype="C1X"  ! ver.2.12 
    if (ver<2.12d0.and.sys==SYS_QZS) mytype="C1C"
    if (ver<2.12d0.and.sys==SYS_SBS) mytype="C1C"
elseif(str=="C2")then
    if(sys==SYS_GPS)then
        if(ver>=2.12d0)then
            mytype="C2W"
        else
            mytype="C2X"
        endif
    elseif(sys==SYS_GLO)then
        mytype="C2C"
    elseif(sys==SYS_QZS)then
        mytype="C2X"
    elseif(sys==SYS_CMP)then
        mytype="C1X"
    endif
elseif(ver>=2.12d0.and.str(2:2)=='A')then
    if(sys==SYS_GPS) mytype=str(1:1)//"1C"
    if(sys==SYS_GLO) mytype=str(1:1)//"1C"
    if(sys==SYS_QZS) mytype=str(1:1)//"1C"
    if(sys==SYS_SBS) mytype=str(1:1)//"1C"
elseif(ver>=2.12d0.and.str(2:2)=='B')then
    if(sys==SYS_GPS) mytype=str(1:1)//"1X"
    if(sys==SYS_QZS) mytype=str(1:1)//"1X"
elseif(ver>=2.12d0.and.str(2:2)=='C')then
    if(sys==SYS_GPS) mytype=str(1:1)//"2X"
    if(sys==SYS_QZS) mytype=str(1:1)//"2X"
elseif(ver>=2.12d0.and.str(2:2)=='D')then
    if(sys==SYS_GLO) mytype=str(1:1)//"2C"
elseif(ver>=2.12d0.and.str(2:2)=='1')then
    if(sys==SYS_GPS) mytype=str(1:1)//"1W"
    if(sys==SYS_GLO) mytype=str(1:1)//"1P"
    if(sys==SYS_GAL) mytype=str(1:1)//"1X"
    if(sys==SYS_CMP) mytype=str(1:1)//"1X"
elseif(ver< 2.12d0.and.str(2:2)=='1')then
    if(sys==SYS_GPS) mytype=str(1:1)//"1C"
    if(sys==SYS_GLO) mytype=str(1:1)//"1C"
    if(sys==SYS_GAL) mytype=str(1:1)//"1X"
    if(sys==SYS_QZS) mytype=str(1:1)//"1C"
    if(sys==SYS_SBS) mytype=str(1:1)//"1C"
elseif(str(2:2)=='2')then
    if(sys==SYS_GPS) mytype=str(1:1)//"2W"
    if(sys==SYS_GLO) mytype=str(1:1)//"2P"
    if(sys==SYS_QZS) mytype=str(1:1)//"2X"
    if(sys==SYS_CMP) mytype=str(1:1)//"1X"
elseif(str(2:2)=='5')then
    if(sys==SYS_GPS) mytype=str(1:1)//"5X"
    if(sys==SYS_GAL) mytype=str(1:1)//"5X"
    if(sys==SYS_QZS) mytype=str(1:1)//"5X"
    if(sys==SYS_SBS) mytype=str(1:1)//"5X"
elseif(str(2:2)=='6')then
    if(sys==SYS_GAL) mytype=str(1:1)//"6X"
    if(sys==SYS_QZS) mytype=str(1:1)//"6X"
    if(sys==SYS_CMP) mytype=str(1:1)//"6X"
elseif(str(2:2)=='7')then
    if(sys==SYS_GAL) mytype=str(1:1)//"7X"
    if(sys==SYS_CMP) mytype=str(1:1)//"7X"
elseif(str(2:2)=='8')then
    if(sys==SYS_GAL) mytype=str(1:1)//"8X"
endif
end subroutine

! set code priority ---------------------------------------------------------
! set code priority for multiple codes in a frequency
! args   : integer*4 sys  I     system (or of SYS_???)
!          integer*4 freq I     frequency (1:L1,2:L2,3:L5,4:L6,5:L7,6:L8,7:L9)
!          char   *pri    I     priority of codes (series of code characters)
!                               (higher priority precedes lower)
! return : none
!----------------------------------------------------------------------------
!subroutine setcodepri(sys, freq, pri)
!implicit none
!include 'file_para.h'
!integer*4, intent(in) :: sys, freq
!character(*), intent(in) :: pri
!if (freq<=0 .or. MAXFREQ<freq) return
!if (and(sys,SYS_GPS)/=0) codepris(1,freq)=pri
!if (and(sys,SYS_GLO)/=0) codepris(2,freq)=pri
!if (and(sys,SYS_GAL)/=0) codepris(3,freq)=pri
!if (and(sys,SYS_QZS)/=0) codepris(4,freq)=pri
!if (and(sys,SYS_SBS)/=0) codepris(5,freq)=pri
!if (and(sys,SYS_CMP)/=0) codepris(6,freq)=pri
!if (and(sys,SYS_IRN)/=0) codepris(7,freq)=pri
!end subroutine

! get code priority ---------------------------------------------------------
! get code priority for multiple codes in a frequency
! args   : integer*4 sys     I     system (SYS_???)
!          character(1) code I     obs code (CODE_???)
!          char   *opt       I     code options (NULL:no option)
! return : priority (15:highest-1:lowest,0:error)
!----------------------------------------------------------------------------
integer*4 function getcodepri(sys, code, opt)
implicit none
include 'file_para.h'
integer*4, intent(in) :: sys, code
character(*), intent(in) :: opt
character(256) p, optstr
character(2) obs
character(8) str
integer*4 i, j, info, posi
integer*4, external :: icond
select case(sys)
case(SYS_GPS)
    i=0; optstr="('-GL'A2)"
case(SYS_GLO)
    i=1; optstr="('-RL'A2)"
case(SYS_GAL)
    i=2; optstr="('-EL'A2)"
case(SYS_QZS)
    i=3; optstr="('-JL'A2)"
case(SYS_SBS)
    i=4; optstr="('-SL'A2)"
case(SYS_CMP)
    i=5; optstr="('-CL'A2)"
case(SYS_IRN)
    i=6; optstr="('-IL'A2)"
case default
    getcodepri=0; return
end select
call code2obs(code, j, obs)
! parse code options
p=trim(opt)
do while(len_trim(p)/=0 .and. index(p,'-')/=0)
    p=p(index(p,'-'):)
    read(p,trim(optstr),iostat=info) str
    if(info/=0 .or. str(1:1)/=obs(1:1))then
        p=p(2:); cycle
    endif
    getcodepri=icond(str(2:2)==obs(2:2),15,0); return
enddo
! search code priority 
! return (p=strchr(codepris(i)(j-1),obs(1)))?14-(integer*4)(p-codepris(i)(j-1)):0
if(j==0)then
    getcodepri=14; return
endif
posi=index(codepris(i+1,j),obs(2:2))
getcodepri=icond(posi/=0,14-(posi-1),0)  ! to be verified
end function
