!! reference and acknowledgement to RTKLIB
!! modified by Kunlun Zhang
! constants/macros ----------------------------------------------------------

module rinex_f90_
use rtkcmn_f90_
implicit none
integer*4, parameter :: NUMSYS      = 7                   ! number of systems 
integer*4, parameter :: MAXRNXLEN   = (16*MAXOBSTYPE+4)   ! max rinex record length 
integer*4, parameter :: MAXPOSHEAD  = 1024                ! max head line position 
integer*4, parameter :: MINFREQ_GLO = -7                  ! min frequency number glonass 
integer*4, parameter :: MAXFREQ_GLO = 13                  ! max frequency number glonass 
integer*4, parameter :: NINCOBS     = 262144/2            ! inclimental number of obs data 

private navsys, syscodes, obscodes, frqcodes
integer*4, parameter :: navsys(8)=(/SYS_GPS,&
    SYS_GLO,SYS_GAL,SYS_QZS,SYS_SBS,SYS_CMP,SYS_IRN,0/)   ! satellite systems 
character(7), parameter :: syscodes='GREJSCI'             ! satellite system codes 
character(4), parameter :: obscodes='CLDS'                ! obs type codes 
character(7), parameter :: frqcodes='1256789'             ! frequency codes 
real*8, parameter :: ura_eph(16)=(/&                      ! ura values (ref [3] 20.3.3.3.1.1) 
    2.4,3.4,4.85,6.85,9.65,13.65,24.0,48.0,96.0,192.0,384.0,768.0,1536.0,&
    3072.0,6144.0,0.0/)
real*8, parameter :: ura_nominal(16)=(/&                  ! ura nominal values 
    2.0,2.8,4.0,5.7,8.0,11.3,16.0,32.0,64.0,128.0,256.0,512.0,1024.0,&
    2048.0,4096.0,8192.0/)

!! type definition -----------------------------------------------------------
type sigind_t                                             ! signal index type
    integer*4 n                                           ! number of index 
    integer*4 frq(MAXOBSTYPE)                             ! signal frequency (1:L1,2:L2,...) 
    integer*4 pos(MAXOBSTYPE)                             ! signal index in obs data (-1:no) 
    integer*4 pri(MAXOBSTYPE)                             ! signal priority (15-0) 
    integer*4 mytype(MAXOBSTYPE)                          ! type (0:C,1:L,2:D,3:S) 
    integer*4 code(MAXOBSTYPE)                            ! obs code (CODE_L??) 
    real*8 shift(MAXOBSTYPE)                              ! phase shift (cycle) 
end type

contains
! set string without tail space ---------------------------------------------
subroutine setstr(dst, src, n)
implicit none
character(*), intent(out) :: dst
character(*), intent(in) :: src
integer*4, intent(in) :: n
if(n<1.or.n>len_trim(src))then
    dst=''; return
endif
dst=trim(src(1:n))
end subroutine

! adjust time considering week handover -------------------------------------
type(gtime_t) function adjweek(t, t0)
implicit none
type(gtime_t), intent(in) :: t,t0
real*8 tt
tt=timediff(t,t0)
if (tt<-302400.d0)then
    adjweek=timeadd(t, 604800.d0); return
endif
if (tt> 302400.d0)then
    adjweek=timeadd(t,-604800.d0); return
endif
adjweek=t
end function

! adjust time considering week handover -------------------------------------
type(gtime_t) function adjday(t, t0)
implicit none
type(gtime_t), intent(in) :: t,t0
real*8 tt
tt=timediff(t,t0)
if (tt<-43200.d0)then
    adjday=timeadd(t, 86400.d0); return
endif
if (tt> 43200.d0)then
    adjday=timeadd(t,-86400.d0); return
endif
adjday=t
end function

! ura index to ura nominal value (m) ----------------------------------------
real*8 function uravalue(sva)
implicit none
integer*4, intent(in) :: sva
uravalue=rcond(0<=sva.and.sva<15,ura_nominal(sva+1),8192.d0)
end function

! ura value (m) to ura index ------------------------------------------------
integer*4 function uraindex(value)
implicit none
real*8, intent(in) :: value
integer*4 i
do i=1,15
    if(ura_eph(i)>=value) exit
enddo
uraindex=i-1
end function

! galileo sisa index to sisa nominal value (m) ------------------------------
real*8 function sisa_value(sisa)
implicit none
integer*4, intent(in) :: sisa
if     (sisa<= 49)then
    sisa_value=sisa*0.01d0
elseif (sisa<= 74)then
    sisa_value=0.5d0+(sisa-50)*0.02d0
elseif (sisa<= 99)then
    sisa_value=1.d0+ (sisa-75)*0.04d0
elseif (sisa<=125)then
    sisa_value=2.d0+(sisa-100)*0.16d0
else
    sisa_value=-1.d0  ! unknown or NAPA 
endif
end function

! galileo sisa value (m) to sisa index --------------------------------------
integer*4 function sisa_index(value)
implicit none
real*8, intent(in) :: value
if (value<0.d0  .or.  value>6.d0)then  ! unknown or NAPA 
    sisa_index=255
elseif (value<=0.5d0)then
    sisa_index=int(value/0.01d0)
elseif (value<=1.d0)then
    sisa_index=int((value-0.5d0)/0.02d0)+50
elseif (value<=2.d0)then
    sisa_index=int((value-1.d0)/0.04d0)+75
else
    sisa_index=int(int(value-2.d0)/0.16d0)+100
endif
end function

! initialize station parameter ----------------------------------------------
subroutine init_sta(sta)
implicit none
type(sta_t), intent(out) :: sta
sta%name   =''
sta%marker =''
sta%antdes =''
sta%antsno =''
sta%rectype=''
sta%recver =''
sta%recsno =''
sta%antsetup = 0
sta%itrf     = 0
sta%deltype  = 0
sta%pos=0.d0
sta%del=0.d0
sta%hgt=0.d0
end subroutine

! decode fname from fpath ---------------------------------------------------
subroutine decodefname(fpath, fname)
implicit none
character(*), intent(in) :: fpath
character(*), intent(out) :: fname
integer*4 ind1, ind2
fname=fpath
ind1=index(fpath,'/',BACK = .TRUE.)
ind2=index(fpath,'\\',BACK = .TRUE.)
if(ind1/=0)then
    fname=trim(fpath(ind1+1:)); return
endif
if(ind2/=0)then
    fname=trim(fpath(ind2+2:)); return
endif
end subroutine

! decode marker name from filepath ------------------------------------------
subroutine decodemarker(fpath, marker)
implicit none
character(*), intent(in) :: fpath
character(*), intent(out) :: marker
integer*4 flen
character(128) fname
marker=''
flen=len_trim(fpath)
if(flen==0) return
if(fpath(flen:flen)=='o'.or.fpath(flen:flen)=='O')then
    call decodefname(fpath,fname)
    call setstr(marker,fname,min(20,len_trim(fname))); return
endif
end subroutine

! decode nav header ---------------------------------------------------------
subroutine decode_navh(buff, nav)
implicit none
character(*), intent(in) :: buff
type(nav_t), intent(out) :: nav
integer*4 i,j,info
character(20) :: label
label=buff(61:)

if(index(label,'ION ALPHA')/=0)then
    read(buff(3:50),"(4D12.4)",iostat=info) nav%ion_gps(1:4)
elseif(index(label,'ION BETA')/=0)then
    read(buff(3:50),"(4D12.4)",iostat=info) nav%ion_gps(5:8)
elseif(index(label,'DELTA-UTC: A0,A1,T,W')/=0)then  ! opt ver.2 
    read(buff(4:59),"(2D19.12,2D9.0)",iostat=info) nav%utc_gps(1:4)
elseif(index(label,'IONOSPHERIC CORR')/=0)then
    if(buff(1:4)=='GPSA') read(buff(6:53),"(4D12.4)",iostat=info) nav%ion_gps(1:4)
    if(buff(1:4)=='GPSB') read(buff(6:53),"(4D12.4)",iostat=info) nav%ion_gps(5:8)
    if(buff(1:4)=='GAL ') read(buff(6:53),"(4D12.4)",iostat=info) nav%ion_gal(1:4)
    if(buff(1:4)=='BDSA') read(buff(6:53),"(4D12.4)",iostat=info) nav%ion_cmp(1:4)
    if(buff(1:4)=='BDSB') read(buff(6:53),"(4D12.4)",iostat=info) nav%ion_cmp(5:8)
elseif(index(label,'TIME SYSTEM CORR')/=0)then      ! opt ver.3 
    !if(buff(1:4)=='GPUT') read(buff(5:50),"(D18.10,D16.9,I7,I5)",iostat=info) nav%utc_gps(1:4)
    !if(buff(1:4)=='GLUT') read(buff(5:50),"(D18.10,D16.9,I7,I5)",iostat=info) nav%utc_glo(1:4)
    if(buff(1:4)=='GPUT') read(buff(5:50),*,iostat=info) nav%utc_gps(1:4)
    if(buff(1:4)=='GLUT') read(buff(5:50),*,iostat=info) nav%utc_glo(1:4)
    if(buff(1:4)=='GAUT') read(buff(5:50),*,iostat=info) nav%utc_gal(1:4)  ! v.3.02 
    if(buff(1:4)=='QZUT') read(buff(5:50),*,iostat=info) nav%utc_qzs(1:4)
    if(buff(1:4)=='BDUT') read(buff(5:50),*,iostat=info) nav%utc_cmp(1:4)
    if(buff(1:4)=='SBUT') read(buff(5:50),*,iostat=info) nav%utc_sbs(1:4)
    if(buff(1:4)=='IRUT') read(buff(5:50),*,iostat=info) nav%utc_irn(1:4)
elseif(index(label,'LEAP SECONDS')/=0)then
    read(buff(1:6),"(I6)",iostat=info) nav%leaps
endif
end subroutine

! decode gnav header --------------------------------------------------------
subroutine decode_gnavh(buff, nav)
implicit none
character(*), intent(in) :: buff
type(nav_t), intent(out) :: nav
character(20) :: label
integer*4 info
label=buff(61:)
!if(index(label,'CORR TO SYTEM TIME')/=0) ;
if(index(label,'LEAP SECONDS')/=0)then
    read(buff(1:6),"(I6)",iostat=info) nav%leaps
endif
end subroutine

! convert rinex obs type ver.2 % ver.3 -------------------------------------
subroutine convcode(ver, sys, str, mytype)
implicit none
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

! decode obs header ---------------------------------------------------------
subroutine decode_obsh(fp, buff, ver, tsys, tobs, obs, nav, sta)
implicit none
integer*4, intent(in) :: fp
character(*), intent(inout) :: buff
real*8, intent(in) :: ver
integer*4, intent(out) :: tsys
character(3), intent(out) :: tobs(:,:)
type(obs_t), intent(out) :: obs
type(nav_t), intent(out) :: nav
type(sta_t), intent(out) :: sta

! default codes for unknown code 
character(8) :: defcodes(7)=(/&
    'CWX    ',&  ! GPS: L125____
    'CC     ',&  ! GLO: L12_____
    'X XXXX ',&  ! GAL: L1_5678_
    'CXXX   ',&  ! QZS: L1256___
    'C X    ',&  ! SBS: L1_5____
    'X  XX  ',&  ! BDS: L1__67__
    '  A   A'/)  ! IRN: L__5___9
real*8 del(3)
integer*4 i,j,k,n,nt,prn,fcn,info,posi
character(1) p
character(20) label
character(4) str
label=buff(61:)
if(index(label,'MARKER NAME')/=0)then
    call setstr(sta%name,buff,60)
elseif(index(label,'MARKER NUMBER')/=0)then  ! opt 
    call setstr(sta%marker,buff,20)
elseif(index(label,'MARKER TYPE')/=0)then
elseif(index(label,'OBSERVER / AGENCY')/=0)then
elseif(index(label,'REC # / TYPE / VERS')/=0)then
    call setstr(sta%recsno ,buff( 1:),20)
    call setstr(sta%rectype,buff(21:),20)
    call setstr(sta%recver ,buff(41:),20)
elseif(index(label,'ANT # / TYPE')/=0)then
    call setstr(sta%antsno,buff( 1:),20)
    call setstr(sta%antdes,buff(21:),20)
elseif(index(label,'APPROX POSITION XYZ')/=0)then
    sta%pos(1)=str2num(buff, 1,14)
    sta%pos(2)=str2num(buff,15,14)
    sta%pos(3)=str2num(buff,29,14)
elseif(index(label,'ANTENNA: DELTA H/E/N')/=0)then
    sta%del(3)=str2num(buff, 1,14)
    sta%del(1)=str2num(buff,15,14)
    sta%del(2)=str2num(buff,29,14)
elseif(index(label,'ANTENNA: DELTA X/Y/Z')/=0)then  ! opt ver.3 
elseif(index(label,'ANTENNA: PHASECENTER')/=0)then  ! opt ver.3 
elseif(index(label,'ANTENNA: B.SIGHT XYZ')/=0)then  ! opt ver.3 
elseif(index(label,'ANTENNA: ZERODIR AZI')/=0)then  ! opt ver.3 
elseif(index(label,'ANTENNA: ZERODIR XYZ')/=0)then  ! opt ver.3 
elseif(index(label,'CENTER OF MASS: XYZ' )/=0)then  ! opt ver.3 
elseif(index(label,'SYS / # / OBS TYPES' )/=0)then  ! ver.3 
    p=buff(1:1)
    if(index(syscodes,p)==0) return
    i=index(syscodes,p)
    tobs(i,:)=""
    read(buff(4:6),"(I3)",iostat=info) n
    j=1; nt=1; k=8
    do while(j<=n)
        if(k>56)then
            read(fp,"(A)",iostat=info) buff
            if(info/=0)exit
            k=8
        endif
        if(nt<MAXOBSTYPE) call setstr(tobs(i,nt),buff(k:),3)
        nt=nt+1; j=j+1; k=k+4
    enddo
    tobs(i,nt)='end'
    
    ! change beidou B1 code: 3.02 
    if (i==6 .and. dabs(ver-3.02d0)<1d-4)then
        do j=1,nt
            if(tobs(i,j)(2:2)=='2') tobs(i,j)(2:2)='1'
        enddo
    endif
    
    ! if unknown code in ver.3, set default code 
    do j=1,nt  !-1
        if(tobs(i,j)(3:3)/=' ') cycle
        if(index(frqcodes,tobs(i,j)(2:2))==0) cycle
        posi=index(frqcodes,tobs(i,j)(2:2))
        tobs(i,j)(3:3)=defcodes(i)(posi+1:posi+1)
    enddo
elseif(index(label,'WAVELENGTH FACT L1/2')/=0)then  ! opt ver.2 
elseif(index(label,'# / TYPES OF OBSERV' )/=0)then  ! ver.2
    read(buff(1:6),"(I6)",iostat=info) n
    i=1; nt=1; j=11
    tobs(1:2,:)=''
    do while(i<=n)
        if(j>59)then
            read(fp,"(A)",iostat=info) buff
            if(info/=0) exit
            j=11
        endif
        if (nt>=MAXOBSTYPE)then
            i=i+1; j=j+6; cycle
        endif
        if(ver<=2.99d0)then
            call setstr(str,buff(j:),2)
            call convcode(ver,SYS_GPS,str,tobs(1,nt))
            call convcode(ver,SYS_GLO,str,tobs(2,nt))
            call convcode(ver,SYS_GAL,str,tobs(3,nt))
            call convcode(ver,SYS_QZS,str,tobs(4,nt))
            call convcode(ver,SYS_SBS,str,tobs(5,nt))
            call convcode(ver,SYS_CMP,str,tobs(6,nt))
        endif
        nt=nt+1; i=i+1; j=j+6
    enddo
    tobs(1,nt)='end'
    tobs(2,nt)='end'
    tobs(3,nt)='end'
    tobs(4,nt)='end'
    tobs(5,nt)='end'
    tobs(6,nt)='end'
elseif(index(label,'SIGNAL STRENGTH UNIT')/=0)then
elseif(index(label,'INTERVAL'            )/=0)then
    read(buff(1:11),*,iostat=info) obs%tint
elseif(index(label,'TIME OF FIRST OBS'   )/=0)then
    if(buff(49:51)=="GPS") tsys=TSYS_GPS
    if(buff(49:51)=="GLO") tsys=TSYS_UTC
    if(buff(49:51)=="GAL") tsys=TSYS_GAL
    if(buff(49:51)=="QZS") tsys=TSYS_QZS  ! ver.3.02
    if(buff(49:51)=="BDT") tsys=TSYS_CMP  ! ver.3.02
    if(buff(49:51)=="IRN") tsys=TSYS_IRN  ! ver.3.03
elseif(index(label,'TIME OF LAST OBS'    )/=0)then  ! opt 
elseif(index(label,'RCV CLOCK OFFS APPL' )/=0)then  ! opt 
elseif(index(label,'SYS / DCBS APPLIED'  )/=0)then  ! opt ver.3 
elseif(index(label,'SYS / PCVS APPLIED'  )/=0)then  ! opt ver.3 
elseif(index(label,'SYS / SCALE FACTOR'  )/=0)then  ! opt ver.3 
elseif(index(label,'SYS / PHASE SHIFTS'  )/=0)then  ! ver.3.01 
elseif(index(label,'GLONASS SLOT / FRQ #')/=0)then  ! ver.3.02 
elseif(index(label,'GLONASS COD/PHS/BIS' )/=0)then  ! ver.3.02
elseif(index(label,'LEAP SECONDS'        )/=0)then  ! opt 
    nav%leaps=int(str2num(buff,1,6))
elseif(index(label,'# OF SALTELLITES'    )/=0)then  ! opt 
elseif(index(label,'PRN / # OF OBS'      )/=0)then  ! opt 
endif
end subroutine

! read rinex header ---------------------------------------------------------
subroutine readrnxh(fp, ver, mytype, sys, tsys, tobs, obs, nav, sta, stat)
implicit none
integer*4, intent(in) :: fp
real*8, intent(out) :: ver
character(*), intent(out) :: mytype
integer*4, intent(out) :: sys, tsys, stat
character(3), intent(out) :: tobs(:,:)
type(obs_t), intent(out) :: obs
type(nav_t), intent(out) :: nav
type(sta_t), intent(out) :: sta
real*8 bias
character(MAXRNXLEN) buff
character(20) label
integer*4 :: i, myblock, sat, info
i=0; myblock=0; sat=0; info=0
label=buff(61:)
ver=2.10; mytype=''; sys=SYS_GPS

do while(.true.)
    read(fp,"(A)",iostat=info) buff
    label=buff(61:)
    if(info/=0) exit
    if (len_trim(buff)<=60)then
        cycle
    elseif(index(label,'RINEX VERSION / TYPE')/=0)then
        read(buff(1:9),*,iostat=info) ver
        read(buff(21:40),"(A)",iostat=info) mytype
        ! satellite system 
        select case(buff(41:41))
        case(' ','G')
            sys=SYS_GPS; tsys=TSYS_GPS
        case('R')
            sys=SYS_GLO; tsys=TSYS_UTC
        case('E')
            sys=SYS_GAL;  tsys=TSYS_GAL  ! v.2.12 
        case('S')
            sys=SYS_SBS;  tsys=TSYS_GPS
        case('J')
            sys=SYS_QZS;  tsys=TSYS_QZS  ! v.3.02 
        case('C')
            sys=SYS_CMP;  tsys=TSYS_CMP  ! v.2.12 
        case('I')
            sys=SYS_IRN;  tsys=TSYS_IRN  ! v.3.03 
        case('M')
            sys=SYS_NONE; tsys=TSYS_GPS  ! mixed 
        case default
        end select
        cycle
    elseif(index(label,'PGM / RUN BY / DATE')/=0)then
        cycle
    elseif(index(label,'COMMENT')/=0)then
        cycle
    endif
    ! file type 
    select case(mytype)
    case('O')
        call decode_obsh(fp,buff,ver,tsys,tobs,obs,nav,sta)
    case('N')
        call decode_navh(buff,nav)
    case('G')
        call decode_gnavh(buff,nav)
    end select
    if (index(label,'END OF HEADER')/=0)then
        stat=1; return
    endif
    i=i+1
    if (i>=MAXPOSHEAD .and. mytype==' ') exit  ! no rinex file 
enddo
stat=0
end subroutine

! decode obs epoch ----------------------------------------------------------
subroutine decode_obsepoch(fp, buff, ver, time, flag, sats, stat)
implicit none
integer*4, intent(in) :: fp
character(*), intent(inout) :: buff
real*8, intent(in) :: ver
type(gtime_t), intent(out) :: time
integer*4, intent(out) :: flag, sats(:), stat
integer*4 i,j,n,info
character(8) :: satid
satid=''
if(ver<=2.99)then  ! ver.2 
    read(buff(30:32),*,iostat=info) n
    if(n<=0.or.info/=0)then
        stat=0; return
    endif
    ! epoch flag: 3:new site,4:header info,5:external event 
    read(buff(29:29),*,iostat=info) flag
    if(info/=0)then ! format check
        stat=0; return
    endif
    if (3<=flag .and. flag<=5)then
        stat=n; return
    endif
    
    call str2time(buff,1,26,time,info)
    if (info/=0)then
        stat=0; return
    endif
    
    i=1; j=33
    do while(i<=n)
        if (j>=68)then
            read(fp,"(A)",iostat=info) buff
            if(info/=0) exit
            j=33
        endif
        if(i<=MAXOBS)then
            satid=adjustl(buff(j:j+2))
            sats(i)=satid2no(satid)
        endif
        i=i+1; j=j+3
    enddo
else  ! ver.3 
    read(buff(33:35),*,iostat=info) n
    if (n<=0.or.info/=0)then
        stat=0; return
    endif
    read(buff(32:32),*,iostat=info) flag
    if(info/=0)then ! format check
        stat=0; return
    endif
    if (3<=flag .and. flag<=5)then
        stat=n; return
    endif
    call str2time(buff,2,28,time,info)
    if (buff(1:1)/='>' .or. info/=0)then
        stat=0; return
    endif
endif
stat=n
end subroutine

! decode obs data -----------------------------------------------------------
subroutine decode_obsdata(fp, buff, ver, mask, index_in, obs, stat)
implicit none
integer*4, intent(in) :: fp, mask
character(*), intent(inout) :: buff
real*8, intent(in) :: ver
type(sigind_t), intent(in) :: index_in(:)
type(obsd_t), intent(out) :: obs
integer*4, intent(out) :: stat
type(sigind_t) ind
real*8 :: val(MAXOBSTYPE)
integer*4 :: lli(MAXOBSTYPE), qual(MAXOBSTYPE)
character(8) :: satid
integer*4 :: i,j,n,m,stat1,p(MAXOBSTYPE),k(16),l(16),prn,sys,info
satid=''
val=0.d0; lli=0
qual=0; stat1=1
if (ver>2.99)then
    satid=buff(1:3)
    obs%sat=satid2no(satid)
endif
call satsys(obs%sat,prn,sys)
if (obs%sat==0)then
    stat1=0
elseif(and(sys,mask)==0)then
    stat1=0
endif
! read obs data fields 
select case(sys)
case(SYS_GLO)
    ind=index_in(2)
case(SYS_GAL)
    ind=index_in(3)
case(SYS_QZS)
    ind=index_in(4)
case(SYS_SBS)
    ind=index_in(5)
case(SYS_CMP)
    ind=index_in(6)
case default
    ind=index_in(1)
end select
i=1; j=icond(ver<=2.99,1,4)
do while(i<=ind%n)
    if(ver<=2.99.and.j>80)then  ! ver.2 
        read(fp,"(A)",iostat=info) buff
        if (info/=0) exit
        j=1
    endif
    if(stat1/=0)then
        !write(*,"(A,A)")buff(j:j+15),"|"
        lli(i)=0; qual(i)=0;
        read(buff(j:j+13),"(D14.3)",iostat=info) val(i); val(i)=val(i)+ind%shift(i)
        if(ichar(buff(j+14:j+14))>=48.and.ichar(buff(j+14:j+14))<=57)then  ! add judgment of char
          read(buff(j+14:j+14),*,iostat=info) lli(i); lli(i)=and(lli(i),3)
        endif
        if(ichar(buff(j+15:j+15))>=48.and.ichar(buff(j+15:j+15))<=57)then
          read(buff(j+15:j+15),*,iostat=info) qual(i)
        endif
    endif
    i=i+1; j=j+16
enddo
if (stat1==0)then
    stat=0; return
endif
obs%P=0.d0; obs%L=0.d0; obs%D=0.d0
obs%SNR=0; obs%LLI=0; obs%qualL=0; obs%qualP=0; obs%code=0
! assign position in obs data 
i=1; n=1; m=1
do while(i<=ind%n)
    p(i)=icond(ver<=2.12,ind%frq(i)-1,ind%pos(i))
    if(ind%mytype(i)==0.and.p(i)==0)then  ! C1? index
        k(n)=i; n=n+1
    endif
    if(ind%mytype(i)==0.and.p(i)==1)then  ! C2? index
        l(m)=i; m=m+1
    endif
    i=i+1
enddo
if(ver<=2.12)then
    ! if multiple codes (C1/P1,C2/P2), select higher priority 
    if(n>2)then
        if(dabs(val(k(1)))<=1d-20 .and. dabs(val(k(2)))<=1d-20)then
            p(k(1))=-1; p(k(2))=-1
        elseif(dabs(val(k(1)))>1d-20 .and. dabs(val(k(2)))<=1d-20)then
            p(k(1))=0; p(k(2))=-1
        elseif(dabs(val(k(1)))<=1d-20 .and. dabs(val(k(2)))>1d-20)then
            p(k(1))=-1; p(k(2))=0
        elseif(ind%pri(k(2))>ind%pri(k(1)))then
            p(k(2))=0; p(k(1))=icond(NEXOBS<1,-1,NFREQ)
        else
            p(k(1))=0; p(k(2))=icond(NEXOBS<1,-1,NFREQ)
        endif
    endif
    if(m>2)then
        if(dabs(val(l(1)))<=1d-20 .and. dabs(val(l(2)))<=1d-20)then
            p(l(1))=-1; p(l(2))=-1
        elseif(dabs(val(l(1)))>1d-20 .and. dabs(val(l(2)))<=1d-20)then
            p(l(1))= 1; p(l(2))=-1
        elseif(dabs(val(l(1)))<=1d-20 .and. dabs(val(l(2)))>1d-20)then
            p(l(1))=-1; p(l(2))= 1
        elseif(ind%pri(l(2))>ind%pri(l(1)))then
            p(l(2))=1; p(l(1))=icond(NEXOBS<2,-1,NFREQ+1)
        else
            p(l(1))=1; p(l(2))=icond(NEXOBS<2,-1,NFREQ+1)
        endif
    endif
endif

! save obs data 
do i=1,ind%n
    if (p(i)<0 .or. p(i)+1>NFREQ+NEXOBS .or. dabs(val(i))<=1d-20) cycle
    select case(ind%mytype(i))
    case(0)
        obs%P(p(i)+1)=val(i)
        obs%code(p(i)+1)=ind%code(i)
        obs%qualP(p(i)+1)=icond(qual(i)>0,qual(i),1)
    case(1)
        obs%L(p(i)+1)=val(i)
        obs%LLI(p(i)+1)=lli(i)
        obs%qualL(p(i)+1)=icond(qual(i)>0,qual(i),1)
    case(2)
        obs%D(p(i)+1)=real(val(i))
    case(3)
        obs%SNR(p(i)+1)=int(val(i)*4.d0+0.5d0)
    end select
enddo
stat=1
end subroutine

! save slips ----------------------------------------------------------------
subroutine saveslips(slips, mydata)
implicit none
integer*4, intent(out) :: slips(:,:)  !(:,NFREQ)
type(obsd_t), intent(in):: mydata
integer*4 i
do i=1, NFREQ
    if(and(mydata%LLI(i),1)/=0)then
        slips(mydata%sat,i)=or(slips(mydata%sat,i),LLI_SLIP)
    endif
enddo
end subroutine

! restore slips -------------------------------------------------------------
subroutine restslips(slips, mydata)
implicit none
integer*4, intent(out) :: slips(:, :)  !(:,NFREQ)
type(obsd_t), intent(out) :: mydata
integer*4 i
do i=1, NFREQ
	if (and(slips(mydata%sat,i),1)/=0) mydata%LLI(i)=or(mydata%LLI(i),LLI_SLIP)
	slips(mydata%sat,i)=0
enddo
end subroutine

! add obs data --------------------------------------------------------------
subroutine addobsdata(obs, mydata1, stat)
implicit none
type(obs_t), intent(out) :: obs
type(obsd_t), intent(in) :: mydata1
integer*4, intent(out) :: stat

type(obsd_t), pointer :: obs_data(:)
if (obs%nmax<=obs%n)then
    if (obs%nmax<=0)then
        obs%nmax=NINCOBS
    else
        !obs%nmax=obs%nmax*2
        obs%nmax=obs%nmax+NINCOBS
    endif
    allocate(obs_data(obs%nmax))
    obs_data(1:obs%n)=obs%mydata(1:obs%n)
    if(obs%n/=0) deallocate(obs%mydata)
    obs%mydata=>obs_data
    nullify(obs_data)
endif

obs%n=obs%n+1
obs%mydata(obs%n)=mydata1
stat=1
end subroutine

! set system mask -----------------------------------------------------------
integer*4 function set_sysmask(opt)  ! -SYS=GR
implicit none
character(*), intent(in) :: opt
character(1) p
integer*4 :: mask, i
mask=SYS_NONE
i=index(opt,'-SYS=')
if (i==0)then
    set_sysmask=SYS_ALL; return
endif
i=i+5
do while(i<len_trim(opt))
    p=opt(i:i)
    select case(p)
    case('G')
        mask=or(mask,SYS_GPS)
    case('R')
        mask=or(mask,SYS_GLO)
    case('E')
        mask=or(mask,SYS_GAL)
    case('J')
        mask=or(mask,SYS_QZS)
    case('C')
        mask=or(mask,SYS_CMP)
    case('I')
        mask=or(mask,SYS_IRN)
    case('S')
        mask=or(mask,SYS_SBS)
    end select
    i=i+1
    if(p==' ') exit
enddo
set_sysmask=mask
end function

! set signal index ----------------------------------------------------------
subroutine set_index(ver, sys, opt, tobs, ind)
implicit none
real*8, intent(in) :: ver
integer*4, intent(in) :: sys
character(*), intent(in) :: opt
character(3), intent(in) :: tobs(MAXOBSTYPE)
type(sigind_t), intent(out) :: ind
character(1) p
real*8 shift
integer*4 :: i,j,k,n,info
logical*1 ltmp

ltmp=.false.
i=1; n=0
do while(tobs(i)/='end')
    if(tobs(i)=='')then  !? to be verified
        ind%frq(i)=0; ind%code(i)=0
        ind%mytype(i)=0;
        if(sys==SYS_GPS)then
            ind%pri(i)=6
        elseif(sys==SYS_SBS.or.sys==SYS_GLO.or.sys==SYS_GAL.or.sys==SYS_QZS.or.sys==SYS_CMP.or.sys==SYS_IRN)then
            ind%pri(i)=14
        elseif(sys==SYS_ALL.or.sys==SYS_LEO.or.sys==SYS_NONE)then
            ind%pri(i)=0
        endif
        ind%pos(i)=-1
        i=i+1; if(i>size(tobs,1)) exit
        n=n+1; cycle
    endif
    call obs2code(tobs(i)(2:3),ind%frq(i),ind%code(i))
    info=index(obscodes,tobs(i)(1:1))
    ind%mytype(i)=icond(info/=0,info-1,0)
    ind%pri(i)=getcodepri(sys,ind%code(i),opt)
    ind%pos(i)=-1
    if (sys==SYS_CMP)then
        if (ind%frq(i)==5) ind%frq(i)=2; ! B2 
        if (ind%frq(i)==4) ind%frq(i)=3; ! B3 
    endif
    i=i+1; n=n+1
    if(i>size(tobs,1)) exit
enddo
! assign index for highest priority code 
do i=1,NFREQ
    k=0
    do j=1,n
        if(k<1)then
            ltmp=.true.
        elseif(ind%pri(j)>ind%pri(k))then
            ltmp=.true.
        else
            ltmp=.false.
        endif
        if(ind%frq(j)==i.and.ind%pri(j)/=0.and.ltmp) k=j
    enddo
    if(k<1) cycle
    do j=1,n
        if(ind%code(j)==ind%code(k)) ind%pos(j)=i-1
    enddo
enddo
! assign index of extended obs data 
! do i=1,NEXOBS
!     do j=1,n
!         if(ind%code(j)/=0 .and. ind%pri(j)/=0 .and. ind%pos(j)<0) exit
!     enddo
!     if(j>n) exit
!     do k=1,n
!         if(ind%code(k)==ind%code(j)) ind%pos(k)=NFREQ+i-1
!     enddo
! enddo
do i=1,n
    if (ind%code(i)==0 .or. ind%pri(i)==0 .or. ind%pos(i)>=0) cycle
enddo
ind%n=n
end subroutine

! read rinex obs data body --------------------------------------------------
subroutine readrnxobsb(fp, opt, ver, tsys, tobs, flag, mydata, sta, stat)
implicit none
integer*4, intent(in) :: fp, tsys
character(*), intent(in) :: opt
real*8, intent(in) :: ver
character(3), intent(in) :: tobs(:,:)
integer*4, intent(out) :: flag, stat
type(obsd_t), intent(out) :: mydata(:)
type(sta_t), intent(in) :: sta
type(gtime_t) :: time
type(sigind_t) :: myindex(7)
type(obs_t) obstmp
character(MAXRNXLEN) buff
integer*4 :: i,n,nsat,sats(MAXOBS),mask,info
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
    ! decode obs epoch 
    if(i==0)then
        call decode_obsepoch(fp,buff,ver,time,flag,sats,nsat)
        if(nsat<=0) cycle
    elseif(flag<=2 .or. flag==6)then
        mydata(n)%time=time
        mydata(n)%sat=sats(i)
        ! decode obs data 
        call decode_obsdata(fp,buff,ver,mask,myindex,mydata(n),info)
        if(info/=0 .and. n<=MAXOBS) n=n+1
    elseif(flag==3 .or. flag==4)then
        ! decode obs header 
        ! call decode_obsh(fp,buff,ver,tsys,tobs,NULL,NULL,sta)
    endif
    i=i+1
    if(i>nsat)then
        stat=n-1; return
    endif
enddo
stat=-1
end subroutine

! read rinex obs ------------------------------------------------------------
subroutine readrnxobs(fp, ts, te, tint, opt, rcv, ver, tsys, tobs, obs, sta, stat)
implicit none
integer*4, intent(in) :: fp, rcv, tsys
type(gtime_t), intent(in) :: ts, te
real*8, intent(in) :: tint, ver
character(*), intent(in) :: opt
character(3), intent(in) :: tobs(:,:)
type(obs_t), intent(out) :: obs
type(sta_t), intent(in) :: sta
integer*4, intent(out) :: stat
type(obsd_t), pointer :: mydata(:)
integer*4 :: i,n,flag,stat1,slips(MAXSAT,NFREQ)
real*8 ti
flag=0; stat1=0; slips=0; ti=0
if(rcv>MAXRCV)then
    stat=0; return
endif
allocate(mydata(MAXOBS))
if(int(obs%tint)/=0.and.int(tint)/=0)then
    obs%tint=iminmul(int(obs%tint),int(tint)); ti=obs%tint  ! change time interval
endif
if(int(obs%tint)/=0.and.int(tint)==0) ti=0
if(int(obs%tint)==0.and.int(tint)/=0)then
    obs%tint=int(tint); ti=int(tint)
endif
if(int(obs%tint)==0.and.int(tint)==0) ti=0

! read rinex obs data body 
do while(.true.)
    call readrnxobsb(fp,opt,ver,tsys,tobs,flag,mydata,sta,n)
    if(.not.(n>=0 .and. stat1>=0)) exit
    do i=1,n
        ! utc % gpst 
        if (tsys==TSYS_UTC) mydata(i)%time=utc2gpst(mydata(i)%time)
        ! save cycle-slip 
        call saveslips(slips,mydata(i))
    enddo
    if(n>0 .and. screent(mydata(1)%time,ts,te,ti)==0) cycle
    do i=1,n
        ! restore cycle-slip 
        call restslips(slips,mydata(i))
        ! save obs data 
        mydata(i)%rcv=1
        call addobsdata(obs,mydata(i),stat1)
        if (stat1<0) exit
    enddo
enddo
deallocate(mydata)
stat=stat1
end subroutine

! decode ephemeris ----------------------------------------------------------
subroutine decode_eph(ver, sat, toc, mydata, eph, stat)
implicit none
real*8, intent(in) :: ver, mydata(:)
integer*4, intent(in) :: sat
type(gtime_t), intent(in) :: toc
type(eph_t), intent(out) :: eph
integer*4, intent(out) :: stat
integer*4 sys, prn

call satsys(sat,prn,sys)
if (and(sys,or(SYS_GPS,or(SYS_GAL,or(SYS_QZS,or(SYS_CMP,SYS_IRN)))))==0)then
    stat=0; return
endif
call init_eph(eph)

eph%sat=sat
eph%toc=toc
eph%f0=mydata(1)
eph%f1=mydata(2)
eph%f2=mydata(3)

eph%A=SQR(mydata(11)); eph%e=mydata( 9); eph%i0  =mydata(16); eph%OMG0=mydata(14)
eph%omg =mydata(18); eph%M0 =mydata( 7); eph%deln=mydata( 6); eph%OMGd=mydata(19)
eph%idot=mydata(20); eph%crc=mydata(17); eph%crs =mydata( 5); eph%cuc =mydata( 8)
eph%cus =mydata(10); eph%cic=mydata(13); eph%cis =mydata(15)
if(sys==SYS_GPS .or. sys==SYS_QZS)then
    eph%iode=int(mydata( 4))      ! IODE 
    eph%iodc=int(mydata(27))      ! IODC 
    eph%toes=   (mydata(12))      ! toe (s) in gps week 
    eph%week=int(mydata(22))      ! gps week 
    eph%toe=adjweek(gpst2time(eph%week,mydata(12)),toc)
    eph%ttr=adjweek(gpst2time(eph%week,mydata(28)),toc)
    
    eph%code=int(mydata(21))      ! GPS: codes on L2 ch 
    eph%svh =int(mydata(25))      ! sv health 
    eph%sva =uraindex(mydata(24)) ! ura (m%index) 
    eph%flag=int(mydata(23))      ! GPS: L2 P data flag 
    
    eph%tgd(1)=   mydata(26)      ! TGD 
    if(sys==SYS_GPS)then
        eph%fit=mydata(29)        ! fit interval (h) 
    else
        eph%fit=rcond(dabs(mydata(29))<=1d-20,1.d0,2.d0) ! fit interval (0:1h,1:>2h) 
    endif
elseif(sys==SYS_GAL)then
    eph%iode=int(mydata( 4))      ! IODnav 
	eph%toes=   (mydata(12))      ! toe (s) in galileo week 
	eph%week=int(mydata(22))      ! gal week = gps week 
	eph%toe=adjweek(gpst2time(eph%week,mydata(12)),toc)
	eph%ttr=adjweek(gpst2time(eph%week,mydata(28)),toc)
	
	eph%code=int(mydata(21))      ! data sources 
								  ! bit 0 set: I/NAV E1-B 
								  ! bit 1 set: F/NAV E5a-I 
								  ! bit 2 set: F/NAV E5b-I 
								  ! bit 8 set: af0-af2 toc are for E5a.E1 
								  ! bit 9 set: af0-af2 toc are for E5b.E1 
	eph%svh =int(mydata(25))      ! sv health 
								  ! bit     0: E1B DVS 
								  ! bit   1-2: E1B HS 
								  ! bit     3: E5a DVS 
								  ! bit   4-5: E5a HS 
								  ! bit     6: E5b DVS 
								  ! bit   7-8: E5b HS 
	eph%sva =sisa_index(mydata(24)) ! sisa (m->index) 
	eph%tgd(1)=mydata(26)         ! BGD E5a/E1 
	eph%tgd(2)=mydata(27)         ! BGD E5b/E1 
elseif(sys==SYS_CMP)then
elseif(sys==SYS_IRN)then
endif
stat=1
end subroutine

subroutine decode_geph(ver, sat, toc_in, mydata, geph, stat)
implicit none
real*8, intent(in) :: ver, mydata(:)
integer*4, intent(in) :: sat
type(gtime_t), intent(in) :: toc_in
type(geph_t), intent(out) :: geph
integer*4, intent(out) :: stat
type(gtime_t) toc, tof
real*8 tow,tod
integer*4 week,dow,prn,sys
toc=toc_in
call satsys(sat,prn,sys)
if(sys/=SYS_GLO)then
    stat=0; return
endif
call init_geph(geph)
geph%sat=sat

! toc rounded by 15 min in utc 
call time2gpst(toc,week,tow)
toc=gpst2time(week,floor((tow+450.d0)/900.d0)*900.d0)
dow=int(floor(tow/86400.d0))

! time of frame in utc 
tod=rcond(ver<=2.99,mydata(3),dmod(mydata(3),86400.d0))  ! tod (v.2), tow (v.3) in utc 
tof=gpst2time(week,tod+dow*86400.d0)
tof=adjday(tof,toc)
geph%toe=utc2gpst(toc)    ! toc (gpst) 
geph%tof=utc2gpst(tof)    ! tof (gpst) 
 
! iode = tb (7bit), tb =index of UTC+3H within current day 
geph%iode=int(dmod(tow+10800.d0,86400.d0)/900.d0+0.5d0)
geph%taun=-mydata(1)       ! -taun 
geph%gamn= mydata(2)       ! +gamman 

geph%pos(1)=mydata(4)*1d3; geph%pos(2)=mydata(8)*1d3;  geph%pos(3)=mydata(12)*1d3
geph%vel(1)=mydata(5)*1d3; geph%vel(2)=mydata(9)*1d3;  geph%vel(3)=mydata(13)*1d3
geph%acc(1)=mydata(6)*1d3; geph%acc(2)=mydata(10)*1d3; geph%acc(3)=mydata(14)*1d3

geph%svh=int(mydata( 7))
geph%frq=int(mydata(11))
geph%age=int(mydata(15))

! some receiver output >128 for minus frequency number 
if (geph%frq>128) geph%frq=geph%frq-256
!if (geph%frq<MINFREQ_GLO .or. MAXFREQ_GLO<geph%frq) {}
stat=1
end subroutine

! read rinex navigation data body -------------------------------------------
subroutine readrnxnavb(fp, opt, ver, sys_in, mytype, eph, geph, stat)
implicit none
integer*4, intent(in) :: fp, sys_in
character(*), intent(in) :: opt
real*8, intent(in) :: ver
integer*4, intent(out) :: mytype, stat
type(eph_t), intent(out) :: eph
type(geph_t), intent(out) :: geph
type(gtime_t) toc
real*8 mydata(64)
integer*4 :: sys,i,j,prn,sat,sp,mask,info,prntmp
character(MAXRNXLEN) buff, p
character(8) :: id

i=1; sat=0; sp=3
id=''
! set system mask 
sys=sys_in
mask=set_sysmask(opt)
do while(.true.)
    read(fp,"(A)",iostat=info) buff
    if(info/=0) exit
    if(i==1)then
        ! decode satellite field 
        if(ver>=3.d0 .or. sys==SYS_GAL .or. sys==SYS_QZS)then ! ver.3 or GAL/QZS 
            id=buff(1:3)
            sat=satid2no(id); sp=4
            if(ver>=3.d0) call satsys(sat,prntmp,sys)
        else
            prn=int(str2num(buff,1,2))
            if (sys==SYS_SBS) sat=satno(SYS_SBS,prn+100)
            if(sys==SYS_GLO)then
                sat=satno(SYS_GLO,prn)
            elseif(93<=prn .and. prn<=97)then
                sat=satno(SYS_QZS,prn+100)
            else
                sat=satno(SYS_GPS,prn)
            endif
        endif
        ! decode toc field 
        call str2time(buff(1+sp:),1,19,toc,info)
        if (info/=0)then
            stat=0; return
        endif
        ! decode data fields 
        p=buff(1+sp+19:)
        do j=1,3
            mydata(i)=str2num(p,1,19)
            i=i+1; p=p(20:)
        enddo
    else
        ! decode data fields 
        p=buff(1+sp:)
        do j=1,4
            mydata(i)=str2num(p,1,19)
            i=i+1; p=p(20:)
        enddo
        ! decode ephemeris 
        if(sys==SYS_GLO .and. i>15)then
            if (and(mask,sys)==0)then
                stat=0; return
            endif
            mytype=1
            call decode_geph(ver,sat,toc,mydata,geph,stat); return
        elseif(i>31)then
            if (and(mask,sys)==0)then
                stat=0; return
            endif
            mytype=0
            call decode_eph(ver,sat,toc,mydata,eph,stat); return
        endif
    endif
enddo
stat=-1
end subroutine

! add ephemeris to navigation data ------------------------------------------
subroutine add_eph(nav, eph, stat)
implicit none
type(nav_t), intent(out) :: nav
type(eph_t), intent(in) :: eph
integer*4, intent(out) :: stat
type(eph_t), pointer :: nav_eph(:)

if(nav%nmax<=nav%n)then
    nav%nmax=nav%nmax+1024
    allocate(nav_eph(nav%nmax))
    nav_eph(1:nav%n)=nav%eph(1:nav%n)
    if(nav%n/=0) deallocate(nav%eph)
    nav%eph=>nav_eph
    nullify(nav_eph)
endif
nav%n=nav%n+1
nav%eph(nav%n)=eph
stat=1
end subroutine
    
subroutine add_geph(nav, geph, stat)
implicit none
type(nav_t), intent(out) :: nav
type(geph_t), intent(in) :: geph
integer*4, intent(out) :: stat
type(geph_t), pointer :: nav_geph(:)

if(nav%ngmax<=nav%ng)then
    nav%ngmax=nav%ngmax+1024
    allocate(nav_geph(nav%ngmax))
    nav_geph(1:nav%ng)=nav%geph(1:nav%ng)
    if(nav%ng/=0) deallocate(nav%geph)
    nav%geph=>nav_geph
    nullify(nav_geph)
endif
nav%ng=nav%ng+1
nav%geph(nav%ng)=geph
stat=1
end subroutine

! read rinex nav/gnav/geo nav -----------------------------------------------
subroutine readrnxnav(fp, opt, ver, sys, nav, stat)
implicit none
integer*4, intent(in) :: fp, sys
character(*), intent(in) :: opt
real*8, intent(in) :: ver
type(nav_t), intent(out) :: nav
integer*4, intent(out) :: stat
type(eph_t) eph
type(geph_t) geph
integer*4 stat1, mytype
    
! read rinex navigation data body 
do while(.true.)
    call readrnxnavb(fp,opt,ver,sys,mytype,eph,geph,stat1)
    if(stat1<0) exit
    ! add ephemeris to navigation data 
    if(stat1/=0)then
        select case(mytype)
        case(1)
            call add_geph(nav,geph,stat1)
        case default
            call add_eph (nav,eph, stat1)
        end select
        if(stat1==0)then
            stat=0; return
        endif
    endif
enddo
if(nav%n>0.or.nav%ng>0)then
    stat=1
else
    stat=0
endif
end subroutine

! read rinex file -----------------------------------------------------------
subroutine readrnxfp(fp, ts, te, tint, opt, flag, index, mytype, obs, nav, sta, stat)
implicit none
integer*4, intent(in) :: fp, flag, index
type(gtime_t), intent(in) :: ts, te
real*8, intent(in) :: tint
character(*), intent(in) :: opt
character(*), intent(out) :: mytype
type(obs_t), intent(out) :: obs
type(nav_t), intent(out) :: nav
type(sta_t), intent(out) :: sta
integer*4, intent(out) :: stat
real*8 ver
integer*4 :: sys,tsys,info
character(3) :: tobs(NUMSYS,MAXOBSTYPE)

tsys=TSYS_GPS
tobs=''
! read rinex header 
call readrnxh(fp,ver,mytype,sys,tsys,tobs,obs,nav,sta,info)
if (info==0)then
    stat=0; return
endif
! read rinex body 
select case(mytype(1:1))
case('O')
    call readrnxobs(fp,ts,te,tint,opt,index,ver,tsys,tobs,obs,sta,stat); return
case('N')
    call readrnxnav(fp,opt,ver,sys,nav,stat); return
case('G')
    call readrnxnav(fp,opt,ver,SYS_GLO,nav,stat); return
case('L')
    call readrnxnav(fp,opt,ver,SYS_GAL,nav,stat); return  ! extension
end select
stat=0
end subroutine

! uncompress and read rinex file --------------------------------------------
subroutine readrnxfile(filepath, ts, te, tint, opt, flag, index, mytype, obs, nav, sta, stat)
implicit none
character(*), intent(in) :: filepath, opt
type(gtime_t), intent(in) :: ts, te
real*8, intent(in) :: tint
integer*4, intent(in) :: flag, index
character(*), intent(out) :: mytype
type(obs_t), intent(out) :: obs
type(nav_t), intent(out) :: nav
type(sta_t), intent(out) :: sta
integer*4, intent(out) :: stat
integer*4 :: fp=FPREAD, stat1, info

call init_sta(sta)
open(unit=fp,file=filepath,status='old',iostat=info)
if(info/=0)then
    stat=-1; return
endif
! read rinex file 
call readrnxfp(fp,ts,te,tint,opt,flag,index,mytype,obs,nav,sta,stat1)
if(sta%name=='')then
    call decodemarker(filepath, sta%name)
endif
close(unit=fp,status='keep')
stat=stat1
end subroutine

! read rinex obs and nav files ------------------------------------------------
! read rinex obs and nav files
! args   : char *file    I      file (wild-card * expanded) ('': stdin)
!          integer*4   rcv     I      receiver number for obs data
!         (gtime_t ts)   I      observation time start (ts.time==0: no limit)
!         (gtime_t te)   I      observation time end   (te.time==0: no limit)
!         (real*8 tint)  I      observation time interval (s) (0:all)
!          char  *opt    I      rinex options (see below,'': no option)
!          obs_t *obs    IO     observation data   (NULL: no input)
!          nav_t *nav    IO     navigation data    (NULL: no input)
!          sta_t *sta    IO     station parameters (NULL: no input)
! return : status (1:ok,0:no data,-1:error)
! notes  : read data are appended to obs and nav struct
!          before calling the function, obs and nav should be initialized.
!          observation data and navigation data are not sorted.
!          navigation data may be duplicated.
!          call sortobs() or uniqnav() to sort data or delete duplicated eph.
!
!          read rinex options (separated by spaces) :
!
!            -GLss(=shift): select GPS signal ss (ss: RINEX 3 code, '1C','2W'...)
!            -RLss(=shift): select GLO signal ss
!            -ELss(=shift): select GAL signal ss
!            -JLss(=shift): select QZS signal ss
!            -CLss(=shift): select BDS signal ss
!            -ILss(=shift): select IRN signal ss
!            -SLss(=shift): select SBS signal ss
!
!                 shift: carrier phase shift to be added (cycle)
!            
!            -SYS=sys(,sys...): select navi systems
!                               (sys=G:GPS,R:GLO,E:GAL,J:QZS,C:BDS,I:IRN,S:SBS)
!
!-----------------------------------------------------------------------------
subroutine readrnxt(filepath, rcv, ts, te, tint, opt, obs, nav, sta, stat)
implicit none
character(*), intent(in) :: filepath, opt
integer*4, intent(in) :: rcv
type(gtime_t), intent(in) :: ts, te
real*8, intent(in) :: tint
type(obs_t), intent(out) :: obs
type(nav_t), intent(out) :: nav
type(sta_t), intent(out) :: sta
integer*4, intent(out) :: stat
integer*4 :: i,n,stat1
character(1) :: mytype
character(150) :: files(2)
stat1=0
mytype=''; files=''
files(1)=filepath; n=1
! read rinex files 
do i=1,n
    if(stat1<0) exit
    call readrnxfile(files(i),ts,te,tint,opt,0,rcv,mytype,obs,nav,sta,stat1)
enddo
stat=stat1
end subroutine

subroutine readrnx(filepath, rcv, opt, obs, nav, sta, stat)
implicit none
character(*), intent(in) :: filepath, opt
integer*4, intent(in) :: rcv
type(obs_t), intent(out) :: obs
type(nav_t), intent(out) :: nav
type(sta_t), intent(out) :: sta
integer*4, intent(out) :: stat
type(gtime_t) :: t
t=gtime_t(0,0.d0)
call readrnxt(filepath,rcv,t,t,0.d0,opt,obs,nav,sta,stat)
end subroutine

end module rinex_f90_
