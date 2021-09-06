!! reference and acknowledgement to RTKLIB
!! modified by Kunlun Zhang
! constants -----------------------------------------------------------------

module rtkcmn_f90_
use rtklib_h_
implicit none

!integer*4, parameter :: POLYCRC32  = 3988292384 !#EDB88320   ! CRC32 polynomial 
!integer*4, parameter :: POLYCRC24Q = 25578747   !#1864CFB    ! CRC24Q polynomial 

private gpst0, gst0, bdt0, leaps
private obscodes, obsfreqs, codepris, timeoffset_

real*8, parameter :: MAX_VAR_EPH = 9d4            ! max variance eph to reject satellite (m^2) 
real*8, parameter :: gpst0(6)=(/1980,1, 6,0,0,0/) ! gps time reference 
real*8, parameter :: gst0 (6)=(/1999,8,22,0,0,0/) ! galileo system time reference 
real*8, parameter :: bdt0 (6)=(/2006,1, 1,0,0,0/) ! beidou time reference 

real*8, parameter :: lam_carr(MAXFREQ)=(/&        ! carrier wave length (m)
    CLIGHT/FREQ1,CLIGHT/FREQ2,CLIGHT/FREQ5,CLIGHT/FREQ6,CLIGHT/FREQ7,&
    CLIGHT/FREQ8,CLIGHT/FREQ9/)

! leap seconds (y,m,d,h,m,s,utc-gpst) 
real*8, parameter :: leaps(19,7)=reshape((/&
    2017,1,1,0,0,0,-18,&
    2015,7,1,0,0,0,-17,&
    2012,7,1,0,0,0,-16,&
    2009,1,1,0,0,0,-15,&
    2006,1,1,0,0,0,-14,&
    1999,1,1,0,0,0,-13,&
    1997,7,1,0,0,0,-12,&
    1996,1,1,0,0,0,-11,&
    1994,7,1,0,0,0,-10,&
    1993,7,1,0,0,0, -9,&
    1992,7,1,0,0,0, -8,&
    1991,1,1,0,0,0, -7,&
    1990,1,1,0,0,0, -6,&
    1988,1,1,0,0,0, -5,&
    1985,7,1,0,0,0, -4,&
    1983,7,1,0,0,0, -3,&
    1982,7,1,0,0,0, -2,&
    1981,7,1,0,0,0, -1,&
    0   ,0,0,0,0,0,  0 &
/),shape(leaps),order=(/2,1/))

! observation code strings 
character(2), save :: obscodes(60)=(/&
    '  ','1C','1P','1W','1Y', '1M','1N','1S','1L','1E',&  !  0- 9 
    '1A','1B','1X','1Z','2C', '2D','2S','2L','2X','2P',&  ! 10-19 
    '2W','2Y','2M','2N','5I', '5Q','5X','7I','7Q','7X',&  ! 20-29 
    '6A','6B','6C','6X','6Z', '6S','6L','8L','8Q','8X',&  ! 30-39 
    '2I','2Q','6I','6Q','3I', '3Q','3X','1I','1Q','5A',&  ! 40-49 
    '5B','5C','9A','9B','9C', '9X','  ','  ','  ','  ' &  ! 50-59
/)

! 1:L1/E1, 2:L2/B1, 3:L5/E5a/L3, 4:L6/LEX/B3, 5:E5b/B2, 6:E5(a+b), 7:S
integer(4), save :: obsfreqs(60)=(/&
    0, 1, 1, 1, 1,  1, 1, 1, 1, 1,&  !  0- 9
    1, 1, 1, 1, 2,  2, 2, 2, 2, 2,&  ! 10-19
    2, 2, 2, 2, 3,  3, 3, 5, 5, 5,&  ! 20-29
    4, 4, 4, 4, 4,  4, 4, 6, 6, 6,&  ! 30-39
    2, 2, 4, 4, 3,  3, 3, 1, 1, 3,&  ! 40-49
    3, 3, 7, 7, 7,  7, 0, 0, 0, 0 &  ! 50-59
/)

! code priority table 
! L1/E1      L2/B1        L5/E5a/L3 L6/LEX/B3 E5b/B2    E5(a+b)  S 
character(16), save :: codepris(7,MAXFREQ)=reshape((/&
    'CPYWMNSL  ','PYWCMNDSLX','IQX       ','          ','          ','          ','          ',&  ! GPS
    'PC        ','PC        ','IQX       ','          ','          ','          ','          ',&  ! GLO
    'CABXZ     ','          ','IQX       ','ABCXZ     ','IQX       ','IQX       ','          ',&  ! GAL
    'CSLXZ     ','SLX       ','IQX       ','SLX       ','          ','          ','          ',&  ! QZS
    'C         ','          ','IQX       ','          ','          ','          ','          ',&  ! SBS
    'IQX       ','IQX       ','IQX       ','IQX       ','IQX       ','          ','          ',&  ! BDS
    '          ','          ','ABCX      ','          ','          ','          ','ABCX      ' &  ! IRN
/),shape(codepris),order=(/2,1/))

real*8, save :: timeoffset_=0.d0    ! time offset (s)

contains
real*8 function SQR(x)
implicit none
real*8, intent(in) :: x
SQR=((x)*(x))
end function

real*8 function SQRT2(x)
implicit none
real*8, intent(in) :: x
SQRT2=rcond(x<=0.d0.or.x/=x,0.d0,dsqrt(x))
end function

! get fname from fpath ---------------------------------------------------
subroutine getfname(fpath, fname)
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

! get fdir from fpath ---------------------------------------------------
subroutine getfdir(fpath, fdir)
implicit none
character(*), intent(in) :: fpath
character(*), intent(out) :: fdir
integer*4 ind1, ind2
fdir=fpath
ind1=index(fpath,'/',BACK = .TRUE.)
ind2=index(fpath,'\\',BACK = .TRUE.)
if(ind1/=0)then
    fdir=trim(fpath(1:ind1)); return
endif
if(ind2/=0)then
    fdir=trim(fpath(1:ind2+1)); return
endif
end subroutine

! confirm file status
! status ==> 0-none, 1-exist
integer*4 function flexist(flname)
implicit none
character(*), intent(in) :: flname
logical*2 :: alive
inquire(file=flname,exist=alive)
if(alive)then
    flexist=1; return
else
    flexist=0; return
endif
end function

! confirm the type of rnxfile
! rnxtype ==> 0-none, 1-stdnav(long),  2-stdobs(long)
!                     3-stdnav(short), 4-stdobs(short)
integer*4 function getrnxtyp(flname)
implicit none
character(*), intent(in) :: flname
integer*4 :: len1
len1=len_trim(adjustl(flname))
getrnxtyp=0
if(len1==38)then
    if(index(flname,"O.rnx")/=0) getrnxtyp=2
elseif(len1==34)then
    if(index(flname,"N.rnx")/=0) getrnxtyp=1
elseif(len1==12)then
    if(flname(8:9)=="0." .and. flname(12:12)=="o") getrnxtyp=4
    if(flname(8:9)=="0." .and. flname(12:12)=="n") getrnxtyp=3
    if(flname(8:9)=="0." .and. flname(12:12)=="p") getrnxtyp=3
endif
end function

! get mjd of rnxname
! mjd ==> 0-error, 1-normal
integer*4 function getrnxmjd(flname)
implicit none
character(*), intent(in) :: flname
integer*4 ityp, info, iyear, idoy
ityp=getrnxtyp(flname)
getrnxmjd=0  ! if(ityp==0)
if(ityp==1 .or. ityp==2)then
    read(flname(13:16),*,iostat=info) iyear
    read(flname(17:19),*,iostat=info) idoy
    getrnxmjd=ydoy2mjd_s(iyear, idoy)
elseif(ityp==3 .or. ityp==4)then
    read(flname(10:11),*,iostat=info) iyear
    read(flname(5:7)  ,*,iostat=info) idoy
    if(0<=iyear.and.iyear<=70) iyear=iyear+2000  ! 1971~2070
    if(70<iyear.and.iyear<=99) iyear=iyear+1900
    getrnxmjd=ydoy2mjd_s(iyear, idoy)
endif
end function

! get the next rnxname
subroutine getrnx_nname(flname, flname2, offset1)
implicit none
character(*), intent(in) :: flname
character(*), intent(out) :: flname2
integer*4, intent(in) :: offset1
integer*4 mjd1, mjd2, ityp
integer*4 iyear2, idoy2

flname2=""
mjd1=getrnxmjd(flname)
if(mjd1==0) return

mjd2=mjd1+offset1
if(mjd2<=0) mjd2=0
call mjd2ydoy_s(mjd2, iyear2, idoy2)

flname2=flname
ityp=getrnxtyp(flname2)
if(ityp==1 .or. ityp==2)then
    write(flname2(13:19),"(I4.4I3.3)") iyear2, idoy2
elseif(ityp==3 .or. ityp==4)then
    write(flname2(10:11),"(I2.2)") mod(iyear2,100)
    write(flname2(5:7)  ,"(I3.3)") idoy2
endif
end subroutine

! conditional operation (integer) --------------------------------------------
integer*4 function icond(key,a,b)
implicit none
logical*4, intent(in) :: key
integer*4, intent(in) :: a,b
if(key)then
    icond=a; return
else
    icond=b; return
endif
end function

! conditional operation (real) -----------------------------------------------
real*8 function rcond(key,a,b)
implicit none
logical*4, intent(in) :: key
real*8, intent(in) :: a,b
if(key)then
    rcond=a; return
else
    rcond=b; return
endif
end function

! average of array (except 0) ------------------------------------------------
real*8 function raver(x)
implicit none
real*8, intent(in) :: x(:)
integer*4 :: i, numZero
numZero=0
do i=1,size(x,dim=1)
    if(dabs(x(i))<=1d-20) numZero=numZero+1
enddo
raver=sum(x)/(size(x,dim=1)-numZero)
end function

! rms of array (except 0) ----------------------------------------------------
real*8 function rrms(x)
implicit none
real*8, intent(in) :: x(:)
real*8, allocatable :: x1(:)
real*8 aver
integer*4 :: i, numZero
numZero=0
allocate(x1(size(x,dim=1)))
x1=x
do i=1,size(x1,dim=1)
    if(dabs(x1(i))<=1d-20) numZero=numZero+1
enddo
aver=raver(x1)
do i=1,size(x1,dim=1)
    if(dabs(x1(i))<=1d-20) x1(i)=aver
enddo
x1=(x1-aver)*(x1-aver)
rrms=dsqrt(sum(x1)/(size(x1,dim=1)-numZero))
deallocate(x1)
end function

! delete values that exceed 3σ -----------------------------------------------
subroutine errorfilter(a, b)  ! (except 0)
implicit none
real*8, intent(in) :: a(:)
real*8, intent(out) :: b(:)
integer*4 i
real*8 aver, rms
aver=raver(a)
rms=rrms(a)
b=a
do i=1,size(a,dim=1)
    if(dabs(b(i)-aver)>=rms*3d0) b(i)=0
enddo
end subroutine

integer*4 function ROUND(x)
implicit none
real*8, intent(in) :: x
ROUND=int(floor(x+0.5d0))
end function

real*8 function ROUNDF(x, n)
implicit none
real*8, intent(in) :: x
integer*4, intent(in) :: n
integer*4 n1, tmp
n1=icond(n<=0,0,n)
ROUNDF=ROUND(x*10**n)*1.d0/(10**n)
end function

real*8 function sqvar(covar)
implicit none
real*8, intent(in) :: covar
sqvar=rcond(covar<0.d0,-dsqrt(-covar),dsqrt(covar))
end function

subroutine InitGlobal()
implicit none
chisqr=(/&  ! chi-sqr(n) (alpha=0.001)
    10.8d0, 13.8d0, 16.3d0, 18.5d0, 20.5d0, 22.5d0, 24.3d0, 26.1d0, 27.9d0, 29.6d0,&
    31.3d0, 32.9d0, 34.5d0, 36.1d0, 37.7d0, 39.3d0, 40.8d0, 42.3d0, 43.8d0, 45.3d0,&
    46.8d0, 48.3d0, 49.7d0, 51.2d0, 52.6d0, 54.1d0, 55.5d0, 56.9d0, 58.3d0, 59.7d0,&
    61.1d0, 62.5d0, 63.9d0, 65.2d0, 66.6d0, 68.0d0, 69.3d0, 70.7d0, 72.1d0, 73.4d0,&
    74.7d0, 76.0d0, 77.3d0, 78.6d0, 80.0d0, 81.3d0, 82.6d0, 84.0d0, 85.4d0, 86.7d0,&
    88.0d0, 89.3d0, 90.6d0, 91.9d0, 93.3d0, 94.7d0, 96.0d0, 97.4d0, 98.7d0, 100.d0,&
    101.d0, 102.d0, 103.d0, 104.d0, 105.d0, 107.d0, 108.d0, 109.d0, 110.d0, 112.d0,&
    113.d0, 114.d0, 115.d0, 116.d0, 118.d0, 119.d0, 120.d0, 122.d0, 123.d0, 125.d0,&
    126.d0, 127.d0, 128.d0, 129.d0, 131.d0, 132.d0, 133.d0, 134.d0, 135.d0, 137.d0,&
    138.d0, 139.d0, 140.d0, 142.d0, 143.d0, 144.d0, 145.d0, 147.d0, 148.d0, 149.d0/)

ura_value=(/&
    2.4d0, 3.4d0, 4.85d0, 6.85d0, 9.65d0, 13.65d0, 24.d0, 48.d0, 96.d0, &
    192.d0, 384.d0, 768.d0, 1536.d0, 3072.d0, 6144.d0/)

! defaults processing options
prcopt_default = prcopt_t(&
    PMODE_SINGLE, 0, 2, or(SYS_GPS,or(SYS_GLO,SYS_GAL)),&  ! mode, soltype:forward, nf:L1+L2, navsys
    10.d0*D2R, 0, 3, 1, 1,&                                  ! elmin, sateph:brdc, ionoopt:dual, tropopt:saas, niter
    30.d0, 30.d0, 0, '', (/0,0,0,0,1,0/))                    ! maxinno, maxgdop, exsats, rnxopt, posopt
! prcopt_default = prcopt_t(&
!     PMODE_SINGLE, 0, 2, SYS_GPS,&                          ! mode, soltype:forward, nf:L1+L2, navsys
!     10.d0*D2R, 0, 3, 1, 1,&                                ! elmin, sateph:brdc, ionoopt:dual, tropopt:saas, niter
!     30.d0, 30.d0, 0, '', (/0,0,0,0,1,0/))                  ! maxinno, maxgdop, exsats, rnxopt, posopt

! defaults solution output options
solopt_default = solopt_t(&
    SOLF_XYZ,TIMES_GPST,0,3,&    ! posf,times,timef,timeu */
    0,1,1,0,0,0,0,&              ! degf,outhead,outopt,outvel,datum,height,geoid */
    0,0,0,0,&                    ! solstatic,sstat,trace,issingle */
    0.d0,&                       ! nmeaintv */
    " ","",&                     ! separator/program name */
    0)                           ! maxsolstd */

end subroutine

subroutine init_anyvar(var)      ! 1-int; 2-double; 3-char*; 0-null
implicit none
type(anyvar), intent(out) :: var
var%mytype=0
var%ivalue=0
var%rvalue=0.d0
var%cvalue=''
end subroutine

type(anyvar) function ivar(value)  ! 1-int
implicit none
integer*4, intent(in) :: value
call init_anyvar(ivar)
ivar%mytype=1
ivar%ivalue=value
end function

type(anyvar) function rvar(value)  ! 2-double
implicit none
real*8, intent(in) :: value
call init_anyvar(rvar)
rvar%mytype=2
rvar%rvalue=value
end function

type(anyvar) function cvar(value)  ! 3-char*
implicit none
character(*), intent(in) :: value
call init_anyvar(cvar)
cvar%mytype=3
cvar%cvalue=value
end function

type(anyvar) function vvar()
implicit none
call init_anyvar(vvar)
end function

! satellite number to satellite system ---------------------------------------
! convert satellite number to satellite system
! args   : integer*4    sat       I   satellite number (1-MAXSAT)
!          integer*4    prn       O   satellite prn/slot number (NULL: no output)
!          integer*4    sys       O
! return : satellite system (SYS_GPS,SYS_GLO,...)
! ----------------------------------------------------------------------------
subroutine satsys(sat, prn, sys)
implicit none 
integer*4,intent(in)  :: sat
integer*4,intent(out) :: prn, sys
sys=SYS_NONE
if(sat<=0 .or. MAXSAT<sat)then
    prn=0
elseif (sat<=NSATGPS)then
    sys=SYS_GPS; prn=sat; prn=prn+MINPRNGPS-1
elseif ((sat-NSATGPS)<=NSATGLO)then
    sys=SYS_GLO; prn=sat-NSATGPS; prn=prn+MINPRNGLO-1
elseif ((sat-NSATGPS-NSATGLO)<=NSATGAL)then
	sys=SYS_GAL; prn=sat-NSATGPS-NSATGLO; prn=prn+MINPRNGAL-1
elseif ((sat-NSATGPS-NSATGLO-NSATGAL)<=NSATQZS)then
	sys=SYS_QZS; prn=sat-NSATGPS-NSATGLO-NSATGAL; prn=prn+MINPRNQZS-1; 
elseif ((sat-NSATGPS-NSATGLO-NSATGAL-NSATQZS)<=NSATCMP)then
	sys=SYS_CMP; prn=sat-NSATGPS-NSATGLO-NSATGAL-NSATQZS; prn=prn+MINPRNCMP-1; 
elseif ((sat-NSATGPS-NSATGLO-NSATGAL-NSATQZS-NSATCMP)<=NSATIRN)then
	sys=SYS_IRN; prn=sat-NSATGPS-NSATGLO-NSATGAL-NSATQZS-NSATCMP; prn=prn+MINPRNIRN-1; 
elseif ((sat-NSATGPS-NSATGLO-NSATGAL-NSATQZS-NSATCMP-NSATIRN)<=NSATLEO)then
	sys=SYS_LEO; prn=sat-NSATGPS-NSATGLO-NSATGAL-NSATQZS-NSATCMP-NSATIRN; prn=prn+MINPRNLEO-1; 
elseif ((sat-NSATGPS-NSATGLO-NSATGAL-NSATQZS-NSATCMP-NSATIRN-NSATLEO)<=NSATSBS)then
	sys=SYS_SBS; prn=sat-NSATGPS-NSATGLO-NSATGAL-NSATQZS-NSATCMP-NSATIRN-NSATLEO; prn=prn+MINPRNSBS-1; 
else
    prn=0
endif
end subroutine

! satellite number to satellite id --------------------------------------------
! convert satellite number to satellite id
! args   : integer*4    sat       I   satellite number
!          char   *id             O   satellite id (Gnn,Rnn,Enn,Jnn,Cnn,Inn or nnn)
! return : none
! ----------------------------------------------------------------------------
subroutine satno2id(sat, id)
implicit none
integer*4,intent(in) :: sat
character(*),intent(out) :: id
integer*4 prn,sys
call satsys(sat, prn, sys)

select case(sys)
case(SYS_GPS)
    write(id,"('G',I2.2)") prn-MINPRNGPS+1
case(SYS_GLO)
    write(id,"('R',I2.2)") prn-MINPRNGLO+1
case(SYS_GAL)
    write(id,"('E',I2.2)") prn-MINPRNGAL+1
case(SYS_QZS)
    write(id,"('J',I2.2)") prn-MINPRNQZS+1
case(SYS_CMP)
    write(id,"('C',I2.2)") prn-MINPRNCMP+1
case(SYS_IRN)
    write(id,"('I',I2.2)") prn-MINPRNIRN+1
case(SYS_LEO)
    write(id,"('L',I2.2)") prn-MINPRNLEO+1
case(SYS_SBS)
    write(id,"(I3.3)") prn
case default
    id=''
end select
end subroutine

! satellite system+prn/slot number to satellite number ------------------------
! convert satellite system+prn/slot number to satellite number
! args   : integer*4    sys       I   satellite system (SYS_GPS,SYS_GLO,...)
!          integer*4    prn       I   satellite prn/slot number
! return : satellite number (0:error)
! ----------------------------------------------------------------------------
integer*4 function satno(sys, prn)
implicit none
integer*4,intent(in) :: sys, prn
if (prn<=0)then
    satno=0; return
endif
select case(sys)
case(SYS_GPS)
    if (prn<MINPRNGPS .or. MAXPRNGPS<prn)then
        satno=0; return
    else
        satno=prn-MINPRNGPS+1; return
    endif
case(SYS_GLO)
    if (prn<MINPRNGLO .or. MAXPRNGLO<prn)then
        satno=0; return
    else
        satno=NSATGPS+prn-MINPRNGLO+1; return
    endif
case(SYS_GAL)
    if (prn<MINPRNGAL .or. MAXPRNGAL<prn)then
        satno=0; return
    else
        satno=NSATGPS+NSATGLO+prn-MINPRNGAL+1; return
    endif
case(SYS_QZS)
    if (prn<MINPRNQZS .or. MAXPRNQZS<prn)then
        satno=0; return
    else
        satno=NSATGPS+NSATGLO+NSATGAL+prn-MINPRNQZS+1; return
    endif
case(SYS_CMP)
    if (prn<MINPRNCMP .or. MAXPRNCMP<prn)then
        satno=0; return
    else
        satno=NSATGPS+NSATGLO+NSATGAL+NSATQZS+prn-MINPRNCMP+1; return
    endif
case(SYS_IRN)
    if (prn<MINPRNIRN .or. MAXPRNIRN<prn)then
        satno=0; return
    else
        satno=NSATGPS+NSATGLO+NSATGAL+NSATQZS+NSATCMP+prn-MINPRNIRN+1; return
    endif
case(SYS_LEO)
    if (prn<MINPRNLEO .or. MAXPRNLEO<prn)then
        satno=0; return
    else
        satno=NSATGPS+NSATGLO+NSATGAL+NSATQZS+NSATCMP+NSATIRN+prn-MINPRNLEO+1; return
    endif
case(SYS_SBS)
    if (prn<MINPRNSBS .or. MAXPRNSBS<prn)then
        satno=0; return
    else
        satno=NSATGPS+NSATGLO+NSATGAL+NSATQZS+NSATCMP+NSATIRN+NSATLEO+prn-MINPRNSBS+1; return
    endif
case default
    satno=0; return
end select
end function

! satellite id to satellite number --------------------------------------------
!* * convert satellite id to satellite number
!* * args   : char   *id       I   satellite id (nn,Gnn,Rnn,Enn,Jnn,Cnn,Inn or Snn)
!* * return : satellite number (0: error)
!* * notes  : 120-142 and 193-199 are also recognized as sbas and qzss
!* *-----------------------------------------------------------------------------
integer*4 function satid2no(id)
implicit none
character(*),intent(in) :: id
integer*4 sys,prn
character(1) code
integer*4 err1, err2
read(id(1:2),"(I2)",iostat=err1) prn
if(err1==0)then
    if (MINPRNGPS<=prn .and. prn<=MAXPRNGPS)then
        sys=SYS_GPS
    elseif (MINPRNSBS<=prn .and. prn<=MAXPRNSBS)then
        sys=SYS_SBS
    elseif (MINPRNQZS<=prn .and. prn<=MAXPRNQZS)then
        sys=SYS_QZS
    else
        satid2no=0; return
    endif
    satid2no=satno(sys,prn); return
endif
read(id(1:3),"(A1,I2)",iostat=err2) code, prn
if(err2/=0)then
    satid2no=0; return
endif
select case(code)
case('G')
    sys=SYS_GPS; prn=prn+MINPRNGPS-1
case('R')
    sys=SYS_GLO; prn=prn+MINPRNGLO-1
case('E')
    sys=SYS_GAL; prn=prn+MINPRNGAL-1
case('J')
    sys=SYS_QZS; prn=prn+MINPRNQZS-1
case('C')
    sys=SYS_CMP; prn=prn+MINPRNCMP-1
case('I')
    sys=SYS_IRN; prn=prn+MINPRNIRN-1
case('L')
    sys=SYS_LEO; prn=prn+MINPRNLEO-1
case('S')
    sys=SYS_SBS; prn=prn+100
case default
    satid2no=0; return
end select
satid2no=satno(sys,prn)
end function

! test excluded satellite -----------------------------------------------------
! test excluded satellite
! args   : integer*4    sat       I   satellite number
!          integer*4    svh       I   sv health flag
!          prcopt_t *opt          I   processing options (NULL: not used)
! return : status (1:excluded,0:not excluded)
!-----------------------------------------------------------------------------
integer*4 function satexclude(sat, var, svh, opt, isopt)
implicit none
integer*4, intent(in) :: sat, svh
real*8, intent(in) :: var
type(prcopt_t), intent(in) :: opt
logical*1, intent(in) :: isopt
integer*4 prn, sys
call satsys(sat,prn,sys)
if (svh<0)then  ! ephemeris unavailable 
    satexclude=1; return
endif
!if(associated(opt))then
if(isopt)then
    if (opt%exsats(sat)==1)then  ! excluded satellite 
        satexclude=1; return
    endif
    if (opt%exsats(sat)==2)then  ! included satellite 
        satexclude=0; return
    endif
    if (and(sys,opt%navsys)==0)then  ! unselected sat sys 
        satexclude=1; return
    endif
endif
if (svh/=0)then
    satexclude=1; return
endif
if (var>MAX_VAR_EPH)then
    satexclude=1; return
endif
satexclude=0
end function

! test SNR mask
! args   : integer*4 base   I   rover or base-station (0:rover,1:base station)
!          integer*4 freq   I   frequency (0:L1,1:L2,2:L3,...)
!          real*8 el        I   elevation angle (rad)
!          real*8 snr       I   C/N0 (dBHz)
!          snrmask_t *mask  I   SNR mask
! return : status (1:masked,0:unmasked)
!-----------------------------------------------------------------------------
integer*4 function testsnr(base, freq, el, snr, mask)
implicit none
integer*4, intent(in) :: base, freq
real*8, intent(in):: el, snr
type(snrmask_t), intent(in) :: mask
real*8 minsnr,a
integer*4 i
if(mask%ena(base+1)==0 .or. freq<0 .or. freq>=NFREQ)then
    testsnr=0; return
endif
a=(el*R2D+5.d0)/10.d0
i=floor(a); a=a-i
if (i<1)then
    minsnr=mask%mask(freq,1)
elseif (i>8)then
    minsnr=mask%mask(freq,9)
else
    minsnr=(1.d0-a)*mask%mask(freq,i)+a*mask%mask(freq,i+1)
endif
testsnr=icond(snr<minsnr,1,0)
end function

! obs type string to obs code -------------------------------------------------
!* * convert obs code type string to obs code
!* * args   : char   *str      I      obs code string ('1C','1P','1Y',...)
!* *          integer*4 *freq  IO     frequency (1:L1,2:L2,3:L5,4:L6,5:L7,6:L8,0:err)
!* *                               (NULL: no output)
!* * return : obs code (CODE_???)
!* * notes  : obs codes are based on reference (6) and qzss extension
!* *---------------------------------------------------------------------------
subroutine obs2code(obs, freq, code)
implicit none
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

! obs code to obs code string -------------------------------------------------
!* * convert obs code to obs code string
!* * args   : character(1) code I obs code (CODE_???)
!* *          integer*4 *freq   IO     frequency (NULL: no output)
!* * (1:L1/E1, 2:L2/B1, 3:L5/E5a/L3, 4:L6/LEX/B3,
!* *  5:E5b/B2, 6:E5(a+b), 7:S)
!* * return : obs code string ('1C','1P','1P',...)
!* * notes  : obs codes are based on reference (6) and qzss extension
!* *---------------------------------------------------------------------------
subroutine code2obs(code, freq, obs)
implicit none
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

! set code priority -----------------------------------------------------------
! set code priority for multiple codes in a frequency
! args   : integer*4 sys  I     system (or of SYS_???)
!          integer*4 freq I     frequency (1:L1,2:L2,3:L5,4:L6,5:L7,6:L8,7:L9)
!          char   *pri    I     priority of codes (series of code characters)
!                               (higher priority precedes lower)
! return : none
!-----------------------------------------------------------------------------
subroutine setcodepri(sys, freq, pri)
implicit none
integer*4, intent(in) :: sys, freq
character(*), intent(in) :: pri
if (freq<=0 .or. MAXFREQ<freq) return
if (and(sys,SYS_GPS)/=0) codepris(1,freq)=pri
if (and(sys,SYS_GLO)/=0) codepris(2,freq)=pri
if (and(sys,SYS_GAL)/=0) codepris(3,freq)=pri
if (and(sys,SYS_QZS)/=0) codepris(4,freq)=pri
if (and(sys,SYS_SBS)/=0) codepris(5,freq)=pri
if (and(sys,SYS_CMP)/=0) codepris(6,freq)=pri
if (and(sys,SYS_IRN)/=0) codepris(7,freq)=pri
end subroutine

! get code priority -----------------------------------------------------------
! get code priority for multiple codes in a frequency
! args   : integer*4 sys     I     system (SYS_???)
!          character(1) code I     obs code (CODE_???)
!          char   *opt       I     code options (NULL:no option)
! return : priority (15:highest-1:lowest,0:error)
!-----------------------------------------------------------------------------
integer*4 function getcodepri(sys, code, opt)
implicit none
integer*4, intent(in) :: sys, code
character(*), intent(in) :: opt
character(256) p, optstr
character(2) obs
character(8) str
integer*4 i, j, info, posi
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

! new matrix ------------------------------------------------------------------
! allocate memory of matrix 
! args   : integer*4    n,m       I   number of rows and columns of matrix
! return : matrix pointer (if n<=0 or m<=0, return NULL)
!-----------------------------------------------------------------------------
subroutine mat(n, m, matrix)
implicit none
integer*4, intent(in) :: n,m
real*8, intent(out) :: matrix(n,m)
matrix=0.d0
end subroutine

! new integer matrix ----------------------------------------------------------
! allocate memory of integer matrix 
! args   : integer*4    n,m       I   number of rows and columns of matrix
! return : matrix pointer (if n<=0 or m<=0, return NULL)
!-----------------------------------------------------------------------------
subroutine imat(n, m, matrix)
implicit none
integer*4, intent(in) :: n,m
integer*4, intent(out) :: matrix(n,m)
matrix=0
end subroutine

! zero matrix -----------------------------------------------------------------
! generate new zero matrix
! args   : integer*4    n,m       I   number of rows and columns of matrix
! return : matrix pointer (if n<=0 or m<=0, return NULL)
!-----------------------------------------------------------------------------
subroutine zeros(n, m, matrix)
implicit none
integer*4, intent(in) :: n,m
real*8, intent(out) :: matrix(n,m)
matrix=0.d0
end subroutine

! identity matrix -------------------------------------------------------------
! generate new identity matrix
! args   : integer*4    n         I   number of rows and columns of matrix
! return : matrix pointer (if n<=0, return NULL)
!-----------------------------------------------------------------------------
subroutine eye(n, matrix)
integer*4, intent(in) :: n
real*8, intent(out) :: matrix(n,n)
integer*8 i,j
matrix=0.d0
forall(i=1:n,j=1:n,i==j)
    matrix(i,j)=1.d0
end forall
end subroutine

! inner product ---------------------------------------------------------------
! inner product of vectors
! args   : real*8 *a,*b     I   vector a,b (n x 1)
!          integer*4 n      I   size of vector a,b
! return : a'*b
!-----------------------------------------------------------------------------
real*8 function dot(a, b, n)
implicit none
integer*4, intent(in) :: n
real*8, intent(in) :: a(n),b(n)
real*8 c(n)
c=a*b
dot=sum(c)
end function
    
! euclid norm -----------------------------------------------------------------
! euclid norm of vector
! args   : real*8 *a        I   vector a (n x 1)
!          integer*4 n      I   size of vector a
! return :  .or.  a  .or. 
!-----------------------------------------------------------------------------
real*8 function norm(a, n)
implicit none
integer*4, intent(in) :: n
real*8, intent(in) :: a(n)
norm=dsqrt(dot(a,a,n))
end function

! outer product of 3d vectors -------------------------------------------------
! outer product of 3d vectors 
! args   : real*8 *a,*b     I   vector a,b (3 x 1)
!          real*8 *c        O   outer product (a x b) (3 x 1)
! return : none
!-----------------------------------------------------------------------------
subroutine cross3(a, b, c)
implicit none
real*8, intent(in) :: a(3),b(3)
real*8, intent(out) :: c(3)
c(1)=a(2)*b(3)-a(3)*b(2)
c(2)=a(3)*b(1)-a(1)*b(3)
c(3)=a(1)*b(2)-a(2)*b(1)
end subroutine

! normalize 3d vector ---------------------------------------------------------
! normalize 3d vector
! args   : real*8 *a        I   vector a (3 x 1)
!          real*8 *b        O   normlized vector (3 x 1)  .or.  b  .or.  = 1
! return : status (1:ok,0:error)
!-----------------------------------------------------------------------------
subroutine normv3(a, b, stat)
implicit none
real*8, intent(in) :: a(3)
real*8, intent(out) :: b(3)
integer*4, intent(out) :: stat
real*8 r
r=norm(a,3)
if(r<=0.d0)then
    stat=0; b=0.d0; return
endif
b(1)=a(1)/r
b(2)=a(2)/r
b(3)=a(3)/r
stat=1
end subroutine

! copy matrix ----------------------------------------------------------------
! copy matrix
! args   : real*8 *A        O   destination matrix A (n x m)
!          real*8 *B        I   source matrix B (n x m)
!          integer*4 n,m    I   number of rows and columns of matrix
! return : none
!-----------------------------------------------------------------------------
subroutine matcpy(A, B, n, m)
implicit none
real*8, intent(out) :: A(n,m)
real*8, intent(in) :: B(n,m)
integer*4, intent(in) :: n,m
A=B
end subroutine

subroutine mattran(n, m, A, B)  ! A(n,m) B(m,n)
implicit none
integer*4, intent(in) :: n,m
real*8, intent(in) :: A(n,m)
real*8, intent(out) :: B(m,n)
B=reshape(A,(/m,n/),order=(/2,1/))
end subroutine

subroutine matmul2(tr, n, k, m, A, B, C)
implicit none
character(2), intent(in) :: tr
integer*4, intent(in) :: n, k, m
real*8, intent(in)  :: A(n,m),B(m,k)
real*8, intent(out) :: C(n,k)
real*8 A2(n,m),B2(m,k)
A2=A; B2=B; C=0.d0

if(tr(1:1)=='N')then
    if(tr(2:2)=='T') B2=reshape(B,(/m,k/),order=(/2,1/))
elseif(tr(1:1)=='T')then
    A2=reshape(A,(/n,m/),order=(/2,1/))
    if(tr(2:2)=='T') B2=reshape(B,(/m,k/),order=(/2,1/))
endif
C=matmul(A2,B2)
end subroutine

!gauss-jordan method to compute the inverse of matrix a
!a -- is the n*n matrix ,it's input and output(inverse of a)
!n--- is the n ,input
!l--- if the matrix a is dsingular ,l=0 ,else l /=0
!is,js-----is a matrix is(n),and js(n)
subroutine matinv(a, n, l)
implicit none
integer*4, intent(in) :: n
integer*4, intent(out) :: l
real*8    a(n,n),t,d
integer*4 i,j,k,is(n),js(n)
l=0
do k=1,n
    d=0.d0
	do i=k,n
        do j=k,n
            if (abs(a(i,j)).gt.d) then
                d=abs(a(i,j))
                is(k)=i
                js(k)=j
            endif
        enddo
    enddo
    if (d<=1d-20) then
        l=-1; !write(*,"(1x,'err**not inv')")
        return
    endif

	do j=1,n
        t=a(k,j)
        a(k,j)=a(is(k),j)
        a(is(k),j)=t
    enddo
	do i=1,n
        t=a(i,k)
        a(i,k)=a(i,js(k))
        a(i,js(k))=t
    enddo
	a(k,k)=1.d0/a(k,k)
	do j=1,n
        if (j.ne.k) then
            a(k,j)=a(k,j)*a(k,k)
        endif
    enddo
	do i=1,n
        if (i.ne.k) then
            do j=1,n
                if (j.ne.k) then
                    a(i,j)=a(i,j)-a(i,k)*a(k,j)
                endif
            enddo
        endif
    enddo
	do i=1,n
        if (i.ne.k) then
            a(i,k)=-a(i,k)*a(k,k)
        endif
    enddo
enddo
do k=n,1,-1
	do j=1,n
        t=a(k,j)
        a(k,j)=a(js(k),j)
        a(js(k),j)=t
    enddo
	do i=1,n
        t=a(i,k)
        a(i,k)=a(i,is(k))
        a(i,is(k))=t
    enddo
enddo
end subroutine

!! solve linear equation -----------------------------------------------------
subroutine solve(tr, A, Y, n, m, X,stat)
implicit none
character(*), intent(in) :: tr
integer*4, intent(in) :: n,m
real*8, intent(in) :: A(n,n),Y(n,m)
real*8, intent(out) :: X(n,m)
integer*4, intent(out) :: stat
real*8 B(n,n)
integer*4 info
B=A
call matinv(B,n,info)
if(info==0)then
    if(tr(1:1)=='N') call matmul2('NN',n,m,n,B,Y,X)
    if(tr(1:1)/='N') call matmul2('TN',n,m,n,B,Y,X)
endif
stat=info
end subroutine

! least square estimation ----------------------------------------------------
! least square estimation by solving normal equation (x=(A*A')^-1*A*y)
! args   : real*8 *A        I   transpose of (weighted) design matrix (n x m)
!          real*8 *y        I   (weighted) measurements (m x 1)
!          integer*4 n,m    I   number of parameters and measurements (n<=m)
!          real*8 *x        O   estmated parameters (n x 1)
!          real*8 *Q        O   esimated parameters covariance matrix (n x n)
! return : status (0:ok,0>:error)
! notes  : for weighted least square, replace A and y by A*w and w*y (w=W^(1/2))
!          matirix stored by column-major order (fortran convention)
!-----------------------------------------------------------------------------
subroutine lsq(A, y, n, m, x, Q, stat)
implicit none
integer*4, intent(in) :: n,m
real*8, intent(in) :: A(n,m),y(m,1)
real*8, intent(out) :: Q(n,n),x(n,1)
integer*4, intent(out) :: stat
real*8 Ay(n,1)
integer*4 info
Q=0.d0
x=0.d0
if(m<n)then
    stat=-1; return
endif
Ay=0.d0
call matmul2('NN',n,1,m,A,y,Ay)
call matmul2('NT',n,n,m,A,A,Q)
call matinv(Q,n,info)
if(info==0) call matmul2('NN',n,1,n,Q,Ay,x)
stat=info
end subroutine

! string to number -----------------------------------------------------------
! convert substring in string to number
! args   : char   *s        I   string ('... nnn.nnn ...')
!          integer*4 i,n    I   substring position and width
! return : converted number (0.d0:error)
!-----------------------------------------------------------------------------
real*8 function str2num(s,i,n)
implicit none
character(*), intent(in) :: s
integer*4, intent(in) :: i,n
integer*4 err
str2num=0.d0
if(i<1 .or. len_trim(s)<i+n-1 .or. n<1) return
read(s(i:i+n-1),*,iostat=err) str2num
if(err/=0) str2num=0.d0
end function

! convert ymd to mjd ---------------------------------------------------------
integer*4 function ymd2mjd(ep)
implicit none
integer*4, intent(in) :: ep(3)
integer*4 year, mon, day
real*8 mjd
year=ep(1); mon=ep(2); day=ep(3)
if(year<100) year=year+2000
if(mon<=2)then
    mon=mon+12
    year=year-1
endif
mjd=year*365.25d0-dmod(year*365.25d0,1.d0)-679006
mjd=mjd+int(30.6001d0*(mon+1))+2-int(year/100.d0)+int(year/400.d0)+day
ymd2mjd=int(mjd)
end function

! convert calendar day/time to time ------------------------------------------
! convert calendar day/time to gtime_t struct
! args   : real*8 *ep       I   day/time {year,month,day,hour,min,sec}
! return : gtime_t struct
! notes  : proper in 1970-2037 or 1970-2099 (64bit time_t)
!-----------------------------------------------------------------------------
type(gtime_t) function epoch2time(ep)
implicit none
real*8, intent(in) :: ep(6)
integer*4 :: doy(12)=(/1,32,60,91,121,152,182,213,244,274,305,335/)
type(gtime_t) :: mytime
integer*4 days,sec,year,mon,day
mytime=gtime_t(0,0.d0)
year=int(ep(1)); mon=int(ep(2)); day=int(ep(3))
if (year<1970 .or. 2099<year .or. mon<1 .or. 12<mon)then
    epoch2time=mytime; return
endif
days=(year-1970)*365+(year-1969)/4+doy(mon)+day-2+icond(mod(year,4)==0.and.mon>=3,1,0)
sec=floor(ep(6))
mytime%time=days*86400+int(ep(4))*3600+int(ep(5))*60+sec
mytime%sec=ep(6)-sec
epoch2time=mytime
end function

! string to time -------------------------------------------------------------
! convert substring in string to gtime_t struct
! args   : char   *s        I   string ('... yyyy mm dd hh mm ss ...')
!          integer*4 i,n    I   substring position and width
!          gtime_t *t       O   gtime_t struct
! return : status (0:ok,0>:error)
!-----------------------------------------------------------------------------
subroutine str2time(s, i, n, t, stat)
implicit none
character(*), intent(in) :: s
integer*4, intent(in) :: i,n
type(gtime_t), intent(out) :: t
integer*4, intent(out) :: stat
real*8 ep(6)
integer*4 err
t=gtime_t(0,0.d0); stat=0
if (i<1 .or. len_trim(s)<i+n-1 .or. n<1)then
    stat=-1; return
endif
read(s(i:i+n-1),*,iostat=err) ep(1:6)
if(err/=0)then
    stat=-1; return
endif
if (ep(1)<100.d0) ep(1)=ep(1)+rcond(ep(1)<80.d0,2000.d0,1900.d0)
t=epoch2time(ep)
end subroutine

! time to calendar day/time --------------------------------------------------
! convert gtime_t struct to calendar day/time
! args   : gtime_t t        I   gtime_t struct
!          real*8 *ep       O   day/time {year,month,day,hour,min,sec}
! return : none
! notes  : proper in 1970-2037 or 1970-2099 (64bit time_t)
!-----------------------------------------------------------------------------
subroutine time2epoch(t, ep)
implicit none
type(gtime_t), intent(in) :: t
real*8, intent(out) :: ep(6)
integer*4 days,sec,mon,day
integer*4, parameter :: mday(48)=(/&
    31,28,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31,&
    31,29,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31/)
days=t%time/86400
sec=int(t%time-days*86400)
day=mod(days,1461)
do mon=1,48,1
    if (day>=mday(mon))then
        day=day-mday(mon)
    else
        exit
    endif
enddo
ep(1)=1970+days/1461*4+(mon-1)/12; ep(2)=mod(mon-1,12)+1; ep(3)=day+1
ep(4)=sec/3600; ep(5)=mod(sec,3600)/60; ep(6)=mod(sec,60)+t%sec
end subroutine

! gps time to time -----------------------------------------------------------
! convert week and tow in gps time to gtime_t struct
! args   : integer*4 week   I   week number in gps time
!          real*8 sec       I   time of week in gps time (s)
! return : gtime_t struct
!-----------------------------------------------------------------------------
type(gtime_t) function gpst2time(week, sec)
implicit none
integer*4, intent(in) :: week
real*8, intent(in) :: sec
type(gtime_t) t
real*8 sec2
sec2=sec; t=epoch2time(gpst0)
if (sec2<-1d9 .or. 1d9<sec2) sec2=0.d0
t%time=t%time+86400*7*week+int(sec2)
t%sec=sec2-int(sec2)
gpst2time=t
end function

! time to gps time -----------------------------------------------------------
! convert gtime_t struct to week and tow in gps time
! args   : gtime_t t        I   gtime_t struct
!          integer*4 *week  IO  week number in gps time (NULL: no output)
! return : time of week in gps time (s)
!-----------------------------------------------------------------------------
subroutine time2gpst(t, week, gpst)
implicit none
type(gtime_t), intent(in) :: t
integer*4, intent(out) :: week
real*8, intent(out) :: gpst
type(gtime_t) t0
integer*4 sec
t0=epoch2time(gpst0); sec=t%time-t0%time
week=sec/(86400*7)
gpst=sec-week*86400*7.d0+t%sec
end subroutine

! convert time to mjd and sod --------------------------------------------------------
subroutine time2mjd(t, mjd, sod)
implicit none
type(gtime_t), intent(in) :: t
integer*4, intent(out) :: mjd
real*8, intent(out) :: sod
integer*4 week
real*8 gpst,ep(6)
call time2epoch(t,ep)
mjd=ymd2mjd(int(ep(1:3)))
call time2gpst(t,week,gpst)
sod=dmod(gpst,86400.d0)
end subroutine

! add time -------------------------------------------------------------------
! add time to gtime_t struct
! args   : gtime_t t        I   gtime_t struct
!          real*8 sec       I   time to add (s)
! return : gtime_t struct (t+sec)
!-----------------------------------------------------------------------------
type(gtime_t) function timeadd(t, sec)
implicit none
type(gtime_t), intent(in) :: t
real*8, intent(in) :: sec
type(gtime_t) t1
real*8 tt
t1=t
t1%sec=t1%sec+sec; tt=floor(t1%sec); t1%time=t1%time+int(tt); t1%sec=t1%sec-tt
timeadd=t1
end function

! time difference ------------------------------------------------------------
! difference between gtime_t structs
! args   : gtime_t t1,t2    I   gtime_t structs
! return : time difference (t1-t2) (s)
!-----------------------------------------------------------------------------
real*8 function timediff(t1,t2)
implicit none
type(gtime_t), intent(in) :: t1,t2
timediff=t1%time-t2%time+t1%sec-t2%sec
end function

! get current time in utc ----------------------------------------------------
! get current time in utc
! args   : none
! return : current time in utc
!-----------------------------------------------------------------------------
type(gtime_t) function timeget()
implicit none
character(8) date
character(10) time_s
character(5) zone
integer*4 values(8)
real*8 :: ep(6)
type(gtime_t) time
ep=0.d0
call date_and_time(date,time_s,zone,values)
ep(1:3)=values(1:3)  ! year, month, day
ep(4:6)=values(5:7)  ! h, m, s
ep(6)=ep(6)+values(8)*1d-3  ! ms
time=epoch2time(ep)
time=timeadd(time,-values(4)*60.d0)
time=timeadd(time,timeoffset_)
timeget=time
end function

! set current time in utc ----------------------------------------------------
! set current time in utc
! args   : gtime_t          I   current time in utc
! return : none
! notes  : just set time offset between cpu time and current time
!          the time offset is reflected to only timeget()
!          not reentrant
!-----------------------------------------------------------------------------
subroutine timeset(t)
implicit none
type(gtime_t), intent(in) :: t
timeoffset_=timeoffset_+timediff(t,timeget())
end subroutine

! gpstime to utc -------------------------------------------------------------
! convert gpstime to utc considering leap seconds
! args   : gtime_t t        I   time expressed in gpstime
! return : time expressed in utc
! notes  : ignore slight time offset under 100 ns
!-----------------------------------------------------------------------------
type(gtime_t) function gpst2utc(t)
implicit none
type(gtime_t), intent(in) :: t
type(gtime_t) :: tu
integer*4 :: i
tu=gtime_t(0,0.d0)
i=1
do while(.true.)
    if(leaps(i,1)<=0)exit
    tu=timeadd(t,leaps(i,7))
    if(timediff(tu,epoch2time(leaps(i,1:6)))>=0.d0)then
        gpst2utc=tu; return
    endif
    i=i+1
enddo
gpst2utc=t
end function

! utc to gpstime -------------------------------------------------------------
! convert utc to gpstime considering leap seconds
! args   : gtime_t t        I   time expressed in utc
! return : time expressed in gpstime
! notes  : ignore slight time offset under 100 ns
!-----------------------------------------------------------------------------
type(gtime_t) function utc2gpst(t)
implicit none
type(gtime_t), intent(in) :: t
integer*4 :: i
i=1
do while(.true.)
    if(leaps(i,1)<=0) exit
    if(timediff(t,epoch2time(leaps(i,1:6)))>=0.d0)then
        utc2gpst=timeadd(t,-leaps(i,7)); return
    endif
    i=i+1
enddo
utc2gpst=t
end function

! time to day and sec -------------------------------------------------------
subroutine time2sec(time, day, sec)
implicit none
type(gtime_t), intent(in) :: time
type(gtime_t), intent(out) :: day
real*8, intent(out) :: sec
real*8 ep(6)
call time2epoch(time,ep)
sec=ep(4)*3600.d0+ep(5)*60.d0+ep(6)
ep(4:6)=0.d0
day=epoch2time(ep)
end subroutine

! utc to gmst ----------------------------------------------------------------
! convert utc to gmst (Greenwich mean sidereal time)
! args   : gtime_t t        I   time expressed in utc
!          real*8 ut1_utc   I   UT1-UTC (s)
! return : gmst (rad)
!-----------------------------------------------------------------------------
real*8 function utc2gmst(t, ut1_utc)
implicit none
type(gtime_t), intent(in) :: t
real*8, intent(in) :: ut1_utc
real*8 :: ep2000(6)
type(gtime_t) tut,tut0
real*8 ut,t1,t2,t3,gmst0,gmst
ep2000=(/2000,1,1,12,0,0/)
tut=timeadd(t,ut1_utc)
call time2sec(tut,tut0,ut)
t1=timediff(tut0,epoch2time(ep2000))/86400.d0/36525.d0
t2=t1*t1; t3=t2*t1
gmst0=24110.54841d0+8640184.812866d0*t1+0.093104d0*t2-6.2d-6*t3
gmst=gmst0+1.002737909350795d0*ut
utc2gmst=dmod(gmst,86400.d0)*PI/43200.d0
end function

! time to string -------------------------------------------------------------
! convert gtime_t struct to string
! args   : gtime_t t        I   gtime_t struct
!          char   *s        O   string ('yyyy/mm/dd hh:mm:ss.ssss')
!          integer*4 n      I   number of decimals
! return : none
!-----------------------------------------------------------------------------
subroutine time2str(t,s,n)
implicit none
type(gtime_t), intent(in) :: t
integer*4, intent(in) :: n
character(*), intent(out) :: s
real*8 ep(6)
type(gtime_t) t1
integer*4 n1,f1,f2
character(60) fmtstr,f6

t1=t; n1=n
if(n1<0) n1=0; if(n1>12) n1=12
if(1.d0-t1%sec<0.5d0/(10.d0**n1))then
    t1%time=t1%time+1; t1%sec=0.d0
endif
call time2epoch(t1,ep)
write(f6,"(I2'.'I2)") icond(n1<=0,2,n1+3),icond(n1<=0,0,n1)
fmtstr="(I4.4'/'I2.2'/'I2.2' 'I2.2':'I2.2':'f"//trim(f6)//")"
write(s,fmtstr) int(ep(1:5)),ep(6)
end subroutine

! get time string -------------------------------------------------------------
! get time string
! args   : gtime_t t        I   gtime_t struct
!          integer*4 n      I   number of decimals
! return : time string
! notes  : not reentrant, do not use multiple in a function
!-----------------------------------------------------------------------------
!extern char *time_str(gtime_t t, integer*4 n)
!{
!    static char buff(64)
!    time2str(t,buff,n)
!    return buff
!}

! time to day of year ---------------------------------------------------------
! convert time to day of year
! args   : gtime_t t        I   gtime_t struct
! return : day of year (days)
!-----------------------------------------------------------------------------
real*8 function time2doy(t)
implicit none
type(gtime_t), intent(in) :: t
real*8 ep(6)
call time2epoch(t,ep)
ep(2:3)=1.d0
ep(4:6)=0.d0
time2doy=timediff(t,epoch2time(ep))/86400.d0+1.d0
end function

! adjust gps week number ------------------------------------------------------
! adjust gps week number using cpu time
! args   : integer*4   week       I   not-adjusted gps week number
! return : adjusted gps week number
!-----------------------------------------------------------------------------
integer*4 function adjgpsweek(week)
implicit none
integer*4, intent(in) :: week
integer*4 w
real*8 gpst
call time2gpst(utc2gpst(timeget()),w,gpst)
if(w<1560) w=1560
adjgpsweek=week+(w-week+512)/1024*1024
end function

! convert degree to deg-min-sec -----------------------------------------------
! convert degree to degree-minute-second
! args   : real*8 deg       I   degree
!          real*8 *dms      O   degree-minute-second {deg,min,sec}
!          integer*4 ndec   I   number of decimals of second
! return : none
!-----------------------------------------------------------------------------
subroutine deg2dms(deg, dms, ndec)
implicit none
real*8, intent(in) :: deg
real*8, intent(out) :: dms(3)
integer*4, intent(in) :: ndec
real*8 sign, a, unit
sign=rcond(deg<0.d0,-1.d0,1.d0)
a=dabs(deg)
unit=0.1d0**ndec
dms(1)=floor(a); a=(a-dms(1))*60.d0
dms(2)=floor(a); a=(a-dms(2))*60.d0
dms(3)=floor(a/unit+0.5)*unit
if(dms(3)>=60.d0)then
    dms(3)=0.d0
    dms(2)=dms(2)+1.d0
    if(dms(2)>=60.d0)then
        dms(2)=0.d0
        dms(1)=dms(1)+1.d0
    endif
endif
dms(1)=dms(1)*sign
end subroutine

! convert deg-min-sec to degree -----------------------------------------------
! convert degree-minute-second to degree
! args   : real*8 *dms      I   degree-minute-second {deg,min,sec}
! return : degree
!-----------------------------------------------------------------------------
real*8 function dms2deg(dms)
implicit none
real*8, intent(in) :: dms(3)
real*8 sign
sign=rcond(dms(1)<0.d0,-1.d0,1.d0)
dms2deg=sign*(dabs(dms(1))+dms(2)/60.d0+dms(3)/3600.d0)
end function

! free observation data -------------------------------------------------------
! free memory for observation data
! args   : obs_t *obs    IO     observation data
! return : none
!-----------------------------------------------------------------------------
subroutine freeobs(obs)
implicit none
type(obs_t), intent(out) :: obs
deallocate(obs%mydata)
nullify(obs%mydata)
obs%n=0
obs%nmax=0
end subroutine

! free navigation data ---------------------------------------------------------
! free memory for navigation data
! args   : nav_t *nav    IO     navigation data
!          integer*4 opt I      option (or of followings)
!                               (0x01: gps/qzs ephmeris, 0x02: glonass ephemeris,
!                                0x04: sbas ephemeris,   0x08: precise ephemeris,
!                                0x10: precise clock     0x20: almanac,
!                                0x40: tec data)
! return : none
!-----------------------------------------------------------------------------
subroutine freenav(nav, opt)
implicit none
type(nav_t), intent(out) :: nav
integer*4, intent(in) :: opt
if(and(opt,01)/=0)then  !#01
    deallocate(nav%eph); nullify(nav%eph)
    nav%n=0; nav%nmax=0
endif
if(and(opt,02)/=0)then  !#02
    deallocate(nav%geph); nullify(nav%geph)
    nav%ng=0; nav%ngmax=0
endif
if(and(opt,32)/=0)then  !#20
    deallocate(nav%alm); nullify(nav%alm)
    nav%na=0; nav%namax=0
endif
end subroutine

! ecef to local coordinate transfromation matrix ------------------------------
! * * compute ecef to local coordinate transfromation matrix
! * * args   : real*8 *pos      I   geodetic position {lat,lon} (rad)
! * *          real*8 *E        O   ecef to local coord transformation matrix (3x3)
! * * return : none
! * * notes  : matirix stored by column-major order (fortran convention)
! * *-----------------------------------------------------------------------------
subroutine xyz2enu(pos, E)  !?
implicit none
real*8, intent(in) :: pos(3)
real*8, intent(out) :: E(3,3)
real*8 sinp,cosp, sinl,cosl
sinp=dsin(pos(1)); cosp=dcos(pos(1))
sinl=dsin(pos(2)); cosl=dcos(pos(2))
E(1,1)=-sinl;      E(1,2)=cosl;       E(1,3)=0.d0
E(2,1)=-sinp*cosl; E(2,2)=-sinp*sinl; E(2,3)=cosp
E(3,1)=cosp*cosl;  E(3,2)=cosp*sinl;  E(3,3)=sinp
end subroutine

! transform ecef vector to local tangental coordinate -------------------------
! * * transform ecef vector to local tangental coordinate
! * * args   : real*8 *pos      I   geodetic position {lat,lon} (rad)
! * *          real*8 *r        I   vector in ecef coordinate {x,y,z}
! * *          real*8 *e        O   vector in local tangental coordinate {e,n,u}
! * * return : none
! * *-----------------------------------------------------------------------------
subroutine ecef2enu(pos, r, e)  !?
implicit none
real*8, intent(in) :: pos(3)
real*8, intent(in) :: r(3,1)
real*8, intent(out) :: e(3,1)
real*8 E1(3,3)
call xyz2enu(pos,E1)
call matmul2('NN',3,1,3,E1,r,e)
end subroutine

! transform ecef to geodetic postion ------------------------------------------
! * * transform ecef position to geodetic position
! * * args   : real*8 *r        I   ecef position {x,y,z} (m)
! * *          real*8 *pos      O   geodetic position {lat,lon,h} (rad,m)
! * * return : none
! * * notes  : WGS84, ellipsoidal height
! * *-----------------------------------------------------------------------------
subroutine ecef2pos(r, pos)
implicit none
real*8, intent(in) :: r(3)
real*8, intent(out) :: pos(3)
real*8 e2,r2,z,zk,v,sinp
e2=FE_WGS84*(2.d0-FE_WGS84)
r2=dot(r,r,2)
v=RE_WGS84
z=r(3); zk=0.d0
do while(dabs(z-zk)>=1d-4)
    zk=z
    sinp=z/dsqrt(r2+z*z)
    v=RE_WGS84/dsqrt(1.d0-e2*sinp*sinp)
    z=r(3)+v*e2*sinp
enddo
pos(1)=rcond(r2>1d-12,datan(z/dsqrt(r2)),rcond(r(3)>0.d0,PI/2.d0,-PI/2.d0))
pos(2)=rcond(r2>1d-12,datan2(r(2),r(1)),0.d0)
pos(3)=dsqrt(r2+z*z)-v
end subroutine

! transform geodetic to ecef position -----------------------------------------
! transform geodetic position to ecef position
! args   : real*8 *pos      I   geodetic position {lat,lon,h} (rad,m)
!          real*8 *r        O   ecef position {x,y,z} (m)
! return : none
! notes  : WGS84, ellipsoidal height
!-----------------------------------------------------------------------------
subroutine pos2ecef(pos, r)
implicit none
real*8, intent(in) :: pos(:)
real*8, intent(out) :: r(:)
real*8 sinp, cosp, sinl, cosl
real*8 e2, v
sinp=dsin(pos(1)); cosp=dcos(pos(1))
sinl=dsin(pos(2)); cosl=dcos(pos(2))
e2=FE_WGS84*(2.d0-FE_WGS84)
v =RE_WGS84/dsqrt(1.d0-e2*sinp*sinp)
r(1)=(v+pos(3))*cosp*cosl
r(2)=(v+pos(3))*cosp*sinl
r(3)=(v*(1.d0-e2)+pos(3))*sinp
end subroutine

! least common multiple (int) -------------------------------------------------
integer*4 function iminmul(a, b)
implicit none
integer*4, intent(in) :: a, b
integer*4 a1, b1, product, tmp
product=a*b
a1=max(a,b)
b1=min(a,b)
do while (mod(a1,b1)/=0)
    tmp=b1
    b1=mod(a1,b1)
    a1=tmp
enddo
iminmul=product/b1
end function

! init ephemeris ------------------------------------------------------------
subroutine init_eph(eph)
implicit none
type(eph_t), intent(out) :: eph
eph=eph_t(0,0,0,0,0,0,0,0,&
        gtime_t(0,0.d0),gtime_t(0,0.d0),gtime_t(0,0.d0),&
        0d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,&
        0d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,&
        0d0,0.d0,0.d0)
end subroutine

subroutine init_geph(geph)
implicit none
type(geph_t), intent(out) :: geph
geph=geph_t(0,0,0,0,0,0,&
        gtime_t(0,0.d0),gtime_t(0,0.d0),&
        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0)
end subroutine

!! compare ephemeris ---------------------------------------------------------
!integer*4 function cmpeph(q1, q2)
!implicit none
!type(eph_t), intent(in) :: q1, q2
!
!if(q1%ttr%time/=q2%ttr%time)then
!    cmpeph=int(q1%ttr%time-q2%ttr%time)
!else
!    if(q1%toe%time/=q2%toe%time)then
!        cmpeph=int(q1%toe%time-q2%toe%time)
!    else
!        cmpeph=q1%sat-q2%sat
!    endif
!endif
!end function
!
!! compare glonass ephemeris -------------------------------------------------
!integer*4 function cmpgeph(q1, q2)
!implicit none
!type(geph_t), intent(in) :: q1, q2
!if(q1%tof%time/=q2%tof%time)then
!    cmpgeph=int(q1%tof%time-q2%tof%time)
!else
!    if(q1%toe%time/=q2%toe%time)then
!        cmpgeph=int(q1%toe%time-q2%toe%time)
!    else
!        cmpgeph=q1%sat-q2%sat
!    endif
!endif
!end function

!! sort and unique ephemeris -------------------------------------------------
!subroutine uniqeph(nav)
!implicit none
!type(nav_t), intent(out) :: nav
!type(eph_t), pointer :: nav_eph(:)
!type(eph_t) eph0
!integer*4 i,j
!if (nav%n<=0) return
!call init_eph(eph0)
!call qsort(nav%eph,nav%n,sizeof(eph0),cmpeph)
!i=2; j=1
!do i=2,nav%n
!    if (nav%eph(i)%sat /=nav%eph(j)%sat .or. nav%eph(i)%iode/=nav%eph(j)%iode)then
!        j=j+1; nav%eph(j)=nav%eph(i)
!    endif
!enddo
!nav%n=j
!allocate(nav_eph(nav%n))
!nav_eph(1:nav%n)=nav%eph(1:nav%n)
!if(nav%n/=0) deallocate(nav%eph)
!nav%eph=>nav_eph
!nullify(nav_eph)
!nav%nmax=nav%n
!end subroutine
!
!! sort and unique glonass ephemeris -----------------------------------------
!subroutine uniqgeph(nav)
!implicit none
!type(nav_t), intent(out) :: nav
!type(geph_t), pointer :: nav_geph(:)
!type(geph_t) geph0
!integer*4 i,j
!if (nav%ng<=0) return
!call init_geph(geph0)
!call qsort(nav%geph,nav%ng,sizeof(geph0),cmpgeph)
!i=1; j=1
!do i=1,nav%ng
!    if (nav%geph(i)%sat     /=nav%geph(j)%sat .or. &
!        nav%geph(i)%toe%time/=nav%geph(j)%toe%time .or. &
!        nav%geph(i)%svh     /=nav%geph(j)%svh)then
!        j=j+1; nav%geph(j)=nav%geph(i)
!    endif
!enddo
!nav%ng=j
!allocate(nav_geph(nav%ng))
!nav_geph(1:nav%ng)=nav%geph(1:nav%ng)
!if(nav%ng/=0) deallocate(nav%geph)
!nav%geph=>nav_geph
!nullify(nav_geph)
!nav%ngmax=nav%ng
!end subroutine

! unique ephemerides ----------------------------------------------------------
! unique ephemerides in navigation data and update carrier wave length
! args   : nav_t *nav    IO     navigation data
! return : number of epochs
!-----------------------------------------------------------------------------
subroutine uniqnav(nav)
implicit none
type(nav_t), intent(out) :: nav
integer*4 i,j

! unique ephemeris 
!call uniqeph (nav)
!call uniqgeph(nav)

! update carrier wave length 
do i=1,MAXSAT
    do j=1,NFREQ
        nav%lam(i,j)=satwavelen(i,j-1,nav)
    enddo
enddo
end subroutine

! screen by time --------------------------------------------------------------
! screening by time start, time end, and time interval
! args   : gtime_t time  I      time
!          gtime_t ts    I      time start (ts.time==0:no screening by ts)
!          gtime_t te    I      time end   (te.time==0:no screening by te)
!          double  tint  I      time interval (s) (0.0:no screen by tint)
! return : 1:on condition, 0:not on condition
! -----------------------------------------------------------------------------
integer*4 function screent(time, ts, te, tint)
implicit none
type(gtime_t), intent(in) :: time, ts, te
real*8, intent(in) :: tint
integer*4 week
real*8 gpst, dwnd
logical(1) l1, l2, l3
call time2gpst(time,week,gpst)
dwnd=min(tint/10.d0,0.3d0)
! l1=(tint<=0.d0 .or. dmod(gpst+DTTOL,tint)<=DTTOL*2.d0)
! l1=(tint<=0.d0 .or. dmod(timediff(time,ts),tint)>=tint-DTTOL .or. dmod(timediff(time,ts),tint)<=DTTOL)
l1=(tint<=0.d0 .or. dmod(timediff(time,ts),tint)>=tint-dwnd .or. dmod(timediff(time,ts),tint)<=dwnd)
l2=(ts%time==0 .or. timediff(time,ts)>=-dwnd)
l3=(te%time==0 .or. timediff(time,te)<  dwnd)
if(l1 .and. l2 .and. l3)then
    screent=1
else
    screent=0
endif
end function

! compute dops ----------------------------------------------------------------
! * * compute DOP (dilution of precision)
! * * args   : integer*4    ns        I   number of satellites
! * *          real*8 *azel     I   satellite azimuth/elevation angle (rad)
! * *          real*8 elmin     I   elevation cutoff angle (rad)
! * *          real*8 *dop      O   DOPs {GDOP,PDOP,HDOP,VDOP}
! * * return : none
! * * notes  : dop(0)-(3) return 0 in case of dop computation error
! * *---------------------------------------------------------------------------

!define SQRT(x)     ((x)<0.d0 .or. (x)/=(x)?0.d0:dsqrt(x))

subroutine dops(ns, azel, elmin, dop)
implicit none
integer*4, intent(in) :: ns
real*8, intent(in) :: azel(ns,2),elmin
real*8, intent(out) :: dop(4)
real*8 H(MAXSAT,4),Q(4,4),cosel,sinel
integer*4 :: i,n,info
i=1; n=1
dop=0.d0
do i=1,min(ns,MAXSAT)
    if (azel(i,2)<elmin .or. azel(i,2)<=0.0) cycle
    cosel=dcos(azel(i,2))
    sinel=dsin(azel(i,2))
    H(n,1)=cosel*dsin(azel(i,1))
    H(n,2)=cosel*dcos(azel(i,1))
    H(n,3)=sinel
    H(n,4)=1.d0
    n=n+1
enddo
if (n-1<4) return

call matmul2("TN",4,4,n-1,H(1:n-1,:),H(1:n-1,:),Q)
call matinv(Q,4,info)
if (info==0)then
    dop(1)=SQRT2(Q(1,1)+Q(2,2)+Q(3,3)+Q(4,4))  ! GDOP 
    dop(2)=SQRT2(Q(1,1)+Q(2,2)+Q(3,3))         ! PDOP 
    dop(3)=SQRT2(Q(1,1)+Q(2,2))                ! HDOP 
    dop(4)=SQRT2(Q(3,3))                       ! VDOP 
endif
end subroutine

! satellite carrier wave length -----------------------------------------------
! get satellite carrier wave lengths
! args   : integer*4 sat    I   satellite number
!          integer*4 frq    I   frequency index (0:L1,1:L2,2:L5/3,...)
!          nav_t  *nav      I   navigation messages
! return : carrier wave length (m) (0.d0: error)
!-----------------------------------------------------------------------------
real*8 function satwavelen(sat, frq, nav)
implicit none
integer*4, intent(in) :: sat, frq
type(nav_t), intent(in) :: nav
real*8 freq_glo(2), dfrq_glo(2)
integer*4 i, sys, prn
freq_glo=(/FREQ1_GLO,FREQ2_GLO/)
dfrq_glo=(/DFRQ1_GLO,DFRQ2_GLO/)
call satsys(sat, prn, sys)
if (sys==SYS_GLO)then
    if (0<=frq .and. frq<=1)then  ! L1,L2 
        do i=1,nav%ng
            if (nav%geph(i)%sat/=sat) cycle
            satwavelen=CLIGHT/(freq_glo(frq+1)+dfrq_glo(frq+1)*nav%geph(i)%frq)
            return
        enddo
    elseif (frq==2)then  ! L3 
        satwavelen=CLIGHT/FREQ3_GLO; return
    endif
elseif (sys==SYS_CMP)then
    if (frq==0)then
        satwavelen=CLIGHT/FREQ1_CMP; return  ! B1 
    endif
    if (frq==1)then
        satwavelen=CLIGHT/FREQ2_CMP; return  ! B2 
    endif
    if (frq==2)then
        satwavelen=CLIGHT/FREQ3_CMP; return  ! B3 
    endif
elseif (sys==SYS_GAL)then
    if (frq==0)then
        satwavelen=CLIGHT/FREQ1; return  ! E1 
    endif
    if (frq==1)then
        satwavelen=CLIGHT/FREQ7; return  ! E5b 
    endif
    if (frq==2)then
        satwavelen=CLIGHT/FREQ5; return  ! E5a 
    endif
    if (frq==3)then
        satwavelen=CLIGHT/FREQ6; return  ! E6 
    endif
    if (frq==5)then
        satwavelen=CLIGHT/FREQ8; return  ! E5ab 
    endif
else  ! GPS,QZS 
    if (frq==0)then
        satwavelen=CLIGHT/FREQ1; return  ! L1 
    endif
    if (frq==1)then
        satwavelen=CLIGHT/FREQ2; return  ! L2 
    endif
    if (frq==2)then
        satwavelen=CLIGHT/FREQ5; return  ! L5 
    endif
    if (frq==3)then
        satwavelen=CLIGHT/FREQ6; return  ! L6/LEX 
    endif
    if (frq==6)then
        satwavelen=CLIGHT/FREQ9; return  ! S 
    endif
endif
satwavelen=0.d0
end function

! geometric distance ---------------------------------------------------------
! compute geometric distance and receiver-to-satellite unit vector
! args   : real*8 *rs       I   satellilte position (ecef at transmission) (m)
!          real*8 *rr       I   receiver position (ecef at reception) (m)
!          real*8 *e        O   line-of-sight vector (ecef)
! return : geometric distance (m) (0>:error/no satellite position)
! notes  : distance includes sagnac effect correction
!-----------------------------------------------------------------------------
subroutine geodist(rs, rr, e, dist)
implicit none
real*8, intent(in) :: rs(3),rr(3)
real*8, intent(out) :: e(3),dist
real*8 r
integer*4 i
if(norm(rs,3)<RE_WGS84)then
    dist=-1.d0; return
endif
e=rs-rr
r=norm(e,3)
e=e/r
dist=r+OMGE*(rs(1)*rr(2)-rs(2)*rr(1))/CLIGHT
end subroutine

! satellite azimuth/elevation angle -------------------------------------------
! compute satellite azimuth/elevation angle
! args   : real*8 *pos      I   geodetic position {lat,lon,h} (rad,m)
!          real*8 *e        I   receiver-to-satellilte unit vevtor (ecef)
!          real*8 *azel     IO  azimuth/elevation {az,el} (rad) (NULL: no output)
!                               (0.d0<=azel(0)<2*pi,-pi/2<=azel(1)<=pi/2)
! return : elevation angle (rad)
!-----------------------------------------------------------------------------
subroutine satazel(pos, e, azel, el)
implicit none
real*8, intent(in) :: pos(3),e(3,1)
real*8, intent(out) :: azel(2),el
real*8 :: az, enu(3)
az=0; el=pi/2.d0
if(pos(3)>-RE_WGS84)then
    call ecef2enu(pos,e,enu)
    az=rcond(dot(enu,enu,2)<1d-12,0.d0,datan2(enu(1),enu(2)))
    if (az<0.d0) az=az+2.d0*PI
    el=dasin(enu(3))
endif
azel(1)=az; azel(2)=el
end subroutine

! ionosphere model ------------------------------------------------------------
! compute ionospheric delay by broadcast ionosphere model (klobuchar model)
! args   : gtime_t t        I   time (gpst)
!          real*8 *ion      I   iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3}
!          real*8 *pos      I   receiver position {lat,lon,h} (rad,m)
!          real*8 *azel     I   azimuth/elevation angle {az,el} (rad)
! return : ionospheric delay (L1) (m)
!-----------------------------------------------------------------------------
real*8 function ionmodel(t, ion_in, pos, azel)
implicit none
type(gtime_t), intent(in) :: t
real*8, intent(in) :: ion_in(8),pos(3),azel(2)
real*8 :: ion_default(8)=(/&
    0.1118d-07,-0.7451d-08,-0.5961d-07, 0.1192d-06,&
    0.1167d+06,-0.2294d+06,-0.1311d+06, 0.1049d+07 /)
real*8 tt,f,psi,phi,lam,amp,per,x,ion(8),gpst
integer*4 week
if (pos(3)<-1d3 .or. azel(2)<=0.d0)then
    ionmodel=0.d0; return
endif
ion=ion_in
if (norm(ion,8)<=0.d0) ion=ion_default

! earth centered angle (semi-circle) 
psi=0.0137d0/(azel(2)/PI+0.11d0)-0.022d0

! subionospheric latitude/longitude (semi-circle) 
phi=pos(1)/PI+psi*dcos(azel(1))
if (phi> 0.416d0) phi= 0.416d0
if (phi<-0.416d0) phi=-0.416d0
lam=pos(2)/PI+psi*dsin(azel(1))/dcos(phi*PI)

! geomagnetic latitude (semi-circle) 
phi=phi+0.064d0*dcos((lam-1.617d0)*PI)

! local time (s) 
call time2gpst(t,week,gpst)
tt=43200.d0*lam+gpst
tt=tt-floor(tt/86400.d0)*86400.d0; ! 0<=tt<86400 

! slant factor 
f=1.d0+16.d0*(0.53-azel(2)/PI)**3.d0

! ionospheric delay 
amp=ion(1)+phi*(ion(2)+phi*(ion(3)+phi*ion(4)))
per=ion(5)+phi*(ion(6)+phi*(ion(7)+phi*ion(8)))
amp=rcond(amp<    0.d0,    0.d0,amp)
per=rcond(per<72000.d0,72000.d0,per)
x=2.d0*PI*(tt-50400.d0)/per

ionmodel=CLIGHT*f*rcond(dabs(x)<1.57d0,5d-9+amp*(1.d0+x*x*(-0.5+x*x/24.d0)),5d-9)
end function

! troposphere model -----------------------------------------------------------
! compute tropospheric delay by standard atmosphere and saastamoinen model
! args   : gtime_t time     I   time
!          real*8 *pos      I   receiver position {lat,lon,h} (rad,m)
!          real*8 *azel     I   azimuth/elevation angle {az,el} (rad)
!          real*8 humi      I   relative humidity
! return : tropospheric delay (m)
!-----------------------------------------------------------------------------
real*8 function tropmodel(time, pos, azel, humi)
implicit none
type(gtime_t), intent(in) :: time
real*8, intent(in) :: pos(3),azel(2),humi
real*8 :: hgt,pres,temp,e,z,trph,trpw,temp0
temp0=15.d0
if (pos(3)<-100.d0 .or. 1d4<pos(3) .or. azel(2)<=0)then
    tropmodel=0.d0; return
endif

! standard atmosphere 
hgt=rcond(pos(3)<0.d0,0.d0,pos(3))
pres=1013.25d0*(1.d0-2.2557d-5*hgt)**5.2568d0
temp=temp0-6.5d-3*hgt+273.16d0
e=6.108d0*humi*dexp((17.15d0*temp-4684.d0)/(temp-38.45d0))

! saastamoninen model 
z=PI/2.d0-azel(2)
trph=0.0022768d0*pres/(1.d0-0.00266d0*dcos(2.d0*pos(1))-0.00028d0*hgt/1d3)/dcos(z)
trpw=0.002277d0*(1255.d0/temp+0.05d0)*e/dcos(z)
tropmodel=trph+trpw
end function

! date function -----------------------------------------------------------
! ydoy2mjd
integer*4 function ydoy2mjd_s(iyear, idoy)
implicit none
integer*4, intent(in) :: iyear, idoy
integer*4 imonth, iday
call ydoy2ymd_s(iyear, idoy, imonth, iday)
ydoy2mjd_s=ymd2mjd_s(iyear, imonth, iday)
end function

! mjd ==> ydoy
subroutine mjd2ydoy_s(jd, iyear, idoy)
implicit none
integer*4 jd, iyear, idoy
iyear = (jd + 678940)/365
idoy = jd - ymd2mjd_s(iyear, 1, 1)
do while (idoy .le. 0)
  iyear = iyear - 1
  idoy = jd - ymd2mjd_s(iyear, 1, 1) + 1
enddo
return
end

! ydoy ==> ymd
subroutine ydoy2ymd_s(iyear, idoy, imonth, iday)
implicit none
integer*4 iyear, idoy, imonth, iday
!! local
integer*4 days_in_month(12), id
data days_in_month/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
days_in_month(2) = 28
if (mod(iyear, 4) .eq. 0 .and. (mod(iyear, 100) .ne. 0 .or. mod(iyear, 400) .eq. 0)) days_in_month(2) = 29
id = idoy
do imonth = 1, 12
  id = id - days_in_month(imonth)
  if (id .gt. 0) cycle
  iday = id + days_in_month(imonth)
  exit
enddo
return
end

! ymd ==> mjd
integer*4 function ymd2mjd_s(iyear, imonth, iday)
implicit none
integer*4 imonth, iday, iyear
!! local
integer*4 iyr, doy_of_month(12), modified_julday
data doy_of_month/0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334/
!! check the input data
if ((iyear .lt. 0 .or. imonth .lt. 0 .or. iday .lt. 0 .or. imonth .gt. 12 .or. iday .gt. 366) &
    .or. (imonth .ne. 0 .and. iday .gt. 31)) then
  write (*, '(a)') '***ERROR(modified_julday, ymd2mjd_s): incorrect date(year,month,day): ', &
    iyear, imonth, iday
  call exit(1)
endif
iyr = iyear
if (imonth .le. 2) iyr = iyr - 1
modified_julday = 365*iyear - 678941 + iyr/4 - iyr/100 + iyr/400 + iday
if (imonth .ne. 0) modified_julday = modified_julday + doy_of_month(imonth)
ymd2mjd_s=modified_julday
return
end
end module rtkcmn_f90_
