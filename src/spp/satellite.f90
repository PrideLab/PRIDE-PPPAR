
! Satellite relevant function -----------------------------------------------
! satellite number to satellite system --------------------------------------
! convert satellite number to satellite system
! args   : integer*4    sat       I   satellite number (1-MAXSAT)
!          integer*4    prn       O   satellite prn/slot number (NULL: no output)
!          integer*4    sys       O
! return : satellite system (SYS_GPS,SYS_GLO,...)
! ---------------------------------------------------------------------------
subroutine satsys(sat, prn, sys)
implicit none
include 'file_para.h'
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

! satellite carrier wave length ---------------------------------------------
! get satellite carrier wave lengths
! args   : integer*4 sat    I   satellite number
!          integer*4 frq    I   frequency index (0:L1,1:L2,2:L5/3,...)
!          nav_t  *nav      I   navigation messages
! return : carrier wave length (m) (0.d0: error)
!----------------------------------------------------------------------------
real*8 function satwavelen(sat, frq, nav)
implicit none
include 'file_para.h'
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

! satellite id to satellite number ------------------------------------------
!* * convert satellite id to satellite number
!* * args   : char   *id       I   satellite id (nn,Gnn,Rnn,Enn,Jnn,Cnn,Inn or Snn)
!* * return : satellite number (0: error)
!* * notes  : 120-142 and 193-199 are also recognized as sbas and qzss
!* *-------------------------------------------------------------------------
integer*4 function satid2no(id)
implicit none
include 'file_para.h'
character(*),intent(in) :: id
integer*4 sys,prn
character(1) code
integer*4 err1, err2, satno
external :: satno
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

! satellite system+prn/slot number to satellite number ----------------------
! convert satellite system+prn/slot number to satellite number
! args   : integer*4    sys       I   satellite system (SYS_GPS,SYS_GLO,...)
!          integer*4    prn       I   satellite prn/slot number
! return : satellite number (0:error)
! ---------------------------------------------------------------------------
integer*4 function satno(sys, prn)
implicit none
include 'file_para.h'
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

! satellite number to satellite id ------------------------------------------
! convert satellite number to satellite id
! args   : integer*4    sat       I   satellite number
!          char   *id             O   satellite id (Gnn,Rnn,Enn,Jnn,Cnn,Inn or nnn)
! return : none
! ---------------------------------------------------------------------------
subroutine satno2id(sat, id)
implicit none
include 'file_para.h'
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

! test excluded satellite ---------------------------------------------------
! test excluded satellite
! args   : integer*4    sat       I   satellite number
!          integer*4    svh       I   sv health flag
!          prcopt_t *opt          I   processing options (NULL: not used)
! return : status (1:excluded,0:not excluded)
!----------------------------------------------------------------------------
integer*4 function satexclude(sat, var, svh, opt, isopt)
implicit none
include 'file_para.h'
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
