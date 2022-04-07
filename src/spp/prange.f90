
! Psendorange with code bias correction -------------------------------------
subroutine prange(obs, nav, azel, iter, opt, var, range)  !жд
implicit none
include 'file_para.h'
type(obsd_t), intent(in) :: obs
type(nav_t), intent(in) :: nav
real*8, intent(in) :: azel(*)
integer*4, intent(in) :: iter
type(prcopt_t), intent(in) :: opt
real*8, intent(out) :: var, range
real*8 lam(NFREQ)
real*8 PC,P1,P2,P1_P2,P1_C1,P2_C2,gamma,SQR,gettgd
integer*4 :: i,j,sys,prn
external :: SQR,gettgd
i=0; j=1
P1_P2=0
P1_C1=0
P2_C2=0
var=0.d0; range=0.d0
lam=nav%lam(obs%sat,:)  ! test values here
call satsys(obs%sat,prn,sys)
if(sys==0) return

! L1-L2 for GPS/GLO/QZS, L1-L5 for GAL/SBS 
if(NFREQ>=3 .and. and(sys,or(SYS_GAL,SYS_SBS))/=0) j=2
if(NFREQ<2 .or. dabs(lam(i+1))<=1d-20 .or. dabs(lam(j+1))<=1d-20) return

gamma=SQR(lam(j+1))/SQR(lam(i+1))  ! f1^2/f2^2 
P1=obs%P(i+1)
P2=obs%P(j+1)
P1_P2=nav%cbias(obs%sat,1)
P1_C1=nav%cbias(obs%sat,2)
P2_C2=nav%cbias(obs%sat,3)

! if no P1-P2 DCB, use TGD instead 
if(dabs(P1_P2)<=1d-20 .and. and(sys,or(SYS_GPS,or(SYS_GAL,SYS_QZS)))/=0)then
    P1_P2=(1.d0-gamma)*gettgd(obs%sat,nav)
endif
if(opt%ionoopt==IONOOPT_IFLC)then            ! dual-frequency 
    if(dabs(P1)<=1d-20 .or. dabs(P2)<=1d-20) return
    if(obs%code(i+1)==CODE_L1C) P1=P1+P1_C1  ! C1%P1 
    if(obs%code(j+1)==CODE_L2C) P2=P2+P2_C2  ! C2%P2 
    PC=(gamma*P1-P2)/(gamma-1.d0)            ! iono-free combination 
else  ! single-frequency 
    if (dabs(P1)<=1d-20) return
    if (obs%code(i+1)==CODE_L1C) P1=P1+P1_C1 ! C1%P1 
    PC=P1-P1_P2/(1.d0-gamma)
endif
!if (opt%sateph==EPHOPT_SBAS) PC=PC-P1_C1  ! sbas clock based C1 
var=SQR(ERR_CBIAS)
range=PC
end subroutine

! get tgd parameter (m) -----------------------------------------------------
real*8 function gettgd(sat, nav)
implicit none
include 'file_para.h'
integer*4, intent(in) :: sat
type(nav_t), intent(in) :: nav
integer*4 i
do i=1,nav%n
    if(nav%eph(i)%sat/=sat) cycle
    gettgd=CLIGHT*nav%eph(i)%tgd(1); return
enddo
gettgd=0.d0
end function
