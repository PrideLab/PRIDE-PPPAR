
! Validate solution ---------------------------------------------------------
subroutine valsol(azel, vsat, n, opt, v, nv, nx1, sol, msg, stat)
implicit none
include 'file_para.h'
integer*4, intent(in) :: n, vsat(n), nv, nx1
real*8, intent(in) :: azel(n,2), v(nv)
type(prcopt_t), intent(in) :: opt
type(sol_t), intent(out) :: sol
character(*), intent(out) :: msg
integer*4, intent(out) :: stat
real*8 azels(MAXOBS,2),dop(4),vv,dot
integer*4 :: i,ns
external :: dot
i=1; ns=1
! chi-square validation of residuals 
vv=dot(v,v,nv)
if(nv>nx1)then
    if(vv>chisqr(nv-nx1))then
        write(msg,"('chi-square error nv=',I4,' vv=',F14.4,' cs=',F14.4)") nv,vv,chisqr(nv-nx1)
        stat=0; return
        !stat=1; return
    endif
endif
! validation of sigma
if(.not. opt%lsa .and. (sol%qr(1)>=400 .or. sol%qr(2)>=400))then
    stat=0; return
endif
! large gdop check 
do i=1,n
    if (vsat(i)==0) cycle
    azels(ns,1)=azel(i,1)
    azels(ns,2)=azel(i,2)
    ns=ns+1
enddo
call dops(ns-1,azels(1:ns-1,:),opt%elmin,dop)
sol%dop=dop
if(dop(1)<=0.d0 .or. dop(1)>opt%maxgdop)then
    write(msg,"('gdop error nv=',I4' gdop=',F10.4)") nv,dop(1)
    stat=0; return
endif
stat=1
end subroutine

!integer*4 function valres()
subroutine valpos(v, R, nv, thres, vflg)
implicit none
include 'file_para.h'
integer*4, intent(in) :: nv
real*8, intent(in) :: v(nv), R(nv,nv), thres
integer*4, intent(out) :: vflg(nv)
real*8 fact
integer*4 i
fact=thres*thres
do i=1, nv
    if(v(i)*v(i)<=fact*R(i,i))then
        vflg(i)=0
    else
        vflg(i)=1
    endif
enddo
end subroutine
