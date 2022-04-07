
! Pseudorange residuals -----------------------------------------------------
subroutine rescode(iter, obs, n, rs, dts, vare, svh, nav, x, opt, v, H, var, vara, valres, valrag, azel, vsat, resp, ns, stat)
implicit none
include 'file_para.h'
integer*4, intent(in) :: iter, n, svh(n), valres(MAXSAT)
type(obsd_t), intent(in) :: obs(n)
real*8, intent(in) :: rs(n,6), dts(n,2), vare(n), x(NX)
type(nav_t), intent(in) :: nav
type(prcopt_t), intent(in) :: opt
real*8, intent(out) :: v(n+4), H(n+4,NX), var(n+4), vara(n), azel(n,2), resp(n)
integer*4, intent(out) :: valrag(MAXSAT), vsat(n), ns, stat
real*8 r,dion,dtrp,vmeas,vion,vtrp,rr(3),pos(3),dtr,e(3),P,lam_L1,el,SQR,rcond,varerr
integer*4 :: i,j,nv,sys,mask(4),prn,info,satexclude,icond
external :: satexclude,icond,SQR,rcond,varerr
logical*1 isopt
nv=1; mask=0
rr(1:3)=x(1:3)
dtr=x(4)
valrag=0
call ecef2pos(rr,pos)
ns=0; i=1
do while(i<=min(n,MAXOBS))
    vsat(i)=0; azel(i,:)=0.d0; resp(i)=0.d0; vara(i)=0.d0
    call satsys(obs(i)%sat,prn,sys)
    if (sys==0)then
        i=i+1; cycle
    endif
    
    ! mark the gross error (valpos)
    if(valres(obs(i)%sat)==-1)then
        i=i+1; cycle
    endif
    
    ! reject duplicated observation data 
    if(i<n .and. i<MAXOBS)then
        if(obs(i)%sat==obs(i+1)%sat)then
            i=i+2; cycle
        endif
    endif
    
    ! geometric distance/azimuth/elevation angle 
    call geodist(rs(i,:),rr,e,r)
    call satazel(pos,e,azel(i,:),el)
    if (r<=0.d0 .or. el<opt%elmin)then
        i=i+1; cycle
    endif
    
    ! psudorange with code bias correction 
    call prange(obs(i),nav,azel(i,:),iter,opt,vmeas,P)
    if (dabs(P)<=1d-20)then
        valrag(i)=-1; i=i+1; cycle
    endif
    
    ! excluded satellite? 
    isopt=.true.
    if (satexclude(obs(i)%sat,vare(i),svh(i),opt,isopt)/=0)then
        i=i+1; cycle
    endif
    
    ! ionospheric corrections 
    call ionocorr(obs(i)%time,nav,obs(i)%sat,pos,azel(i,:),icond(iter>0,opt%ionoopt,IONOOPT_BRDC),dion,vion,info)
    if (info==0)then
        i=i+1; cycle
    endif
    
    ! GPS-L1 -> L1/B1 
    lam_L1=nav%lam(obs(i)%sat,1)
    if (lam_L1>0.d0)then
        dion=dion*SQR(lam_L1/lam_carr(1))
    endif
    
    ! tropospheric corrections 
    call tropcorr(obs(i)%time,nav,pos,azel(i,:),icond(iter>0,opt%tropopt,TROPOPT_SAAS),dtrp,vtrp,info)
    if (info==0)then
        i=i+1; cycle
    endif
    
    ! pseudorange residual 
    v(nv)=P-(r+dtr-CLIGHT*dts(i,1)+dion+dtrp)
    
    ! design matrix 
    do j=1,NX
        !H(nv,j)=rcond(j<4,-e(j),rcond(j==4,1.d0,0.d0))
        if(j<4)then
            H(nv,j)=-e(j)
        else
            H(nv,j)=rcond(j==4,1.d0,0.d0)
        endif
    enddo
    
    ! time system and receiver bias offset correction 
    if(sys==SYS_GLO)then
        v(nv)=v(nv)-x(5); H(nv,5)=1.d0; mask(2)=1
    elseif(sys==SYS_GAL)then
        v(nv)=v(nv)-x(6); H(nv,6)=1.d0; mask(3)=1
    elseif(sys==SYS_CMP)then
        v(nv)=v(nv)-x(7); H(nv,7)=1.d0; mask(4)=1
    else
        mask(1)=1
    endif
    vsat(i)=1; resp(i)=v(nv); ns=ns+1
    
    ! error variance 
    var(nv)=varerr(opt,azel(i,2),sys)+vare(i)+vmeas+vion+vtrp  ! 
    vara(i)=var(nv) ! all
    nv=nv+1; i=i+1
enddo

! constraint to avoid rank-deficient 
do i=1,4
    if(mask(i)/=0) cycle
    v(nv)=0.d0
    do j=1,NX
        H(nv,j)=rcond(j==i+3,1.d0,0.d0)
    enddo
    var(nv)=0.01d0; nv=nv+1
enddo
stat=nv-1
end subroutine
