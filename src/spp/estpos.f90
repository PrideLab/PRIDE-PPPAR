
! Estimate receiver position ------------------------------------------------
subroutine estpos(obs, n, rs, dts, vare, svh, nav, opt, sol, azel, vsat, resp, msg, stat)
implicit none
include 'file_para.h'
integer*4, intent(in) :: n, svh(n)
type(obsd_t), intent(in) :: obs(n)
real*8, intent(in) :: rs(n,6), dts(n,2), vare(n)
type(nav_t), intent(in) :: nav
type(prcopt_t), intent(in) :: opt
type(sol_t), intent(out) :: sol
real*8, intent(out) :: azel(n,2), resp(n)
integer*4, intent(out) :: vsat(n), stat
character(*), intent(out) :: msg
real*8 :: x(NX),dx(NX),Q(NX,NX),v(n+4,1),H(n+4,NX),var(n+4,1),sig,HT(NX,n+4)
integer*4 i,j,k,info,stat1,nv,ns
real*8 thres, fact, vara(n), rr0(3), maxresp  ! valpos
integer*4 valres(MAXSAT), nslim, nt, nf  ! number of true / false
integer*4 navsy2,ii,prn,sy2,method,re1
integer*4 valrag(MAXSAT), maxindex
real*8, external :: norm, rcond
integer*4, external :: icond
type(gtime_t), external :: timeadd
method=1; re1=0  ! 1: 1d-4, 2: nf<nt
rr0=sol%rr(1:3)
110 valres=0
120 thres=3d0; fact=thres*thres
x=0.d0; nslim=0; nt=0; nf=0; valrag=0
x(1:3)=rr0  ! 0.d0
do i=1,MAXITR
    ! pseudorange residuals 
100 call rescode(i-1,obs,n,rs,dts,vare,svh,nav,x,opt,v,H,var,vara,valres,valrag,azel,vsat,resp,ns,nv)
    
    nslim=3; navsy2=SYS_NONE
    do ii=1,n
        if(valrag(ii)==0)then
            call satsys(obs(ii)%sat,prn,sy2)
            navsy2=or(navsy2,sy2)
        endif
    enddo
    
    if(and(navsy2,SYS_GPS)/=0) nslim=nslim+1
    if(and(navsy2,SYS_GLO)/=0) nslim=nslim+1
    if(and(navsy2,SYS_GAL)/=0) nslim=nslim+1
    if(and(navsy2,SYS_CMP)/=0) nslim=nslim+1
    
    if (ns<nslim)then
        write(msg,"('lack of valid sats ns=',I4)") ns; exit
    endif
    ! weight by variance 
    do j=1,nv
        sig=dsqrt(var(j,1))
        v(j,1)=v(j,1)/sig
        do k=1,NX
            H(j,k)=H(j,k)/sig
        enddo
    enddo
    ! least square estimation 
    call mattran(n+4,NX,H,HT)
    call lsq(HT,v,NX,nv,dx,Q,info)
    if(info/=0)then
        write(msg,"('lsq error info=',I4)") info; exit
    endif
    do j=1,NX
        x(j)=x(j)+dx(j)
    enddo
    if(i==5) then ! del outlier
        maxresp=0.d0
        maxindex=0
        do j=1,n 
           if(dabs(resp(j))>=maxresp) then
                maxresp=dabs(resp(j))
                maxindex=j
            endif
        enddo
        if(maxindex/=0 .and. maxresp>=10) valres(obs(maxindex)%sat)=-1
    endif
    if(method==2)then  ! nf<nt
        nt=0; nf=0
        do j=1,n
            if(dabs(resp(j))>=1d-20 .and. resp(j)**2> fact*vara(j)) nf=nf+1
            if(dabs(resp(j))>=1d-20 .and. resp(j)**2<=fact*vara(j)) nt=nt+1
        enddo
        if(0<nf .and. nf<nt)then  ! .and. nf<=4
            do j=1,n
                if(dabs(resp(j))>=1d-20 .and. resp(j)**2>fact*vara(j))then
                    valres(obs(j)%sat)=-1
                endif
            enddo
            goto 100
        endif
    endif
    if(norm(dx,NX)<1d-4)then
        if(method==1)then  ! 1d-4
            nt=0; nf=0
            do j=1,n
                if(dabs(resp(j))>=1d-20 .and. resp(j)**2> fact*vara(j)) nf=nf+1
                if(dabs(resp(j))>=1d-20 .and. resp(j)**2<=fact*vara(j)) nt=nt+1
            enddo
            if(0<nf .and. nf<nt)then  ! .and. nf<=4
                do j=1,n
                    if(dabs(resp(j))>=1d-20 .and. resp(j)**2>fact*vara(j))then
                        valres(obs(j)%sat)=-1
                    endif
                enddo
                goto 100
            endif
        endif
        sol%mytype=0
        sol%time=timeadd(obs(1)%time,-x(4)/CLIGHT)
        sol%dtr(1)=x(4)/CLIGHT  ! receiver clock bias (s) 
        sol%dtr(2)=x(5)/CLIGHT  ! glo-gps time offset (s) 
        sol%dtr(3)=x(6)/CLIGHT  ! gal-gps time offset (s) 
        sol%dtr(4)=x(7)/CLIGHT  ! bds-gps time offset (s) 
        do j=1,6
            sol%rr(j)=rcond(j<=3,x(j),0.d0)
        enddo
        do j=1,3
            sol%qr(j)=Q(j,j)
        enddo
        sol%qr(4)=Q(1,2)     ! cov xy 
        sol%qr(5)=Q(2,3)     ! cov yz 
        sol%qr(6)=Q(1,3)     ! cov zx 
        sol%ns=ns
        ! validate solution 
        call valsol(azel,vsat,n,opt,v,nv,NX,sol,msg,stat1)
        sol%stat=icond(opt%sateph==EPHOPT_SBAS,SOLQ_SBAS,SOLQ_SINGLE)
        stat=stat1
        if(stat1==1) return
    endif
enddo
if(method==1)then
    if(re1<=1)then
        nf=0
        do j=1,n
            if(dabs(resp(j))>=1d-20 .and. resp(j)**2>fact*vara(j))then
                valres(obs(j)%sat)=-1; nf=nf+1
            endif
        enddo
        re1=re1+1
        if(nf>0) goto 120
    else
        method=2; goto 110
    endif
    method=2; goto 110
endif
if (i>MAXITR) write(msg,"('iteration divergent i=',I4)") i-1
stat=0
end subroutine
