
! Decode observation data ---------------------------------------------------
! decode obs header ---------------------------------------------------------
subroutine decode_obsh(fp, buff, ver, tsys, tobs, tint, sta)
implicit none
include 'file_para.h'
integer*4, intent(in) :: fp
character(*), intent(inout) :: buff
real*8, intent(in) :: ver
integer*4, intent(out) :: tsys
character(3), intent(out) :: tobs(NUMSYS,MAXOBSTYPE)  !tobs(:,:)
real*8, intent(out) :: tint
type(sta_t), intent(out) :: sta

! default codes for unknown code 
character(8) :: defcodes(7)=(/&
    'CWX    ',&  ! GPS: L125____
    'CC     ',&  ! GLO: L12_____
    'X XXXX ',&  ! GAL: L1_5678_
    'CXXX   ',&  ! QZS: L1256___
    'C X    ',&  ! SBS: L1_5____
    'XIXIIX ',&  ! BDS: L125678_
    '  A   A'/)  ! IRN: L__5___9
real*8 del(3)
integer*4 i,j,k,n,nt,prn,fcn,info,posi
character(1) p
character(20) label
character(4) str
real*8, external :: str2num
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
    if(index(systyps,p)==0) return
    i=index(systyps,p)
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
        if(index(frqtyps,tobs(i,j)(2:2))==0) cycle
        posi=index(frqtyps,tobs(i,j)(2:2))
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
    read(buff(1:11),*,iostat=info) tint
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
    ! nav%leaps=int(str2num(buff,1,6))
elseif(index(label,'# OF SALTELLITES'    )/=0)then  ! opt 
elseif(index(label,'PRN / # OF OBS'      )/=0)then  ! opt 
endif
end subroutine

! decode obs epoch ----------------------------------------------------------
subroutine decode_obsepoch(fp, buff, ver, time, flag, sats, stat)
implicit none
include 'file_para.h'
integer*4, intent(in) :: fp
character(*), intent(inout) :: buff
real*8, intent(in) :: ver
type(gtime_t), intent(out) :: time
integer*4, intent(out) :: flag, sats(MAXOBS), stat
integer*4 i,j,n,info
character(8) :: satid
integer*4, external :: satid2no
satid=''
if(ver<=2.99)then  ! ver.2 
    read(buff(30:32),*,iostat=info) n
    if(n<=0.or.info/=0)then
        stat=0; return
    endif
    ! epoch flag: 3:new site,4:header info,5:external event 
    read(buff(29:29),*,iostat=info) flag
    if(info/=0)then  ! format check
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
    if(info/=0)then  ! format check
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

! cycle obs data -----------------------------------------------------------
subroutine cycle_obsdata(fp, ver, index_in)
implicit none
include 'file_para.h'
integer*4, intent(in) :: fp
real*8, intent(in) :: ver
type(sigind_t), intent(in) :: index_in(7)
! local
character(MAXRNXLEN) :: buff
type(sigind_t) ind
character(8) :: satid
integer*4 :: i,j,sys,info
integer*4, external :: satid2no,icond
satid=''
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
    i=i+1; j=j+16
enddo
end subroutine

! decode obs data -----------------------------------------------------------
subroutine decode_obsdata(fp, buff, ver, mask, index_in, obs, stat)
implicit none
include 'file_para.h'
integer*4, intent(in) :: fp, mask
character(*), intent(inout) :: buff
real*8, intent(in) :: ver
type(sigind_t), intent(in) :: index_in(7)
type(obsd_t), intent(out) :: obs
integer*4, intent(out) :: stat
type(sigind_t) ind
real*8 :: val(MAXOBSTYPE)
integer*4 :: lli(MAXOBSTYPE), qual(MAXOBSTYPE)
character(8) :: satid
integer*4 :: i,j,n,m,stat1,p(MAXOBSTYPE),k(16),l(16),prn,sys,info
integer*4, external :: satid2no,icond
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
        read(buff(j:j+13),*,iostat=info) val(i); val(i)=val(i)+ind%shift(i)
        !read(buff(j:j+13),"(D14.3)",iostat=info) val(i); val(i)=val(i)+ind%shift(i)
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
        if(obs%P(p(i)+1)<=1d-20)     obs%P(p(i)+1)=val(i)
        if(obs%code(p(i)+1)<=1d-20)  obs%code(p(i)+1)=ind%code(i)
        if(obs%qualP(p(i)+1)<=1d-20) obs%qualP(p(i)+1)=icond(qual(i)>0,qual(i),1)
    case(1)
        if(obs%L(p(i)+1)<=1d-20)     obs%L(p(i)+1)=val(i)
        if(obs%LLI(p(i)+1)<=1d-20)   obs%LLI(p(i)+1)=lli(i)
        if(obs%qualL(p(i)+1)<=1d-20) obs%qualL(p(i)+1)=icond(qual(i)>0,qual(i),1)
    case(2)
        if(obs%D(p(i)+1)<=1d-20)     obs%D(p(i)+1)=real(val(i))
    case(3)
        if(obs%SNR(p(i)+1)<=1d-20)   obs%SNR(p(i)+1)=int(val(i)*4.d0+0.5d0)
    end select
enddo
stat=1
end subroutine

! decode marker name from filepath ------------------------------------------
subroutine decodemarker(fpath, marker)
implicit none
include 'file_para.h'
character(*), intent(in) :: fpath
character(*), intent(out) :: marker
integer*4 flen
character(128) fname
marker=''
flen=len_trim(fpath)
if(flen==0) return
if(fpath(flen:flen)=='o'.or.fpath(flen:flen)=='O')then
    call getfname(fpath,fname)
    call setstr(marker,fname,min(20,len_trim(fname))); return
endif
end subroutine
