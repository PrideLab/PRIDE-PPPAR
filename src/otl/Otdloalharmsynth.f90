!******************************************************************************************
!     Copyright (C) 2024 ZHANG Chuanyin
!
!     This module was written by ZHANG Chuanyin and is used in the PRIDE PPP-AR project.
!     All rights reserved.
!
!     This module is protected by copyright law. Unauthorized copying, modification, publication,
!     transmission, display, performance, sale, distribution of this module's code, in whole or
!     in part, is strictly prohibited without the prior written permission of the original author.
!     You may use this module code in accordance with the license agreement of the original author
!     and PRIDE Lab.
!
! Website: www.zcyphygeodesy.com
! Contact Information: zhangchy@casm.ac.cn
!******************************************************************************************
subroutine Otdloalharmsynth(plon,plat,phgt,mjds,sods,mjde,sode,interval,DIM,lfn_outfil)
  implicit none
  include '../header/const.h'
  ! args
  real*8, intent(in) :: plon,plat,phgt ! position
  real*8, intent(in) :: sods,sode      ! second of day
  real*8, intent(in) :: interval       ! sample time
  integer*4, intent(in) :: mjds,mjde,DIM ! modified julian day
  integer*4, intent(in) :: lfn_outfil  ! output file
  ! common
  real*8 RAD,GRS(6)
  ! local
  character*100 :: line
  character*4   :: fh
  integer*8 :: i,j,n,nn,mjd0,mjd1
  integer*8 :: stat
  real*8 sod0,sod1,jds,enu(3)
  real*8 flv(4000,3),BLH(3),tmp,st(7)
  real*8 rln(3),NFD(5),gr,t
  real*8,allocatable :: fes(:,:)
  real*8,allocatable :: pnm(:),dpt1(:),dpt2(:),cnm(:),snm(:)
  ! function called
  real*8, external ::timdif
  !---------------------------------------------------------------------
  ! initialize common
  GRS(1)= 3.986004415d14; GRS(2)=6378137.d0; GRS(3)=1.0826359d-3
  GRS(4) = 7.292115d-5; GRS(5)=1.d0/298.25641153d0
  RAD=PI/180.d0
  ! initialize variables
  mjd0=mjds;sod0=sods;mjd1=mjde;sod1=sode
  if(timdif(mjd1,sod1,mjd0,sod0) .lt. 0) goto 902
  call timinc(mjd0, sod0, -interval*DIM/2, mjd0, sod0)
  call timinc(mjd1, sod1,  interval*DIM/2, mjd1, sod1)
  allocate(fes((MAXDEG+2)**2*40,7), STAT=stat)
  if(stat .ne. 0) goto 902
  allocate(pnm((MAXDEG+2)**2),dpt1((MAXDEG+2)**2),dpt2((MAXDEG+2)**2), &
    cnm((MAXDEG+2)**2),snm((MAXDEG+2)**2), STAT=stat)
  if(stat .ne. 0) goto 902
  fes=0.d0
  BLH(1)=plat;BLH(2)=plon;BLH(3)=phgt
  !
  !! read FES2004
  open(unit=8,file="FES2004S1.dat", status="old",iostat=stat)
  if(stat .ne. 0) then
    write (*, '(a)') '***ERROR(Otdloalharmsynth): open file FES2004S1.dat '
    goto 902 
  endif
  do i=1,3
    read(8,'(a)') line
  enddo
  i=0
  do
    read(8,*,end=903,iostat=stat) st(1),fh,(st(j),j=2,3),tmp,tmp,tmp,tmp,(st(j),j=4,7)
    if(stat .ne. 0) exit
    if(st(2)<MAXDEG+1)then
      i=i+1;st(1)=st(1)*1.d3;fes(i,1:7)=st(1:7);
    endif
  enddo
903 continue
  close(8)
  !
  !! read Love_load
  nn=i
  flv=0.d0
  open(unit=8,file="Love_load_cm.dat", status="old",iostat=stat)
  if (stat .ne. 0) then
    write (*, '(a)') '***ERROR(Otdloalharmsynth): open file Love_load_cm.dat '
    goto 902 
  endif
  do i=1,6
    read(8,'(a)') line
  enddo
  n=0
  do while(n<3600)
    n=n+1
    read(8,*,end=904,iostat=stat)i,(flv(n,j),j=1,3)
    if(stat .ne. 0) exit
  enddo
904 continue
  close(8)

  !
  !! calculating  
  call BLH_RLAT(GRS,BLH,rln)
  call GNormalfd(BLH,NFD,GRS)
  gr=NFD(2)
  t=dsin(rln(2)*RAD)
  call BelPnmdt(pnm,dpt1,dpt2,MAXDEG,t)
  do while(timdif(mjd1,sod1,mjd0,sod0) .gt. -0.01)
    jds=mjd0+sod0/86400.d0
    call OLoadDFlu(jds,cnm,snm,MAXDEG,fes,nn)
    call LTideFlupnm(rln,MAXDEG,cnm,snm,flv,enu,GRS,pnm,dpt1,dpt2,gr)
    write(lfn_outfil,'(I6,F10.1,3F14.4)') mjd0,sod0,(enu(i),i=1,3)
    call timinc(mjd0, sod0, interval, mjd0, sod0)
  enddo
902 continue
  deallocate(fes,dpt1,dpt2,cnm,snm)
end