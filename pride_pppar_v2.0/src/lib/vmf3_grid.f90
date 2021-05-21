!!
!! vmf3_fit.f90
!!
!!    Copyright (C) 2021 by Wuhan University
!!
!!    This program is an open source software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License (version 3) as
!!    published by the Free Software Foundation.
!!
!!    This program is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License (version 3) for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with this program.  If not, see <https://www.gnu.org/licenses/>.
!!
!! Contributor: Jianghui Geng, Songfeng Yang
!! 
!!
!!
!! purpose  : configure VMF3 for static station. read new vmf format 
!!            by J. Boehm at http://vmf.geo.tuwien.ac.at/codes/
!! parameter:
!!    input : flnvmf -- vmf3 file
!!    output: SITE -- station struct
!!
!! note     : vmf3 files should be concatenated to one file
!
subroutine vmf3_grid(flnvmf,SITE)
implicit none
include '../header/const.h'
include '../header/station.h'

character*20 flnvmf
type(station) SITE
!
!! local
logical*1 lfirst,lfill(4),ltim
integer*4 i,lfnh,lfnv,ilat,ilon,iy,imon,id,ih,im,ierr,j
real*8 lat,lon,med,blh(3,4),xyh(3,4),oc(4,4),temp(4),coef(3,4),xy(2),pres,sec
character*256 line
real*8 geod_tmp(3),ellh(360),ylen,y1,y2,blh_tmp(2),blh_tmp0(2),xyh_tmp(2),xy_tmp(2)
!
!! function used
integer*4 get_valid_unit,modified_julday
real*8 ffun

data lfirst/.true./
save lfirst,lfnh,lfnv

if(lfirst) then
  lfirst=.false.
!! grid vmf3 file
  lfnh=get_valid_unit(10)
  open(lfnh,file=flnvmf,status='old',iostat=ierr)
  if(ierr.ne.0) then
    write(*,'(2a)') '***ERROR(vmf3_grid): open file ',flnvmf(1:len_trim(flnvmf))
    call exit(1)
  endif
!! grid ellipsoidal height file
  lfnv=get_valid_unit(10)
  open(lfnv,file='orography_ell_1x1',status='old',iostat=ierr)
  if(ierr.ne.0) then
    write(*,'(a)') '***ERROR(vmf3_grid): open file orography_ell_1x1 '
    call exit(1)
  endif
endif
rewind lfnh
rewind lfnv
!
!! determine which 4 grid points should be used
lat=SITE%geod(1)*180.d0/PI            ! to degrees
lon=SITE%geod(2)*180.d0/PI
if(lat<=89.5d0 .and. lat>-89.5d0) then
  ilat=int((89.5d0-lat)/1.d0)+1         ! latitude  grid interval 1 deg
  ilon=int((lon-0.5d0)/1.d0)+1          ! longitude grid interval 1 deg
  blh(1,1:2)=89.5d0-(ilat-1)*1.d0
  blh(1,3:4)=89.5d0-ilat*1.d0
  blh(2,1)  =(ilon-1)*1.d0+0.5d0
  if(blh(2,1)>=360.d0) blh(2,1)=blh(2,1)-360.d0
  blh(2,3)  =blh(2,1)
  blh(2,2)  =ilon*1.d0+0.5d0
  if(blh(2,2)>=360.d0) blh(2,2)=blh(2,2)-360.d0
  blh(2,4)  =blh(2,2)
  !
  !! read ellipsoidal heights
  do i=1,ilat
    do j=1,360
      read(lfnv,*,err=100) ellh(j)
    enddo
  enddo
  blh(3,1)=ellh(ilon)
  blh(3,2)=ellh(ilon+1)
  do j=1,360
    read(lfnv,*,err=100) ellh(j)
  enddo
  blh(3,3)=ellh(ilon)
  blh(3,4)=ellh(ilon+1)
  !
  !! transform bl to xy using mercator projection
  med=(blh(2,1)+blh(2,2))/2.d0*PI/180.d0
  do i=1,4
    xyh(1,i)=blh(1,i)*PI/180.d0         ! to radian
    xyh(2,i)=blh(2,i)*PI/180.d0
    xyh(3,i)=blh(3,i)                   ! in meter
    call mercator(xyh(1,i),med,0.d0,0.d0,xyh(1,i))
  enddo
  !
  !! for this station ...
  call mercator(SITE%geod,med,0.d0,0.d0,xy)
  !
  !! read vmf3
  SITE%nvm=0
  ltim=.false.
  lfill(1:4)=.false.
  do while(.true.)
    read(lfnh,'(a)',iostat=ierr) line
    if(ierr.ne.0) exit
    if(line(1:1).eq.'!') then
      if(index(line,'Epoch').ne.0) then
        if(any(lfill(1:4).eqv..true.)) then
          write(*,'(a,i4,4i3,f7.2)') '***ERROR(vmf3_grid): lack of grid points ',iy,imon,id,ih,im,sec
          call exit(1)
        endif
        ltim=.true.
        oc=0.d0
        lfill(1:4)=.false.
        read(line(23:),*) iy,imon,id,ih,im,sec
        SITE%nvm=SITE%nvm+1
        if(SITE%nvm.gt.13) then           ! at most for 3 days
          write(*,'(a)') '***ERROR(vmf3_grid): too many epochs for vmf3 '
          call exit(1)
        endif
      endif
      cycle
    endif
    if(.not.ltim) cycle
  !
  !! a new epoch for vmf3
    read(line,*) lat,lon,temp(1:4)
    do i=1,4
      if(lat.eq.blh(1,i).and.lon.eq.blh(2,i)) then
        lfill(i)=.true.
        oc(1:4,i)=temp(1:4)
  !
  !! transfer hydrostatic delays to station height at grid points
        pres=oc(3,i)*ffun(blh(1,i)*PI/180.d0,blh(3,i)*1.d-3)/0.22768d-2
        pres=pres-1013.25d0*5.225d0*(1.d0-0.226d-4*blh(3,i))**4.225d0*0.226d-4*(SITE%geod(3)*1.d3-blh(3,i))
        oc(3,i)=0.22768d-2*pres/ffun(blh(1,i)*PI/180.d0,SITE%geod(3))
        if(all(lfill(1:4).eqv..true.)) exit
        cycle
      endif
    enddo
  !
  !! compute coefficients for interpolation
    if(all(lfill(1:4).eqv..true.)) then
      SITE%jdv (SITE%nvm)=modified_julday(id,imon,iy)
      SITE%sodv(SITE%nvm)=ih*3600.d0+im*60.d0+sec
      call vmf1_fit(xyh,oc,coef)
      do i=1,4
        SITE%vm1(i,SITE%nvm)=coef(1,i)+coef(2,i)*xy(1)*1.d-3+coef(3,i)*xy(2)*1.d-3
      enddo
      ltim=.false.
      lfill(1:4)=.false.
    endif
  enddo
else
  ilon=int((lon-0.5d0)/1.d0)+1          ! longitude grid interval 1 deg
  if(lat>89.5d0)then
    ilat=1
    blh(1,1:2)=89.5d0
    blh(1,3:4)=89.5d0
    geod_tmp(1)=89.5d0
  elseif(lat<=-89.5)then
    ilat=180
    blh(1,1:2)=-89.5d0
    blh(1,3:4)=-89.5d0
    geod_tmp(1)=-89.5d0
  endif
  blh(2,1)  =(ilon-1)*1.d0+0.5d0
  if(blh(2,1)>=360.d0) blh(2,1)=blh(2,1)-360.d0
  if(lon<0.5d0) blh(2,1)=359.5d0
  blh(2,3)  =blh(2,1)
  blh(2,2)  =ilon*1.d0+0.5d0
  if(blh(2,2)>=360.d0) blh(2,2)=blh(2,2)-360.d0
  if(lon<0.5d0) blh(2,2)=0.5d0
  blh(2,4)  =blh(2,2)
  !
  !! read ellipsoidal heights
  do i=1,ilat
    do j=1,360
      read(lfnv,'(f7.2)',err=100) ellh(j)
    enddo
  enddo
  blh(3,1)=ellh(ilon)
  blh(3,2)=ellh(ilon+1)
  blh(3,3)=ellh(ilon)
  blh(3,4)=ellh(ilon+1)
  !
  !! transform bl to xy using mercator projection
  med=(blh(2,1)+blh(2,2))/2.d0*PI/180.d0
  do i=1,4
    xyh(1,i)=blh(1,i)*PI/180.d0         ! to radian
    xyh(2,i)=blh(2,i)*PI/180.d0
    xyh(3,i)=blh(3,i)                   ! in meter
    call mercator(xyh(1,i),med,0.d0,0.d0,xyh(1,i))
  enddo
  !
  !! for this station ...
  geod_tmp(2)=SITE%geod(2)
  geod_tmp(3)=SITE%geod(3)
  call mercator(geod_tmp,med,0.d0,0.d0,xy)
  !
  !! read vmf3
  SITE%nvm=0
  ltim=.false.
  lfill(1:4)=.false.
  do while(.true.)
    read(lfnh,'(a)',iostat=ierr) line
    if(ierr.ne.0) exit
    if(line(1:1).eq.'!') then
      if(index(line,'Epoch').ne.0) then
        if(any(lfill(1:4).eqv..true.)) then
          write(*,'(a,i4,4i3,f7.2)') '***ERROR(vmf3_grid): lack of grid points ',iy,imon,id,ih,im,sec
          call exit(1)
        endif
        ltim=.true.
        oc=0.d0
        lfill(1:4)=.false.
        read(line(23:),*) iy,imon,id,ih,im,sec
        SITE%nvm=SITE%nvm+1
        if(SITE%nvm.gt.13) then           ! at most for 3 days
          write(*,'(a)') '***ERROR(vmf3_grid): too many epochs for vmf3 '
          call exit(1)
        endif
      endif
      cycle
    endif
    if(.not.ltim) cycle
  !
  !! a new epoch for vmf3
    read(line,*) lat,lon,temp(1:4)
    do i=1,2
      if(lat.eq.blh(1,i).and.lon.eq.blh(2,i)) then
        lfill(i)=.true.
        lfill(i+2)=.true.
        oc(1:4,i)=temp(1:4)
        oc(1:4,i+2)=temp(1:4)
  !
  !! transfer hydrostatic delays to station height at grid points
        pres=oc(3,i)*ffun(blh(1,i)*PI/180.d0,blh(3,i)*1.d-3)/0.22768d-2
        pres=pres-1013.25d0*5.225d0*(1.d0-0.226d-4*blh(3,i))**4.225d0*0.226d-4*(SITE%geod(3)*1.d3-blh(3,i))
        oc(3,i)=0.22768d-2*pres/ffun(blh(1,i)*PI/180.d0,SITE%geod(3))
    
        pres=oc(3,i+2)*ffun(blh(1,i+2)*PI/180.d0,blh(3,i+2)*1.d-3)/0.22768d-2
        pres=pres-1013.25d0*5.225d0*(1.d0-0.226d-4*blh(3,i+2))**4.225d0*0.226d-4*(SITE%geod(3)*1.d3-blh(3,i+2))
        oc(3,i+2)=0.22768d-2*pres/ffun(blh(1,i+2)*PI/180.d0,SITE%geod(3))
        if(all(lfill(1:4).eqv..true.)) exit
        cycle
      endif
    enddo
  !
  !! compute coefficients for interpolation
    if(all(lfill(1:4).eqv..true.)) then
      SITE%jdv (SITE%nvm)=modified_julday(id,imon,iy)
      SITE%sodv(SITE%nvm)=ih*3600.d0+im*60.d0+sec
      blh_tmp(1)=blh(1,2)
      blh_tmp(2)=blh(2,2)
      blh_tmp0(2)=geod_tmp(1)
      blh_tmp0(2)=geod_tmp(2)
      if(blh(2,2)<blh(2,1) .and. geod_tmp(2)<blh(2,1))then
        blh_tmp(2)=blh(2,2)+360
        blh_tmp0(2)=geod_tmp(2)+360
      endif
      if(blh(2,2)<blh(2,1) .and. geod_tmp(2)>=blh(2,1))then
        blh_tmp(2)=blh(2,2)+360
        blh_tmp0(2)=geod_tmp(2)
      endif
      do i=1,2
        blh_tmp(i)=blh_tmp(i)*PI/180.d0         ! to radian
        blh_tmp0(i)=blh_tmp0(i)*PI/180.d0
      enddo
      med=(blh_tmp(1)+blh_tmp(2))/2.d0*PI/180.d0
      call mercator(blh_tmp,med,0.d0,0.d0,xyh_tmp)
      call mercator(blh_tmp0,med,0.d0,0.d0,xy_tmp)
      ylen=xyh_tmp(2)-xyh(1,1)
      y1=(xyh_tmp(2)-xy_tmp(2))/ylen
      y2=(xy_tmp(2)-xyh(1,1))/ylen
      do i=1,4
        SITE%vm1(i,SITE%nvm)=y1*oc(i,1)+y2*oc(i,2)
      enddo
      ltim=.false.
      lfill(1:4)=.false.
    endif
  enddo  
endif

return
100 write(*,'(a)') '***ERROR(vmf3_grid): read file orography_ell_1x1 '
call exit(1)
end
