!!
!! vmf1_fit.f90
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
!! purpose  : configure VMF1 for static station. read new vmf format 
!!            by J. Boehm at http://www.hg.tuwien.ac.at/~ecmwf1/
!! parameter:
!! input    : flnvmf -- vmf1 file
!! output   : SITE   -- station struct
!! note     : vmf1 files should be concatenated to one file
!!

subroutine vmf1_grid(flnvmf, SITE)


implicit none
include '../header/const.h'
include '../header/station.h'

character*20 flnvmf
type(station) SITE
!
!! local
logical*1 lfirst,lfill(4),ltim
integer*4 i,lfnh,lfnv,ilat,ilon,iy,imon,id,ih,im,ellh(145),ierr
real*8 lat,lon,med,blh(3,4),xyh(3,4),oc(4,4),temp(4),coef(3,4),xy(2),pres,sec
character*256 line
!
!! function used
integer*4 get_valid_unit,modified_julday
real*8 ffun

data lfirst/.true./
save lfirst,lfnh,lfnv

if(lfirst) then
  lfirst=.false.
!! grid vmf1 file
  lfnh=get_valid_unit(10)
  open(lfnh,file=flnvmf,status='old',iostat=ierr)
  if(ierr.ne.0) then
    write(*,'(2a)') '***ERROR(vmf1_grid): open file ',flnvmf(1:len_trim(flnvmf))
    call exit(1)
  endif
!! grid ellipsoidal height file
  lfnv=get_valid_unit(10)
  open(lfnv,file='orography_ell',status='old',iostat=ierr)
  if(ierr.ne.0) then
    write(*,'(a)') '***ERROR(vmf1_grid): open file orography_ell '
    call exit(1)
  endif
endif
rewind lfnh
rewind lfnv
read(lfnv,'(a)') line      ! skip single header line
!
!! determine which 4 grid points should be used
lat=SITE%geod(1)*180.d0/PI            ! to degrees
lon=SITE%geod(2)*180.d0/PI
ilat=int((90.d0-lat)/2.0d0)+1         ! latitude  grid interval 2 deg
ilon=int(lon/2.5d0)+1                 ! longitude grid interval 2.5 deg
blh(1,1:2)=90.d0-(ilat-1)*2.d0
blh(1,3:4)=90.d0-ilat*2.d0
blh(2,1)  =(ilon-1)*2.5d0
if(blh(2,1).eq.360.d0) blh(2,1)=0.d0
blh(2,3)  =blh(2,1)
blh(2,2)  =ilon*2.5d0
if(blh(2,2).eq.360.d0) blh(2,2)=0.d0
blh(2,4)  =blh(2,2)
!
!! read ellipsoidal heights
do i=1,ilat
  read(lfnv,'(14(10i8,/),5i8)',err=100) ellh(1:145)
enddo
blh(3,1)=ellh(ilon)
blh(3,2)=ellh(ilon+1)
read(lfnv,'(14(10i8,/),5i8)',err=100) ellh(1:145)
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
!! read vmf1
SITE%nvm=0
ltim=.false.
lfill(1:4)=.false.
do while(.true.)
  read(lfnh,'(a)',iostat=ierr) line
  if(ierr.ne.0) exit
  if(line(1:1).eq.'!') then
    if(index(line,'Epoch').ne.0) then
      if(any(lfill(1:4).eqv..true.)) then
        write(*,'(a,i4,4i3,f7.2)') '***ERROR(vmf1_grid): lack of grid points ',iy,imon,id,ih,im,sec
        call exit(1)
      endif
      ltim=.true.
      oc=0.d0
      lfill(1:4)=.false.
      read(line(23:),*) iy,imon,id,ih,im,sec
      SITE%nvm=SITE%nvm+1
      if(SITE%nvm.gt.13) then           ! at most for 3 days
        write(*,'(a)') '***ERROR(vmf1_grid): too many epochs for vmf1 '
        call exit(1)
      endif
    endif
    cycle
  endif
  if(.not.ltim) cycle
!
!! a new epoch for vmf1
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

return
100 write(*,'(a)') '***ERROR(vmf1_grid): read file orography_ell '
call exit(1)
end
