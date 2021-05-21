!
!! calstech.f90
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
!! Contributor: Shaoming Xin, Jianghui Geng, Songfeng Yang
!! 
!!
!!
! introduce ionosphere observations
! cal the high order ionosphere correction
subroutine calstech(jd,sod,LCF,isat,SAT,IM,OB,stec)
implicit none
include '../header/const.h'
include '../header/satellite.h'
include '../header/rnxobs.h'
include '../header/ionex.h'
!
!! ion model
include 'lsq.h'
include 'lsqcfg.h'

type(lsqcfg)   LCF
type(satellite)SAT(1:*)
type(ionex)    IM
type(rnxobr)   OB
integer*4      jd
real*8         sod
real*8         stec(4)
!
!! local
logical*1      lfirst
integer*4      i,j,ib,ie,isat,iion,ipar,jdint
real*8         tint,ippmap,lon1,lon2,sec2rad,atim(maxsat),aiond(maxsat),aionrms(maxsat)
real*8         maxion,freqt,iond,iondat1,iondat2,ionrms,ionrms1,ionrms2
real*8         u1,u2,u3,udb,udt
!
!! function used
integer*4      pointer_string
real*8         taiutc,timdif

integer*4 frequency_glo_nu,prn_int
real*8 :: FREQ1_R(-50:50),FREQ2_R(-50:50)

data lfirst,maxion/.true.,30.d0/
save lfirst,maxion,sec2rad,atim,aiond,aionrms

call frequency_glonass(FREQ1_R,FREQ2_R)
!
!! initialization
if(lfirst) then
  lfirst=.false.
  sec2rad=pi/43200.d0
  atim=0.d0
  aiond=0.d0
  aionrms=0.d0
endif
!
!! which ionosphere product
ib=0
ie=0
lon1=0.d0
lon2=0.d0
if(IM%ctm(1:3).eq.'UTC') then
  call timinc(jd,sod,19.d0-taiutc(jd),jdint,tint)
else ! GPS time
  jdint=jd
  tint =sod
endif
!! look for closest maps
do j=1,IM%bmap
  if(timdif(jdint,tint,IM%mjd(j),IM%sod(j)).le.0.d0) exit
enddo
  
  if(j.eq.1) then
    ib=1
    ie=0
    if(IM%lrot) lon1=timdif(jdint,tint,IM%mjd(ib),IM%sod(ib))*sec2rad
    !write(*,'(a,i5,f8.1)') '###WARNING(lsq_add_tecmap): time outside of map ',jdint,tint
  else if(j.eq.IM%bmap+1) then
    ib=IM%bmap
    ie=0
    if(IM%lrot) lon1=timdif(jdint,tint,IM%mjd(ib),IM%sod(ib))*sec2rad
    !write(*,'(a,i5,f8.1)') '###WARNING(lsq_add_tecmap): time outside of map ',jdint,tint
  else
    ib=j-1
    ie=j
    if(IM%lrot) then
      lon1=timdif(jdint,tint,IM%mjd(ib),IM%sod(ib))*sec2rad
      lon2=timdif(IM%mjd(ie),IM%sod(ie),jdint,tint)*sec2rad
    endif
  endif
!
!! modified SLM as used by CODE
if(IM%cmf.eq.'COSZ') then
  ippmap=1.d0/dcos(OB%zeni(isat))
else
  ippmap=1.d0/dcos(dasin(6371.d0/6877.7d0*dsin(0.9782d0*(pi/2.d0-OB%elev(isat)))))
endif
!
!! interpolate ionosphere
if(ie.ne.0) then
  call interp_ionex(IM,ib,OB%ilat(isat),OB%ilon(isat)+lon1,iondat1,ionrms1)
  call interp_ionex(IM,ie,OB%ilat(isat),OB%ilon(isat)-lon2,iondat2,ionrms2)
  if(iondat1.eq.0.d0.or.iondat2.eq.0.d0.or.ionrms1.eq.0.d0.or.ionrms2.eq.0.d0) then
    iond=0.d0
    ionrms=0.d0
  else
    iond=iondat1+(iondat2-iondat1)*timdif(jdint,tint,IM%mjd(ib),IM%sod(ib))/IM%dintv ! vertical TEC
    ionrms=ionrms1+(ionrms2-ionrms1)*timdif(jdint,tint,IM%mjd(ib),IM%sod(ib))/IM%dintv
  endif
else
  call interp_ionex(IM,ib,OB%ilat(isat),OB%ilon(isat)+lon1,iond,ionrms)
endif
!
!! slant delay and rms
if(iond.eq.0.d0.or.ionrms.eq.0.d0) then
  iond=10.d0   ! a priori assumed
  ionrms=10.d0
else
  do i=1,2
    if(OB%prn(isat)(1:1) .eq. 'G')then
      if(i.eq.1) freqt=freq1_G
      if(i.eq.2) freqt=freq2_G
    elseif(OB%prn(isat)(1:1) .eq. 'R')then
      read(OB%prn(isat),'(1x,i2)') prn_int
      frequency_glo_nu=OB%glschn(prn_int)
      if(i.eq.1) freqt=freq1_R(frequency_glo_nu)
      if(i.eq.2) freqt=freq2_R(frequency_glo_nu)
    elseif(OB%prn(isat)(1:1) .eq. 'E')then
      if(i.eq.1) freqt=freq1_E
      if(i.eq.2) freqt=freq2_E
    elseif(OB%prn(isat)(1:1) .eq. 'C')then
      if(i.eq.1) freqt=freq1_C
      if(i.eq.2) freqt=freq2_C
    elseif(OB%prn(isat)(1:1) .eq. 'J')then
      if(i.eq.1) freqt=freq1_J
      if(i.eq.2) freqt=freq2_J
    endif
    u1=0.d0;u2=0.d0;u3=0.d0;udb=0.d0;udt=0.d0
    u1=40.3d0*1.d16/freqt**2
    u2=7527.d0*1.d16*vlight*(OB%ib0(isat)*1.d-9*dcos(OB%itheta(isat)))/freqt**3
    u3=2437.d0*1.d16*(((20.d0-6.d0)*1.d12)/((4.55d0-1.38d0)*1.d18)*(iond*1.d16-4.55d0*1.d18)+20.d0*1.d12)*0.66d0/freqt**4
    udb=2.9344d0*1.d5*(1.d0/dsqrt(1.d0-0.826d0*dcos(OB%elev(isat))**2)-1.d0)/(1.d-24*freqt**4)*(ippmap*iond)
    udt=58.689d0*1.d16*(1.d0/dsqrt(1.d0-0.826d0*dcos(OB%elev(isat))**2)-1.d0)/(1.d-12*freqt**4)*(ippmap*iond)
    stec(i)=ippmap*iond*(u2/2.d0+u3/3.d0-udb+udt)      ! from TECU to meters on L1
    stec(i+2)=ippmap*iond*(u2+u3+udb+udt)  ! from TECU to meters on P1
  enddo
endif
!
return
end
