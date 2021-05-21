!
!! lsq_dop.f90
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
!
subroutine lsq_dop(jd,sod,LCF,OB,nbb)
implicit none
include '../header/const.h'
include '../header/rnxobs.h'
include 'lsqcfg.h'
integer*4 jd
real*8 sod
type(lsqcfg) LCF
type(rnxobr) OB
integer*4 i,j,k,NROW,NCOLB,NCOLA
integer*4 isat,isat_real,ipar
real*8 wrng(MAXSAT),amat(MAXSAT,4),amatt(4,MAXSAT),nbb(4,4)
real*8 SUM,det
! R
integer*4 frequency_glo_nu,prn_int
real*8 :: FREQ1_R(-50:50),FREQ2_R(-50:50)

call frequency_glonass(FREQ1_R,FREQ2_R)

isat_real=0
amat=0.d0
amatt=0.d0
nbb=0.d0
do isat=1,LCF%nprn
  if(OB%omc(isat,1).eq.0.d0.or.OB%omc(isat,3).eq.0.d0) cycle
  isat_real=isat_real+1
  do ipar=1,OB%npar
    if(OB%pname(ipar)(1:5) .eq. 'STAPX')then
      do i = 1, 3
        if(LCF%prn(isat)(1:1).eq.'G')then
          amat(isat_real,i)=OB%amat(ipar+i-1,isat)*VLIGHT/FREQ1_G
          amatt(i,isat_real)=OB%amat(ipar+i-1,isat)*VLIGHT/FREQ1_G
        elseif(LCF%prn(isat)(1:1).eq.'R')then
          read(LCF%prn(isat),'(1x,i2)') prn_int
          frequency_glo_nu=OB%glschn(prn_int)
          amat(isat_real,i)=OB%amat(ipar+i-1,isat)*VLIGHT/FREQ1_R(frequency_glo_nu)
          amatt(i,isat_real)=OB%amat(ipar+i-1,isat)*VLIGHT/FREQ1_R(frequency_glo_nu)
        elseif(LCF%prn(isat)(1:1).eq.'E')then
          amat(isat_real,i)=OB%amat(ipar+i-1,isat)*VLIGHT/FREQ1_E
          amatt(i,isat_real)=OB%amat(ipar+i-1,isat)*VLIGHT/FREQ1_E
        elseif(LCF%prn(isat)(1:1).eq.'C')then
          amat(isat_real,i)=OB%amat(ipar+i-1,isat)*VLIGHT/FREQ1_C
          amatt(i,isat_real)=OB%amat(ipar+i-1,isat)*VLIGHT/FREQ1_C
        elseif(LCF%prn(isat)(1:1).eq.'J')then
          amat(isat_real,i)=OB%amat(ipar+i-1,isat)*VLIGHT/FREQ1_J
          amatt(i,isat_real)=OB%amat(ipar+i-1,isat)*VLIGHT/FREQ1_J
        end if
      enddo
    endif
  end do
  amat(isat_real,4)=1.d0
  amatt(4,isat_real)=1.d0
end do

if(isat_real .gt. 4)then
  NROW=4
  NCOLB=4
  NCOLA=isat_real
  DO i=1,NROW
    DO j=1,NCOLB
      SUM = 0.D0
      DO k=1,NCOLA
    	SUM = SUM + amatt(i,k)*amat(k,j)
      enddo
      nbb(i,j) = SUM
    enddo
  enddo
  call matinv(nbb,4,4,det)
endif
end