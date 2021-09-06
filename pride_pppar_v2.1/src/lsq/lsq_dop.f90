!
!! lsq_dop.f90
!!
!!    Copyright (C) 2021 by Wuhan University
!!
!!    This program belongs to PRIDE PPP-AR which is an open source software:
!!    you can redistribute it and/or modify it under the terms of the GNU
!!    General Public License (version 3) as published by the Free Software Foundation.
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
integer*4 isat,isat_real
real*8 amat(MAXSAT,4),amatt(4,MAXSAT),nbb(4,4),sinel,cosel
real*8 SUM,det

isat_real=0
amat=0.d0
amatt=0.d0
nbb=0.d0
do isat=1,LCF%nprn
  if(OB%omc(isat,1).eq.0.d0.or.OB%omc(isat,3).eq.0.d0) cycle
  isat_real=isat_real+1
  sinel=dsin(OB%elev(isat))
  cosel=dcos(OB%elev(isat))
  amat(isat_real,1)=cosel*dsin(OB%azim(isat))
  amat(isat_real,2)=cosel*dcos(OB%azim(isat))
  amat(isat_real,3)=sinel
  amat(isat_real,4)=1.d0
  amatt(1,isat_real)=cosel*dsin(OB%azim(isat))
  amatt(2,isat_real)=cosel*dcos(OB%azim(isat))
  amatt(3,isat_real)=sinel
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