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
!! purpose  : determine vmf1 interpolation coefficients
!! parameter:
!! input    : xyh  -- projected coordinates
!!            oc   -- vmf1 values
!! output   : coef -- vmf1 coefficients
!!

subroutine vmf1_fit(xyh, oc, coef)
implicit none
include '../header/const.h'

real*8 xyh(3,4),oc(4,4),coef(3,4)
!
!! local
integer*4 i,j,k,m
real*8 a(3),c(3,3),l(3),det
!
!! four set of coefficients
do m=1,4
  c=0.d0
  l=0.d0
  do i=1,4
    a(1)=1.d0
    a(2)=xyh(1,i)*1.d-3     ! to km
    a(3)=xyh(2,i)*1.d-3
!
!! add to normal equation
    do j=1,3
      do k=j,3
        c(j,k)=c(j,k)+a(j)*a(k)
      enddo
      l(j)=l(j)+a(j)*oc(m,i)
    enddo
  enddo
!
!! make a full normal matrix
  do j=1,3
    do k=j+1,3
      c(k,j)=c(j,k)
    enddo
  enddo
!
!! invert matrix
  call matinv(c,3,3,det)
  if(det.eq.0.d0) then
    write(*,'(a)') '***ERROR(vmf1_grid): matrix singularity '
    call exit(1)
  endif
!
!! retrieve parameter estimates
  do j=1,3
    coef(j,m)=0.d0
    do k=1,3
      coef(j,m)=coef(j,m)+c(j,k)*l(k)
    enddo
  enddo
enddo

return
end
