!
!! lagrange_coeff.f90
!!
!!    Copyright (C) 2023 by Wuhan University
!!
!!    This program belongs to PRIDE PPP-AR which is an open source software: 
!!    you can redistribute it and/or modify it under the terms of the GNU 
!!    General Public License (version 3) as published by the Free Software Foundation.
!!
!!    This program is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!!    GNU General Public License (version 3) for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with this program. If not, see <https://www.gnu.org/licenses/>.
!!
!! Contributor: Maorong Ge, Shuyin Mao
!! 
!!
!!
!! purpose  : compute coefficients of the lagrange interpolation
!!
!! parameter: norder  -- degree of the interpolation
!!            dt      -- time for interpolation
!!            lvel    -- whether to compute velocity coefficients
!!            coeff   -- lagrange coefficients
!!
!
  subroutine lagrange_coeff(norder, dt, lvel, coeff)
  implicit none

  integer*4 norder
  real*8 dt
  logical*1 lvel
  real*8 coeff(1:*)

  !
  !! local 
  integer*4, parameter ::  MAXORD=30
  integer*4 i,j,k
  logical*1 first
  real*8 s
  real*8 coeff0(0:MAXORD,0:MAXORD)

  data first/.true./
  save first, coeff0

  if(first) then
    first = .false.
    do i = 0, MAXORD
      do j = 0, i 
        coeff0(i,j) = 1.d0
        do k = 0, i
          if(k.eq.j) cycle
          coeff0(i,j)=coeff0(i,j)*(j-k)*1.d0 
        enddo
      enddo
    enddo
  endif

  do i = 0, norder
    coeff(i+1) = 1.d0
    !
    !! compute pos coefficient
    do j = 0, norder
      if (j.eq.i) cycle
      coeff(i+1)=coeff(i+1)*(dt-j)
    enddo
    coeff(i+1)=coeff(i+1)/coeff0(norder,i)
    !
    !! compute velocity coefficient
    if (lvel) then
      coeff(norder+1+i+1) = 0.d0
      do j = 0, norder
        if (j.eq.i) cycle
        s = 1.d0
        do k = 0, norder
          if (k.eq.j .or. k.eq.i) cycle
          s = s*(dt-k)
        enddo
        coeff(norder+1+i+1)=coeff(norder+1+i+1)+s
      enddo
      coeff(norder+1+i+1)=coeff(norder+1+i+1)/coeff0(norder,i)
    endif
  enddo

  return
  end
