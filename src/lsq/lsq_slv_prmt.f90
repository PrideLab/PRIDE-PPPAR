!
!! lsq_slv_prmt.f90
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
!! Contributor: Maorong Ge, Jianghui Geng
!!
!!
!!
!! purpose   : solve for parameter in LSQ
!!
!! parameters: NM,PM -- parameter & normal equation struct
!!
!
subroutine lsq_slv_prmt(LCF, NM, PM)
  implicit none
  include '../header/const.h'
  include '../header/station.h'
  include 'lsqcfg.h'
  include 'lsq.h'

  type(lsqcfg) LCF
  type(prmt) PM(1:*)
  type(norm) NM
!
!! local
  integer*4 i, j
  real*8 det
!
!! copy upper half
  do i = 1, NM%imtx
    do j = i + 1, NM%imtx
      NM%norx(j, i) = NM%norx(i, j)
    end do
  end do
!
!! inverse
  call matinv(NM%norx, NM%nmtx, NM%imtx, det)
  if (det .eq. 0.d0) then
    write (*, '(a)') '***ERROR(lsq_slv_prmt): matrix inversion, singularity '
    call exit(1)
  end if
!
!! solve parameter
  do i = 1, NM%imtx
    det = 0.d0
    do j = 1, NM%imtx
      det = det + NM%norx(i, j)*NM%norx(j, NM%imtx + 1)
    end do
    PM(NM%iptp(i))%xcor = det
    NM%ltpl = NM%ltpl - PM(NM%iptp(i))%xcor*NM%norx(i, NM%imtx + 1)
  end do
!
!! sigma0
  NM%sig0 = dsqrt(dabs(NM%ltpl)/(NM%nobs - NM%nuk))
!
!! sigma of parameters
  do i = 1, NM%imtx
    PM(NM%iptp(i))%xsig = dsqrt(NM%norx(i, i))*NM%sig0
  end do

  return
end
