!
!! sort_invx.f90
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
!! Contributor: Jianghui Geng
!! 
!!
!!
!! purpose  : sort inversed normal matrix (lower triangular part)
!! parameter:
!!    input : MD%-- diffrenced ambiguities
!!            maxdev -- maximum deviation
!!            maxsig -- maximum sigma
!!    output: QN%-- inversed normal matrix
!
subroutine sort_invx(MD, QN, invx, bias, q22)
  implicit none
  include 'ambssd.h'
  include 'invnor.h'

  type(ambssd) MD(1:*)
  type(invm) QN
  real*8 bias(1:*), q22(1:*), invx(1:*)
!
!! local
  integer*4 i, j, k, ifnd, ntot
  real*8 dump
  type(ambssd) MDX
!
!! find fixable MD ambiguities
  ntot = QN%nxyz + QN%ndam
  do i = QN%nxyz + 1, ntot
    if (MD(i - QN%nxyz)%id .eq. 1) cycle
    ifnd = 1
    do j = i + 1, ntot
      if (MD(j - QN%nxyz)%id .eq. 1) then
        ifnd = 0
        MDX = MD(i - QN%nxyz)
        MD(i - QN%nxyz) = MD(j - QN%nxyz)
        MD(j - QN%nxyz) = MDX
!
!! move elements in whole matrix
        do k = 1, i - 1                         ! column
          dump = invx(QN%idq(k) + i)
          invx(QN%idq(k) + i) = invx(QN%idq(k) + j)
          invx(QN%idq(k) + j) = dump
        enddo
        do k = j + 1, ntot                      ! row
          dump = invx(QN%idq(i) + k)
          invx(QN%idq(i) + k) = invx(QN%idq(j) + k)
          invx(QN%idq(j) + k) = dump
        enddo
        do k = i + 1, j - 1                       ! row
          dump = invx(QN%idq(i) + k)
          invx(QN%idq(i) + k) = invx(QN%idq(k) + j)
          invx(QN%idq(k) + j) = dump
        enddo
        dump = invx(QN%idq(i) + i)
        invx(QN%idq(i) + i) = invx(QN%idq(j) + j)
        invx(QN%idq(j) + j) = dump
        exit
      endif
    enddo
    if (ifnd .ne. 0) exit
  enddo
!
!! yield Q22 & bias
  k = 0
  do j = ntot - QN%ncad + 1, ntot
    bias(j - ntot + QN%ncad) = MD(j - QN%nxyz)%rnl
    do i = j, ntot
      k = k + 1
      q22(k) = invx(QN%idq(j) + i)
    enddo
  enddo

  return
end
