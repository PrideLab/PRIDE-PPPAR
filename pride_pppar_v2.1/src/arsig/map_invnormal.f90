!
!! map_invnormal.f90
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
!! Contributor: Jianghui Geng
!! 
!!
!!
!! purpose  : map inversed normal matrix from one-way to multiple-difference ones
!! parameter:
!!    input : MD -- multiple difference struct (SD)
!!    output: QN -- inversed normal matrix (lower triangular part)
!
subroutine map_invnormal(MD, QN, invx)
  implicit none
  include '../header/const.h'
  include 'ambssd.h'
  include 'invnor.h'

  real*8 invx(1:*)
  type(ambssd) MD(1:*)
  type(invm) QN
!
!! local
  integer*4 i, j, k, l, ir, ic, ierr
  real*8 op(2),rl,rr
  real*8, pointer :: hp(:, :)

  data op/1.d0, -1.d0/
  save op

  allocate (hp(QN%indp, QN%indp + QN%nxyz), stat=ierr)
  if (ierr .ne. 0) then
    write (*, '(a)') '***ERROR(map_invnormal): memory allocation hp '
    call exit(1)
  endif
!
!! M * Q21
  do i = 1, QN%indp
    if(MD(i)%sys .eq. 'G') rl=(FREQ1_G+FREQ2_G)/VLIGHT
    if(MD(i)%sys .eq. 'E') rl=(FREQ1_E+FREQ2_E)/VLIGHT
    if((MD(i)%sys .eq. 'C') .or. (MD(i)%sys .eq. '3')) rl=(FREQ1_C+FREQ2_C)/VLIGHT
    if(MD(i)%sys .eq. 'J') rl=(FREQ1_J+FREQ2_J)/VLIGHT
    do j = 1, QN%nxyz
      hp(i, j) = 0.d0
      do k = 1, 2
        hp(i, j) = hp(i, j) + op(k)*invx(QN%idq(j) + MD(i)%ipt(k))*rl
      enddo
    enddo
  enddo
  do j = 1, QN%nxyz
    do i = 1, QN%ntot - QN%nxyz
      invx(QN%idq(j) + QN%nxyz + i) = 0.d0
      if (i .le. QN%indp) invx(QN%idq(j) + QN%nxyz + i) = hp(i, j)
    enddo
  enddo
!
!! M * Q22 * Mt
  do i = 1, QN%indp
    if(MD(i)%sys .eq. 'G') rl=(FREQ1_G+FREQ2_G)/VLIGHT
    if(MD(i)%sys .eq. 'E') rl=(FREQ1_E+FREQ2_E)/VLIGHT
    if((MD(i)%sys .eq. 'C') .or. (MD(i)%sys .eq. '3')) rl=(FREQ1_C+FREQ2_C)/VLIGHT
    if(MD(i)%sys .eq. 'J') rl=(FREQ1_J+FREQ2_J)/VLIGHT
    do j = 1, i
      if(MD(j)%sys .eq. 'G') rr=(FREQ1_G+FREQ2_G)/VLIGHT
      if(MD(j)%sys .eq. 'E') rr=(FREQ1_E+FREQ2_E)/VLIGHT
      if((MD(j)%sys .eq. 'C') .or. (MD(j)%sys .eq. '3')) rr=(FREQ1_C+FREQ2_C)/VLIGHT
      if(MD(j)%sys .eq. 'J') rr=(FREQ1_J+FREQ2_J)/VLIGHT
      hp(i, j) = 0.d0
      do k = 1, 2
        do l = 1, 2
          ir = max(MD(i)%ipt(k), MD(j)%ipt(l))
          ic = min(MD(i)%ipt(k), MD(j)%ipt(l))
          hp(i, j) = hp(i, j) + op(k)*invx(QN%idq(ic) + ir)*op(l)*rl*rr
        enddo
      enddo
    enddo
  enddo
  do j = QN%nxyz + 1, QN%ntot
    do i = j, QN%ntot
      invx(QN%idq(j) + i) = 0.d0
      if (i .le. QN%nxyz + QN%indp) invx(QN%idq(j) + i) = hp(i - QN%nxyz, j - QN%nxyz)
    enddo
  enddo
  deallocate (hp)

  return
end
