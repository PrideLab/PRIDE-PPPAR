!
!! map_invnormal.f90
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
!! Contributor: Jianghui Geng, Jihang Lin
!!
!!
!!
!! purpose  : map inversed normal matrix from one-way to multiple-difference ones
!! parameter:
!!    input : SD -- multiple difference struct (SD)
!!    output: QN -- inversed normal matrix (lower triangular part)
!
subroutine map_invnormal(SD, QN, invx)
  implicit none
  include '../header/const.h'
  include 'ambssd.h'
  include 'invnor.h'

! common
  integer*4     idxfrq(MAXSYS, 2)
  common        idxfrq
! parameter
  real*8        invx(1:*)
  type(ambssd)  SD(1:*)
  type(invm)    QN
! local
  integer*4     i0, i, j, k, l, ir, ic, ierr
  real*8        op(2), rl, rr
  real*8, pointer :: hp(:, :)
  real*8        f1(MAXSYS), f2(MAXSYS)

  data op/1.d0, -1.d0/
  save op

!
!! frequency
  do i0 = 1, MAXSYS
    f1(i0) = FREQ_SYS(idxfrq(i0, 1), i0)
    if (f1(i0) .eq. 0.d0) goto 100
    f2(i0) = FREQ_SYS(idxfrq(i0, 2), i0)
    if (f2(i0) .eq. 0.d0) goto 100
  end do
!
!! initialize
  allocate (hp(QN%indp, QN%indp + QN%nxyz), stat=ierr)
  if (ierr .ne. 0) then
    write (*, '(a)') '***ERROR(map_invnormal): memory allocation hp '
    call exit(1)
  end if
!
!! M * Q21
  do i = 1, QN%indp
    if (SD(i)%sys .ne. '3') then
      i0 = index(GNSS_PRIO, SD(i)%sys)
    else
      i0 = index(GNSS_PRIO, 'C')
    end if
    rl = (f1(i0) + f2(i0))/VLIGHT
    do j = 1, QN%nxyz
      hp(i, j) = 0.d0
      do k = 1, 2
        hp(i, j) = hp(i, j) + op(k) * invx(QN%idq(j) + SD(i)%ipt(k)) * rl
      end do
    end do
  end do
  do j = 1, QN%nxyz
    do i = 1, QN%ntot - QN%nxyz
      invx(QN%idq(j) + QN%nxyz + i) = 0.d0
      if (i .le. QN%indp) invx(QN%idq(j) + QN%nxyz + i) = hp(i, j)
    end do
  end do
!
!! M * Q22 * Mt
  do i = 1, QN%indp
    if (SD(i)%sys .ne. '3') then
      i0 = index(GNSS_PRIO, SD(i)%sys)
    else
      i0 = index(GNSS_PRIO, 'C')
    end if
    rl = (f1(i0) + f2(i0))/VLIGHT
    do j = 1, i
      if (SD(j)%sys .ne. '3') then
        i0 = index(GNSS_PRIO, SD(j)%sys)
      else
        i0 = index(GNSS_PRIO, 'C')
      end if
      rr = (f1(i0) + f2(i0))/VLIGHT
      hp(i, j) = 0.d0
      do k = 1, 2
        do l = 1, 2
          ir = max(SD(i)%ipt(k), SD(j)%ipt(l))
          ic = min(SD(i)%ipt(k), SD(j)%ipt(l))
          hp(i, j) = hp(i, j) + op(k)*invx(QN%idq(ic) + ir)*op(l)*rl*rr
        end do
      end do
    end do
  end do
  do j = QN%nxyz + 1, QN%ntot
    do i = j, QN%ntot
      invx(QN%idq(j) + i) = 0.d0
      if (i .le. QN%nxyz + QN%indp) invx(QN%idq(j) + i) = hp(i - QN%nxyz, j - QN%nxyz)
    end do
  end do
  deallocate (hp)
  return

100 continue
  write (*, '(a,5(1x,a,2i1))') '***ERROR(lsq_add_obs): invalid frequency number:', &
    'G', idxfrq(index(GNSS_PRIO, 'G'), 1:2), &
    'R', idxfrq(index(GNSS_PRIO, 'R'), 1:2), &
    'E', idxfrq(index(GNSS_PRIO, 'E'), 1:2), &
    'C', idxfrq(index(GNSS_PRIO, 'C'), 1:2), &
    'J', idxfrq(index(GNSS_PRIO, 'J'), 1:2)
  call exit(1)
end subroutine
