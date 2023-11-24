!
!! prn_matbld.f90
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
!! Contributor: Songfeng Yang, Jihang Lin
!! 
!!
!!
subroutine prn_matbld(prn_mat)
  implicit none
  include '../header/const.h'
  character*3 prn_mat(MAXSAT)
  integer*4 i0, i, j
  j = 0
  do i0 = 1, MAXSYS
    do i = 1, MAXSAT_SYS(i0)
      j = j + 1
      if (j .gt. MAXSAT) exit
      write (prn_mat(j), '(a1,i0.2)') GNSS_PRIO(i0:i0), i
    end do
  end do
  return
end subroutine
