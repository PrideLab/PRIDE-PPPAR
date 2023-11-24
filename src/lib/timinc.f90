!
!! timinc.f90
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
!! Contributor: Maorong Ge, Jianghui Geng, Songfeng Yang
!! 
!!
subroutine timinc(jd, sec, delt, jd1, sec1)
  implicit none

  integer*4 jd, inc, nss, jd1
  real*8 sec, delt, sec1
  data nss/86400/

  sec1 = sec + delt
  inc = int(sec1/nss)
  sec1 = sec1 - inc*nss
  jd1 = jd + inc
  if (sec1 .ge. 0) return
  jd1 = jd1 - 1
  sec1 = nss + sec1
  return
end

