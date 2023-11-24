!
!! istrue.f90
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
logical*1 function istrue(flag, char)
  implicit none
  include 'data_flag.h'
  integer*4 flag, i, j, k
  character*(*) char

  istrue = .FALSE.

  if (char(1:2) .eq. 'ok' .or. char(1:2) .eq. 'OK') then
    if (iand(flag, z'FFFF0000') .eq. 0) istrue = .TRUE.
  else if (char(1:3) .eq. 'no4') then
    i = ibset(flag, flag_no4)
    istrue = i .eq. flag
  else if (char(1:3) .eq. 'lli') then
    i = ibset(flag, flag_lli)
    istrue = i .eq. flag
  else if (char(1:3) .eq. 'amb' .or. char(1:3) .eq. 'AMB') then
    if (iand(flag, z'FFFF0000') .eq. 0 .and. flag .ne. 0) istrue = .TRUE.
  else if (char(1:3) .eq. 'gap') then
    i = ibset(flag, flag_gap)
    istrue = i .eq. flag
  else if (char(1:4) .eq. 'good' .or. char(1:4) .eq. 'GOOD') then
    if (flag .eq. 0) istrue = .TRUE.
  else if (char(1:4) .eq. 'shrt') then
    i = ibset(flag, flag_shrt)
    istrue = i .eq. flag
  else if (char(1:5) .eq. 'lgbad') then
    i = ibset(flag, flag_lgbad)
    istrue = i .eq. flag
  else if (char(1:5) .eq. 'lwbad') then
    i = ibset(flag, flag_lwbad)
    istrue = i .eq. flag
  else if (char(1:5) .eq. 'pcbad') then
    i = ibset(flag, flag_pcbad)
    istrue = i .eq. flag
  else if (char(1:5) .eq. 'pc1ms') then
    i = ibset(flag, flag_pc1ms)
    istrue = i .eq. flag
  else if (char(1:5) .eq. 'bigsd') then
    i = ibset(flag, flag_bigsd)
    istrue = i .eq. flag
  else if (char(1:6) .eq. 'lwjump') then
    i = ibset(flag, flag_lwjump)
    istrue = i .eq. flag
  else if (char(1:6) .eq. 'lgjump') then
    i = ibset(flag, flag_lgjump)
    istrue = i .eq. flag
  else if (char(1:6) .eq. 'lwconn') then
    i = ibset(flag, flag_lwconn)
    istrue = i .eq. flag
  else if (char(1:6) .eq. 'lowele') then
    i = ibset(flag, flag_lowele)
    istrue = i .eq. flag
  else if (char(1:6) .eq. 'nodata') then
    i = ibset(flag, flag_nodata)
    istrue = i .eq. flag
  else if (char(1:7) .eq. 'lccheck') then
    i = ibset(flag, flag_lccheck)
    istrue = i .eq. flag
  else if (char(1:8) .eq. 'range_ok') then
    i = ibset(flag, flag_lwjump)
    j = ibset(flag, flag_lgjump)
    k = ibset(flag, flag_pcbad)
    if (iand(flag, z'FFFF0000') .eq. 0 .or. &
        (k .eq. flag .and. flag .ne. i .and. flag .ne. j)) istrue = .TRUE.
  else
    call exit(1)
  endif

  return
end
