!
!! set_flag.f90
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
!! Contributor: Maorong Ge, Jianghui Geng, Songfeng Yang
!! 
!!
!
integer*4 function set_flag(flag, char)
  implicit none
  include 'data_flag.h'
  integer*4 flag, i, j
  character*(*) char

  j = flag
  if (char(1:2) .eq. 'ok' .or. char(1:2) .eq. 'OK') then
    do i = 16, 31
      j = ibclr(j, i)
    enddo
  else if (char(1:4) .eq. 'good' .or. char(1:4) .eq. 'GOOD') then
    do i = 0, 31
      j = ibclr(j, i)
    enddo
  else if (char(1:6) .eq. 'lgjump') then
    j = ibset(j, flag_lgjump)
  else if (char(1:6) .eq. 'lwjump') then
    j = ibset(j, flag_lwjump)
  else if (char(1:5) .eq. 'lgbad') then
    j = ibset(j, flag_lgbad)
  else if (char(1:5) .eq. 'lwbad') then
    j = ibset(j, flag_lwbad)
  else if (char(1:5) .eq. 'pcbad') then
    j = ibset(j, flag_pcbad)
  else if (char(1:5) .eq. 'pc1ms') then
    j = ibset(j, flag_pc1ms)
  else if (char(1:6) .eq. 'lggood') then
    j = ibclr(j, flag_lgjump)
  else if (char(1:6) .eq. 'lwgood') then
    j = ibclr(j, flag_lwjump)
  else if (char(1:6) .eq. 'nodata') then
    j = ibset(j, flag_nodata)
  else if (char(1:6) .eq. 'lwconn') then
    j = ibset(j, flag_lwconn)
  else if (char(1:3) .eq. 'lli') then
    j = ibset(j, flag_lli)
  else if (char(1:3) .eq. 'no4') then
    j = ibset(j, flag_no4)
  else if (char(1:3) .eq. 'gap') then
    j = ibset(j, flag_gap)
  else if (char(1:4) .eq. 'shrt') then
    j = ibset(j, flag_shrt)
  else if (char(1:6) .eq. 'lowele') then
    j = ibset(j, flag_lowele)
  else if (char(1:5) .eq. 'bigsd') then
    j = ibset(j, flag_bigsd)
  else if (char(1:7) .eq. 'lccheck') then
    j = ibset(j, flag_lccheck)
  else
    call exit(1)
  endif
  set_flag = j
  return
end

