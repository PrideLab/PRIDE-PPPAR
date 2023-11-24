!
!! run_tim.f90
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
!!
!! purpose  : show running time
!! parameter:
!
character*10 function run_tim()
  implicit none

!
!! local
  integer*4 date_time(8)
  character*8 :: real_clock1 = ' '
  character*10 :: real_clock2 = ' '
  character*5 :: real_clock3 = ' '

  call date_and_time(real_clock1, real_clock2, real_clock3, date_time)
  write (run_tim, '(a1,2(i2.2,1h:),i2.2,a1)') ' ', date_time(5:7), ' '

  return
end
