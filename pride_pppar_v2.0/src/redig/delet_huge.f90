!
!! delet_huge.f90
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
!!
!! purpose  : delete huge residuals
!! parameter:
!!    input : lupd -- update rhd or not
!!            xres -- threshold value for huge residuals
!!            nepo -- # of residuals
!!            resi -- residuals
!!            trsi -- time tags for each epoch
!!    output: flag -- flags of each epoch
!
subroutine delet_huge(lupd, xres, nepo, flag, resi, trsi)
  implicit none
  include '../header/const.h'
  include 'data_flag.h'

  logical*1 lupd
  integer*4 nepo, flag(1:*)
  real*8 xres, resi(1:*)
  character*27 trsi(1:*)
!
!! local
  integer*4 i, j
!
!! function used
  logical*1 istrue

  do i = 1, nepo
    if (.not. istrue(flag(i), 'OK')) cycle
    if (dabs(resi(i)) .gt. xres) then
      if (istrue(flag(i), 'AMB')) then
        if (i .lt. nepo) then
          call find_flag(i + 1, nepo, flag, 'OK', j)
          if (j .gt. 0 .and. istrue(flag(j), 'GOOD')) flag(j) = flag(i)
        endif
      endif
      lupd = .true.
      flag(i) = DELBAD
      write (*, '(4x,a4,a27,a,i6,a)') 'TIM ', trsi(i), ' Epo ', i, ' deleted as huge residual'
    endif
  enddo

  return
end
