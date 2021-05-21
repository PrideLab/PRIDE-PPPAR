!
!! remov_shrt.f90
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
!! purpose  : remove short arc
!! parameter:
!!    input : shrt -- short limit
!!            nepo -- # of epochs of this residual series
!!            flag -- flags of each epoch
!!            trsi -- time tags for each epoch
!!    output:
!
subroutine remov_shrt(lupd, shrt, nepo, flag, trsi)
  implicit none
  include '../header/const.h'
  include 'data_flag.h'

  logical*1 lupd
  integer*4 nepo, shrt, flag(1:*)
  character*27 trsi(1:*)
!
!! local
  logical*1 lfnd
  integer*4 i, j, k, iepo, nok
!
!! function called
  logical*1 istrue

  k = 1
  lfnd = .true.
  do while (lfnd)
    call find_flag(k, nepo, flag, 'OLDAMB', i)
    if (i .gt. 0) then
      if (i .eq. nepo) then
        lupd = .true.
        flag(i) = DELSHT
        write (*, '(4x,a4,a27,a,i6,a)') 'TIM ', trsi(i), ' Epo ', i, ' removed as short'
        exit
      else
        call find_flag(i + 1, nepo, flag, 'OLDAMB', k)
        if (k .lt. 0) k = nepo + 1
        call find_flag(k - 1, i, flag, 'OK', j)
        nok = 0
        do iepo = i, j
          if (istrue(flag(iepo), 'OK')) nok = nok + 1
        enddo
        if (nok .le. shrt) then
          do iepo = i, j
            lupd = .true.
            if (istrue(flag(iepo), 'OK')) then
              flag(iepo) = DELSHT
              write (*, '(4x,a4,a27,a,i6,a)') 'TIM ', trsi(iepo), ' Epo ', iepo, ' removed as short'
            endif
          enddo
        endif
      endif
    else
      lfnd = .false.
    endif
    if (k .eq. nepo + 1) exit
  enddo

  return
end
