!
!! check_slip.f90
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
!! Contributor: Maorong Ge, Jianghui Geng, Songfeng Yang
!! 
!!
!!
!! purpose  : check residual time series to identify slips or blunders
!! parameter:
!!    input : 
!!            jump -- residual jump threshold value
!!            nepo -- # of residuals
!!            flag -- flags for each residual
!!            resi -- residuals
!!            trsi -- time tag for each residual
!!    output: lfnd -- jump found or not
!
subroutine check_slip(jump, nepo, flag, resi, trsi, lfnd)
  implicit none
  include '../header/const.h'
  include 'data_flag.h'

  logical*1 lfnd
  integer*4 nepo, flag(1:*)
  real*8 jump, resi(1:*)
  character*27 trsi(1:*)
!
!! local
  integer*4 i, j, k, frt, ndif, imax, ipt(nepo)
  real*8 jmp, difmax, mean, sig, signew, dif(nepo)
!
!! function called
  logical*1 istrue, chitst

  lfnd = .false.
!
!! form difference between adjacent epochs
  dif(1) = 0.d0
  ipt(1) = 0
  do i = nepo, 2, -1
    dif(i) = 0.d0
    ipt(i) = 0
    if (.not. istrue(flag(i), 'OK')) cycle
    do j = i - 1, 1, -1
      if (.not. istrue(flag(j), 'OK')) cycle
      dif(i) = resi(i) - resi(j)
      ipt(i) = j
      exit
    enddo
  enddo
!
!! mean of residual time series
100 continue
  ndif = 0
  mean = 0.d0
  imax = 0
  difmax = -1.d0
  do i = nepo, 1, -1
    if (ipt(i) .ne. 0) frt = i
    if (istrue(flag(i), 'GOOD')) then
      ndif = ndif + 1
      mean = mean + dif(i)
      if (dabs(dif(i)) .gt. difmax) then
        imax = i
        difmax = abs(dif(i))
      endif
    endif
  enddo
  if (ndif .eq. 0) return
!
!! sigma of residual time series
  mean = mean/ndif
  sig = 0.d0
  do i = frt, nepo
    if (ipt(i) .ne. 0 .and. istrue(flag(i), 'GOOD')) then
      sig = sig + (dif(i) - mean)**2
    endif
  enddo
  sig = dsqrt(sig/ndif)
!
!! new sigma when max residual difference is removed
  mean = (mean*ndif - dif(imax))/(ndif - 1)
  signew = 0.d0
  do i = frt, nepo
    if (imax .eq. i) cycle
    if (ipt(i) .ne. 0 .and. istrue(flag(i), 'GOOD')) then
      signew = signew + (dif(i) - mean)**2
    endif
  enddo
  signew = dsqrt(signew/(ndif - 1))
!
!! whether to remove the largest residual difference
  if (.not. chitst(-1, ndif - 1, sig, signew, 0.99d0) .and. dabs(dif(imax)) .gt. jump) then
    write (*, '(4x,a4,a27,a,i6,a)') 'TIM ', trsi(imax), ' Epo ', imax, ' marked as biggest jump'
    flag(imax) = NEWAMB
    lfnd = .true.
    goto 100
  endif
!
!! new jump setting
  jmp = jump
  if (jmp .lt. 3.d0*sig) jmp = 3.d0*sig
!
!! check jump
  do i = frt, nepo
    if (dabs(dif(i)) .gt. jmp) then
      if (istrue(flag(i), 'GOOD')) then
        write (*, '(4x,a4,a27,a,i6,a)') 'TIM ', trsi(i), ' Epo ', i, ' marked as jump'
        flag(i) = NEWAMB
        lfnd = .true.
      endif
    else
      if (istrue(flag(i), 'NEWAMB')) then
        write (*, '(4x,a4,a27,a,i6,a)') 'TIM ', trsi(i), ' Epo ', i, ' reset as good'
        flag(i) = GOOD
      endif
    endif
  enddo
!
!! find bad
  do i = nepo, frt, -1
    if (ipt(i) .eq. 0 .or. ipt(ipt(i)) .eq. 0) cycle
    if (istrue(flag(i), 'NEWAMB') .and. istrue(flag(ipt(i)), 'NEWAMB')) then
      if (dsign(1.d0, dif(i))*dsign(1.d0, dif(ipt(i))) .lt. 0.d0) then
        write (*, '(4x,a4,a27,a,i6,a)') 'TIM ', trsi(ipt(i)), ' Epo ', ipt(i), ' deleted as bad'
        flag(ipt(i)) = DELBAD
        dif(i) = dif(i) + dif(ipt(i))
        k = ipt(i)
        ipt(i) = ipt(ipt(i))
        dif(k) = 0.d0
        ipt(k) = 0
        if (dabs(dif(i)) .lt. jmp) then
          write (*, '(4x,a4,a27,a,i6,a)') 'TIM ', trsi(i), ' Epo ', i, ' reset as good'
          flag(i) = GOOD
        endif
      endif
    endif
  enddo
!
!! handling special conditions
  do i = frt, nepo
    if (istrue(flag(i), 'AMB')) then
      k = ipt(i)
      j = i + 1
      do while (ipt(j) .eq. 0 .and. j .le. nepo)
        j = j + 1
      enddo
      if (k .lt. frt) then
        write (*, '(4x,a4,a27,a,i6,a)') 'TIM ', trsi(k), ' Epo ', k, ' deleted as redundant ambiguity'
        flag(k) = DELBAD
      else if (j .gt. nepo) then
        write (*, '(4x,a4,a27,a,i6,a)') 'TIM ', trsi(i), ' Epo ', i, ' deleted as last ambiguity'
        flag(i) = DELBAD
      else if (istrue(flag(j), 'AMB')) then
        if (istrue(flag(i), 'OLDAMB')) then
          write (*, '(4x,a4,a27,a,i6,a)') 'TIM ', trsi(i), ' Epo ', i, ' deleted as redundant ambiguity'
          flag(i) = DELBAD
          flag(j) = OLDAMB
        else if (istrue(flag(i), 'NEWAMB')) then
          write (*, '(4x,a4,a27,a,i6,a)') 'TIM ', trsi(i), ' Epo ', i, ' deleted as redundant ambiguity'
          flag(i) = DELBAD
        endif
      endif
    endif
  enddo

  return
end
