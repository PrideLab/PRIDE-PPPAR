!
!! check_pcres.f90
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
!! purpose  : check pc residual time series to identify blunders
!! parameter:
!!    input : 
!!            jump -- residual jump threshold value
!!            nepo -- # of residuals
!!            flag -- flags for each residual
!!            resi -- residuals
!!            trsi -- time tag for each residual
!!    output: lfnd -- jump found or not
!
subroutine check_pcres(jump, nepo, flag, resp, trsi) !, lfnd)
  implicit none
  include '../header/const.h'
  include 'data_flag.h'

  logical*1 lfnd
  integer*4 nepo, flag(1:*)
  real*8 jump, resp(1:*)
  character*27 trsi(1:*)
!
!! local
  integer*4 i, j, k, nres, imax
  real*8 jmp, resmax, mean, sig, signew
!
!! function called
  logical*1 istrue, chitst

  lfnd = .false.
!
!! mean of PC residual time series
100 continue
  nres = 0
  mean = 0.d0
  imax = 0
  resmax = -1.d0
  do i = nepo, 1, -1
    if (istrue(flag(i), 'GOOD')) then
      nres = nres + 1
      mean = mean + resp(i)
      if (dabs(resp(i)) .gt. resmax) then
        imax = i
        resmax = abs(resp(i))
      endif
    endif
  enddo
  if (nres .eq. 0) return
!
!! sigma of residual time series
  mean = mean/nres
  sig = 0.d0
  do i = 1, nepo
    if (istrue(flag(i), 'GOOD')) then
      sig = sig + (resp(i) - mean)**2
    endif
  enddo
  sig = dsqrt(sig/nres)
!
!! new sigma when max residual difference is removed
  mean = (mean*nres - resp(imax))/(nres - 1)
  signew = 0.d0
  do i = 1, nepo
    if (imax .eq. i) cycle
    if (istrue(flag(i), 'GOOD')) then
      signew = signew + (resp(i) - mean)**2
    endif
  enddo
  signew = dsqrt(signew/(nres - 1))
 !write (*,*) "mean, sig, signew, max:",mean,sig,signew,resp(imax)
!
!! whether to remove the largest residual difference
  if (.not. chitst(-1, nres - 1, sig, signew, 0.99d0) .and. dabs(resp(imax)) .gt. jump) then
    !write (*, '(4x,a4,a27,a,i6,a)') 'TIM ', trsi(imax), ' Epo ', imax, ' marked as biggest PC outlier'
    flag(imax) = DELBAD
    lfnd = .true.
    goto 100
  endif
!
!! new jump setting
  !write(*,*) mean,sig,signew,resp(imax),jump,3.d0*sig
  jmp = jump
  !if (jmp .lt. 3.d0*sig) jmp = 3.d0*sig
!
!! check jump
  do i = 1, nepo
    if (dabs(resp(i)) .gt. jmp) then
      if (istrue(flag(i), 'GOOD')) then
        write (*, '(4x,a4,a27,a,i6,a)') 'TIM ', trsi(i), ' Epo ', i, ' marked as outlier'
        flag(i) = DELBAD
        lfnd = .true.
      endif
    else
      if (istrue(flag(i), 'DELBAD')) then
        write (*, '(4x,a4,a27,a,i6,a)') 'TIM ', trsi(i), ' Epo ', i, ' reset as good'
        flag(i) = GOOD
      endif
    endif
  enddo
  return
end
