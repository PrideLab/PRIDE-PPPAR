!
!! check_arcslip.f90
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
subroutine check_arcslip(jump, nepo, flag, resi, trsi, lfnd)
  implicit none
  include '../header/const.h'
  include 'data_flag.h'

  logical*1 lfnd
  integer*4 nepo, flag(1:*)
  real*8 jump, resi(1:*)
  character*27 trsi(1:*)
!
!! local
  integer*4 i, j, k, frt, ndif, imax, ipt(nepo), iepo
  integer*4 flg(nepo), ifirst, ilast, naf, nbe
  real*8 mean, sig, signew, jmp, difmax, dif(nepo), maf, mbe
!
!! function called
  logical*1 istrue, chitst

  lfnd = .false.
  do i = 1, nepo
    flg(i) = 2
    if (istrue(flag(i), 'AMB')) flg(i) = 1
    if (istrue(flag(i), 'GOOD')) flg(i) = 0
  enddo

  iepo = 1
  do while(iepo .le. nepo)
    if (flg(iepo) .ne. 1) then
      iepo = iepo + 1
      cycle
    endif
    10 continue
    ifirst = iepo
    ilast = ifirst
    do j = ifirst+1, nepo
      if(flg(j) .eq. 1) exit
      if(flg(j) .eq. 2) cycle
      ilast = j
      if(j .eq. nepo) ilast = j
    enddo
    if(ilast .eq. ifirst) goto 50
    !
    !! form difference between adjacent epochs
    do i = ilast, ifirst+1, -1
      dif(i) = 0.d0
      ipt(i) = 0
      if(flg(i) .ge. 1) cycle
      do j = i-1, ifirst, -1
        if(flg(j) .eq. 2) cycle
        dif(i) = resi(i) - resi(j)
        ipt(i) = j
        exit
      enddo
    enddo
    !!
    !!! mean of residual time series in one arc
    ndif = 0
    mean = 0.d0
    sig = 0.d0
    imax = 0
    difmax = -1.d0
    do i = ilast, ifirst+1, -1
      if(flg(i) .eq. 0 .and. ipt(i) .ne. 0) then
        ndif = ndif + 1
        mean = mean + dif(i)
        sig = sig + dif(i)**2
        if (dabs(dif(i)) .gt. difmax) then
          imax = i
          difmax = dif(i)
        endif
      endif
    enddo
    if (ndif .le. 1) goto 50
    signew = sqrt((sig-difmax**2)/(ndif-1) - ((mean-difmax)/(ndif-1))**2)
    mean = mean/ndif
    sig = sqrt(sig/ndif - mean**2)
    if (.not. chitst(-1, ndif - 1, sig, signew, 0.99d0) .and. dabs(dif(imax)) .gt. jump) then
      write (*, '(4x,a4,a27,a,i6,a)') 'TIM ', trsi(imax), ' Epo ', imax, ' marked as biggest jump'
      flag(i) = NEWAMB
      flg(imax) = 1
      goto 10
    endif
    !
    !! new jump setting
    jmp = jump/2
    if(jmp .lt. 3*sig) jmp = 3*sig
    !
    !! check jump
    do i = ilast, ifirst+1, -1
      if (dabs(dif(i)) .gt. jmp .and. flg(i) .eq. 0) then
        nbe = 0
        naf = 0
        mbe = 0.d0
        maf = 0.d0
        do j = ifirst, i
          if(flg(j) .eq.2 .or. j .eq. i) cycle
          mbe = mbe + resi(j)
          nbe = nbe + 1
          if(nbe .eq. 50) exit
        enddo
        do j = i, ilast
          if (flg(j) .eq. 2) cycle
          maf = maf + resi(j)
          naf = naf + 1
          if(naf .eq. 50) exit
        enddo
        if (nbe .eq. 0 .or. naf .eq. 0) cycle
        mbe = mbe/nbe
        maf = maf/naf
        if (dabs(maf-mbe) .gt. jmp/2) then
          write (*, '(4x,a4,a27,a,i6,a)') 'TIM ', trsi(i), ' Epo ', i, ' marked as jump'
          flg(i) = 1
          flag(i) = NEWAMB
          goto 10
        endif
      endif
    enddo
    !
    !! find bad
    do i = ilast, ifirst+1, -1
      if (ipt(i) .eq. 0 .or. ipt(ipt(i)) .eq. 0) cycle
      if (istrue(flag(i), 'NEWAMB') .and. istrue(flag(ipt(i)), 'NEWAMB')) then
        if (dsign(1.d0, dif(i))*dsign(1.d0, dif(ipt(i))) .lt. 0.d0) then
          write (*, '(4x,a4,a27,a,i6,a)') 'TIM ', trsi(ipt(i)), ' Epo ', ipt(i), ' deleted as bad'
          flg(ipt(i)) = 2
          flag(ipt(i)) = DELBAD
          dif(i) = dif(i) + dif(ipt(i))
          k = ipt(i)
          ipt(i) = ipt(ipt(i))
          dif(k) = 0.d0
          ipt(k) = 0
          if (dabs(dif(i)) .lt. jmp) then
            write (*, '(4x,a4,a27,a,i6,a)') 'TIM ', trsi(i), ' Epo ', i, ' reset as good'
            flg(i) = 0
            flag(i) = GOOD
          endif
        endif
      endif
    enddo
    50 continue
    iepo = ilast+1
  enddo

  return
end
