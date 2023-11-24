!
!! lfbmwa.f90
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
!! Contributor: Jing Zeng
!!
!
subroutine lfbmwa(nepo, ti, nw, elev, flagall, limit)
  implicit none
  integer*4 nepo, flagall(1:*)
  real*8 ti(1:*), nw(1:*), elev(1:*), limit

! function called
  logical*1 ok, istrue
  integer*4 set_flag

  integer*4 iepo, i, j, flg(nepo), ifirst, ilast, imaxdiff, iprevious
  integer*4 nba, nfo, numba, numfo, ifoe, ibas, jfoe, jbas
  real*8 mfo(nepo), mba(nepo), sigfo(nepo), sigba(nepo), maxdiff, mbe, maf, weight, kobs
  real*8 vmaxba, vmaxfo, coeffo(2), coefba(2), rmsfo, rmsba, mnw, signw
  logical*1 flgba, flgfo

! simply check for jump
  do i = 1, nepo
    flg(i) = 2
    if (istrue(flagall(i), 'ok')) then
      flg(i) = 1
      if (.not. istrue(flagall(i), 'lwjump') .and. &
          .not. istrue(flagall(i), 'gap') .and. &
          .not. istrue(flagall(i), 'bigsd')) flg(i) = 0
    endif
  enddo

! forward and backward
  nfo = 30
  nba = 30
  mfo = 0.d0
  mba = 0.d0
  sigfo = 0.d0
  sigba = 0.d0
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
      if(j .eq. nepo) ilast=j
    enddo
    if (ilast-ifirst .lt. nfo) goto 50
    numba = 0
    numfo = 0
    ibas = -1
    ifoe = -1
    do i = ifirst, ilast
      if (flg(i) .ne. 2) numba = numba+1
      if (numba .eq. nba) ibas = i ! backward start point
      j = ilast-i+ifirst
      if (flg(j) .ne. 2) numfo = numfo+1
      if (numfo .eq. nfo) ifoe = j ! forward end point
    enddo
    if (ifoe .eq. -1) goto 50
    ! forward and turboedit
    do i = ifirst, ifoe
      if (flg(i) .eq. 2) cycle
      kobs = 0.d0
      numfo = 0
      mfo(i) = 0.d0
      sigfo(i) = 0.d0
      do j = i, ilast
        if (flg(j) .eq. 2) cycle
        weight = 1.d0
        if (elev(j) .lt. 10.d0) weight = 0.3d0
        mfo(i) = mfo(i)+nw(j)*weight
        sigfo(i) = sigfo(i)+nw(j)**2*weight
        kobs = kobs + weight
        numfo = numfo + 1
        if (numfo .eq. nfo) exit
      enddo
      if (kobs .eq. 0.d0) cycle
      mfo(i) = mfo(i)/kobs
      sigfo(i) = sqrt(sigfo(i)/kobs-mfo(i)**2)
    enddo
    ! reverse for forward turboedit
    numfo = 0
    do i = ifoe, ifirst, -1
      if (flg(i) .eq. 2) cycle
      numfo = numfo+1
      if (numfo .le. 15) cycle ! reduce influence of gross error
      iprevious = -1
      do j = i+1, ilast
        if(flg(j) .eq. 0) then
          iprevious = j
          exit
        endif
      enddo
      if (flg(i) .eq. 2) cycle
      if (dabs(nw(i)-mfo(i)) .gt. 3*limit .or. iprevious .ne. -1 .and. &
          dabs(nw(i)-mfo(i)) .gt. 4.2*limit .and. dabs(nw(i)-nw(iprevious)) &
          .gt. 0.5) then
        flg(iprevious) = 1    ! reverse!!!
        goto 10
      endif
    enddo
    if (ibas .eq. -1) goto 50
    ! backward and bfmw
    maxdiff = -1
    imaxdiff = -1
    numfo = 0    ! number conter for backward turboedit
    do i = ibas, ilast
      if (flg(i) .eq. 2) cycle
      kobs = 0.d0
      numba = 0
      mba(i) = 0.d0
      sigba(i) = 0.d0
      do j = i, ifirst, -1
        if (flg(j) .eq. 2) cycle
        weight = 1.d0
        if (elev(j) .lt. 10.d0) weight = 0.3d0
        mba(i) = mba(i) + nw(j)*weight
        sigba(i) = sigba(i) + nw(j)**2*weight
        kobs = kobs + weight
        numba = numba + 1
        if (numba .eq. nba) exit
      enddo
      if (kobs .eq. 0.d0) cycle
      mba(i) = mba(i)/kobs
      sigba(i) = sqrt(sigba(i)/kobs - mba(i)**2)
      if (mba(i)*mfo(i) .ne. 0.d0 .and. i .le. ifoe .and. &
          dabs(mba(i)-mfo(i)) .ge. maxdiff .and. flg(i) .ne. 2) then
        maxdiff = dabs(mba(i)-mfo(i))
        imaxdiff = i
      endif
      iprevious = -1 ! backward turboedit
      do j = i-1, ifirst, -1
        if (flg(j) .eq. 0) then
          iprevious = j
          exit
        endif
      enddo
      numfo = numfo + 1
      if (numfo .le. 15) cycle
      if (dabs(nw(i)-mba(i)) .gt. 3*limit .or. iprevious .ne. -1 .and. &
          dabs(nw(i)-mba(i)) .gt. 4.2*limit .and. dabs(nw(i)-nw(iprevious)) &
          .gt. 0.5) then
        flg(i) = 1
        goto 10
      endif
    enddo
    flgfo = .false.
    flgba = .false.
    if (imaxdiff .ne. -1) then
      call linearfit(ti(imaxdiff), mba(imaxdiff), flg(imaxdiff), 1, ilast-imaxdiff, nba, coefba, rmsba, flgba)
      call linearfit(ti(ifirst), mfo(ifirst), flg(ifirst), -1, imaxdiff-ifirst+1, nfo, coeffo, rmsfo, flgfo)
      if (flgfo .or. flgba) then
        flg(imaxdiff) = 1
        goto 10
      endif
    endif
    50 continue
    iepo = ilast + 1
  enddo
! falg cycle silp of MW
  do i = 1, nepo
    if (flg(i) .eq. 1 .and. .not. istrue(flagall(i), 'lwjump')) &
      flagall(i) = set_flag(flagall(i), 'lwjump')
  enddo
  return
end

subroutine linearfit(x, y, flag, index, endp, ndgr, coeff, rms, flgfit)  ! y = kx+b
  implicit none
  integer*4 i, ndgr, flag(1:*), index, endp, nump, ifirst, ilast
  real*8 dot, x(1:*), y(1:*), vmax, coeff(2), v, rms
  real*8 sumx, sumy, sumxx, sumxy, yfit
  logical*1 flgfit

  ifirst = -1
  ilast = -1
  sumx = 0.d0
  sumy = 0.d0
  sumxx = 0.d0
  sumxy = 0.d0
  flgfit = .false.
  if (x(endp) .le. x(1)) return
  nump = 0
  if (index .eq. 1) then
    do i=1,endp
      if (flag(i) .eq. 2) cycle
      if (ifirst .eq. -1) ifirst = i
      sumx = sumx + x(i)
      sumy = sumy + y(i)
      sumxx = sumxx + x(i)**2
      sumxy = sumxy + x(i)*y(i)
      nump = nump+1
      ilast = i
      if (nump .eq. ndgr) exit
    enddo
  else if (index .eq. -1) then
    do i=endp, 1, -1
      if (flag(i) .eq. 2) cycle
      if (ilast .eq. -1) ilast = i
      sumx = sumx + x(i)
      sumy = sumy + y(i)
      sumxx = sumxx + x(i)**2
      sumxy = sumxy + x(i)*y(i)
      nump = nump+1
      ifirst = i
      if (nump .eq. ndgr) exit
    enddo
  endif
  if (ilast .eq. -1 .or. ifirst .eq. -1) return
  coeff(1) = (sumy*sumxx-sumx*sumxy)/(nump*sumxx-(sumx)**2)   ! b
  coeff(2) = (nump*sumxy-sumx*sumy)/(nump*sumxx-(sumx)**2)    ! k
  vmax = -1.d0
  rms = 0.d0
  nump = 0
  if (index .eq. 1) then
    do i = 1, endp
      if (flag(i) .eq. 2) cycle
      yfit = coeff(1)+coeff(2)*x(i)
      v = dabs(yfit-y(i))
      if (v .ge. vmax) vmax=v
      rms = rms+v**2
      nump = nump+1
      if (nump .eq. ndgr) exit
    enddo
  else if (index .eq. -1) then
    do i = endp, 1, -1
      if (flag(i) .eq. 2) cycle
      yfit = coeff(1)+coeff(2)*x(i)
      v = yfit-y(i)
      if (v .ge. vmax) vmax=v
      rms = rms+v**2
      nump = nump+1
      if (nump .eq. ndgr) exit
    enddo
  endif
  rms = sqrt(rms/nump)
  if(rms .lt. 2.d-2 .and. dabs(coeff(2)) .gt. 0.5d0/(x(ilast)-x(ifirst))) flgfit=.true. !0.5cycle
  return
end
