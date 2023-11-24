!
!! stpir.f90
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
subroutine mstpir(nepo, flagall, ti, pg)
  implicit none
  integer*4 nepo, flagall(1:*)
  real*8 ti(1:*), pg(1:*)

  logical*1 istrue
  integer*4 set_flag, ifirst, jlast, jnext, isecond, j0, j1
  real*8 ddpg(nepo, 2), ddmpg, ddsigb, dpg
  integer*4 iepo, flg(nepo), jepo, kobs, nwin, wine, step, iws, iwe

  do iepo = 1, nepo
    flg(iepo) = 2
    ddpg(iepo, 1:2) = 999
    if (istrue(flagall(iepo), 'ok')) then
      flg(iepo) = 1
      if (.not. istrue(flagall(iepo), 'lgjump') .and. &
          .not. istrue(flagall(iepo), 'lwjump') .and. &
          .not. istrue(flagall(iepo), 'gap')) flg(iepo) = 0
      !if (.not. istrue(flagall(iepo), 'lgjump') .and. &
      !    .not. istrue(flagall(iepo), 'gap')) flg(iepo) = 0
    endif
  enddo
  
  iepo = 1
  nwin = 50
  step = 20
  do while(iepo .le. nepo)
    if (flg(iepo) .ne. 1) then
      iepo = iepo + 1
      cycle
    endif
    ifirst = iepo
    10 continue
    isecond = ifirst
    ddmpg = 0.d0
    ddsigb = 0.d0
    kobs = 0
    do jepo = ifirst+1, nepo
      if(flg(jepo) .eq. 1) exit
      if(flg(jepo) .eq. 2) cycle
      isecond = jepo
      if(jepo .eq. nepo) isecond=nepo
    enddo
    if (isecond-ifirst .lt. 3) goto 50
    do jepo = ifirst+1, isecond-1
      ddpg(jepo, 1:2) = 999
      jlast = jepo
      jnext = jepo
      if(flg(jepo) .eq. 2) cycle
      do j0 = jepo-1, ifirst, -1
        if (flg(j0) .lt. 2) then
          jlast = j0
          exit
        endif
        if(flg(j0) .eq. 2) jlast = jlast - 1
      enddo
      if (jlast .lt. ifirst) cycle
      do j1 = jepo+1, isecond
        if (flg(j1) .eq. 0) then
          jnext = j1
          exit
        endif
        if(flg(j1) .eq. 2) jnext = jnext + 1
      enddo
      if (jnext .gt. isecond) cycle
      ddpg(jepo, 1) = pg(jnext)-2*pg(jepo)+pg(jlast)
      ddpg(jepo, 2) = jlast
    enddo
    if ((isecond-ifirst) .le. nwin) then
      kobs = 0
      ddmpg = 0.d0
      ddsigb = 0.d0
      do j0 = ifirst, isecond
          ddmpg = ddmpg + ddpg(j0, 1)
          ddsigb = ddsigb + ddpg(j0, 1)*ddpg(j0, 1)
          kobs = kobs + 1
      enddo
      if (kobs .eq. 0) goto 50
      ddmpg = ddmpg/kobs
      ddsigb = ddsigb/kobs - ddmpg*ddmpg
      ddsigb = sqrt(ddsigb)
      do jepo = ifirst+1, isecond-1
        if(ddpg(jepo, 1) .eq. 999) cycle
        if(dabs(ddpg(jepo, 1)-ddmpg) .gt. 4*ddsigb) then
          jlast = ddpg(jepo, 2)
          if (flg(jlast) .eq. 3) then
            flg(jepo) = 1
            goto 10
          else if (flg(jlast) .ne. 3) then
            flg(jepo) = 3
          endif
        endif
      enddo
    else
      jepo = ifirst
      do while(jepo .lt. isecond-nwin)
        iws = jepo
        iwe = jepo+nwin
        if (iwe .gt. isecond) iwe = isecond
        kobs = 0
        ddmpg = 0.d0
        ddsigb = 0.d0
        do j0 = iws, iwe
          if(ddpg(j0, 1) .eq. 999) cycle
          ddmpg = ddmpg + ddpg(j0, 1)
          ddsigb = ddsigb + ddpg(j0, 1)*ddpg(j0, 1)
          kobs = kobs + 1
        enddo
        if (kobs .eq. 0) goto 20
        ddmpg = ddmpg/kobs
        ddsigb = ddsigb/kobs - ddmpg*ddmpg
        ddsigb = sqrt(ddsigb)
        do j0 = iws, iwe
          if(ddpg(j0, 1) .eq. 999) cycle
          if(dabs(ddpg(j0, 1)-ddmpg) .gt. 4*ddsigb) then
            jlast = ddpg(j0, 2)
            if (flg(jlast) .eq. 3) then
              flg(j0) = 1
              if (j0 .ne. ifirst) goto 10
            else if (flg(jlast) .ne. 3) then
              flg(j0) = 3
            endif
          endif
        enddo
        20 continue
        jepo = jepo+step
      enddo
    endif
    50 continue
    iepo = isecond + 1
  enddo

  do iepo = 1, nepo
    if (flg(iepo) .eq. 1 .and. .not. istrue(flagall(iepo), 'lgjump')) &
      flagall(iepo) = set_flag(flagall(iepo), 'lgjump')
  enddo

  return
end
