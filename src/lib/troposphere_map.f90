!
!! troposphere_map.f90
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
!! Contributor: Maorong Ge, Jianghui Geng, Songfeng Yang, Jing Zeng
!! 
!!
!!
!! purpose  : map functions of troposhperic delay correction. Here we directly call several
!!            ztd model and mapping function code from GAMIT
!! parameter: jd -- julian day (GPST)
!!            elev -- elevation in radian
!!            SITE%-- station information struct
!!            damp,wmap -- map functions of hydrostatic and wet troposphere delay
!
subroutine troposphere_map(jd, sod, elev, SITE, dmap, wmap)
  implicit none
  include '../header/const.h'
  include '../header/station.h'

  integer*4 jd
  real*8 sod, elev, dmap, wmap
  type(station) SITE
!
!! local
  integer*4 i, iyr, idoy, jdutc
  real*8 radian, wmf(2), hmf(2), tlen, al0, al1, ah, aw, sodutc
!
!! function called
  real*8 cfa, timdif, taiutc

  radian = datan(1.d0)/45.d0 ! radian per degree
  call timinc(jd, sod, 19.d0 - taiutc(jd), jdutc, sodutc)
!
!! mapping functions
  if (SITE%map(1:3) .eq. 'VM1' .or. SITE%map(1:3) .eq. 'VM3') then           ! vienna mapping function 1 or 3
    do i = 1, SITE%nvm - 1
      if (timdif(SITE%jdv(i), SITE%sodv(i), jdutc, sodutc) .le. 0.d0 .and. &
          timdif(SITE%jdv(i + 1), SITE%sodv(i + 1), jdutc, sodutc) .ge. 0.d0) then
!! linear interpolation
        tlen = timdif(SITE%jdv(i + 1), SITE%sodv(i + 1), SITE%jdv(i), SITE%sodv(i))
        al0 = timdif(SITE%jdv(i + 1), SITE%sodv(i + 1), jdutc, sodutc)/tlen
        al1 = timdif(jdutc, sodutc, SITE%jdv(i), SITE%sodv(i))/tlen
        ah = al0*SITE%vm1(1, i) + al1*SITE%vm1(1, i + 1)
        aw = al0*SITE%vm1(2, i) + al1*SITE%vm1(2, i + 1)
!
!!
        if(SITE%map(1:3).eq.'VM1') then           ! vienna mapping function 1
          call vmf1_ht(ah, aw, jdutc + sodutc/86400.d0, SITE%geod(1), SITE%geod(3)*1.d3, PI/2.d0 - elev, dmap, wmap)
        else if(SITE%map(1:3).eq.'VM3') then      ! vienna mapping function 3
          call vmf3_ht(ah, aw, jdutc + sodutc/86400.d0, SITE%geod(1), SITE%geod(2), SITE%geod(3)*1.d3, PI/2.d0-elev, dmap, wmap)
        endif
        exit
      endif
    enddo
  else if (SITE%map(1:3) .eq. 'GMF') then      ! global mapping function
    call global_map(jdutc*1.d0, SITE%geod(1), SITE%geod(2), SITE%geod(3)*1.d3, PI/2.d0 - elev, dmap, wmap)
  else if (SITE%map(1:3) .eq. 'NIE') then      ! Niell mapping function
    call mjd2doy(jdutc, iyr, idoy)
    call nmfh2p1(idoy*1.d0, SITE%geod(1)/radian, SITE%geod(3)*1.d3 - SITE%undu, elev/radian, hmf)
    dmap = hmf(1)
    call nmfw2(SITE%geod(1)/radian, elev/radian, wmf)
    wmap = wmf(1)
  else if (SITE%map(1:3) .eq. 'CFA') then      ! Chao mapping function
    dmap = cfa(SITE%p0, SITE%t0, SITE%hr0, 'R', SITE%geod(1), SITE%geod(3), elev)
    wmap = dmap
  else if (SITE%map(1:3) .eq. 'NON') then
    dmap = 0.d0
    wmap = 0.d0
  else
    write (*, '(2a)') '***ERROR(troposphere_map): unknow mapping function ', SITE%map(1:3)
    call exit(1)
  endif

  return
end
