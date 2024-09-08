!
!! elevation.f90
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
!! Contributor: Maorong Ge, Jianghui Geng, Songfeng Yang, Jihang Lin
!!
!!
!!
!! purpose  : compute elevation angle of satellie
!! parameter:
!!   input  : neph,ephem ----- ephemeris number and ephemeris
!!            iprn ----------- satellite number
!!            jd0 ti --------- epoch time (julday)
!!            x,y,z ---------- receiver position
!!            elev ----------- the elevation angle
!!            dist ----------- the distance between satellite and reciver
!!            dtsat ---------- clock correction, units: meter
!!            v -------------- true anomaly
!!            jeph ----------- the prn of ephemris which is chosen
!!            lstick --------- use the last one (jeph) if it is not zero, to avoid discontinuity by two brdeph
!!   output :
!!
!
subroutine elevation(neph, ephem, iprn, jd0, ti, x, y, z, elev, dist, dtsat, v, jeph, lstick)
  implicit none
  include '../header/const.h'
  include '../header/brdeph.h'

  type(brdeph) ephem(1:*)
  integer*4 neph, nepo, jeph, jd0
  character*3 iprn
  real*8 ti, x, y, z, elev, v
  logical*1 lstick
  integer*4 ite

! local
  real*8 dtsat, xsat(6), tsend, tsend1, dist, r12, r1, cosz, dtref

! function called
  integer*4 modified_julday
  real*8 dot

  if (x*y*z .eq. 0.d0) then
    dist = -1.d0
    return
  end if

  tsend = ti - 0.075d0
  ite = 0
  if (index("GECJ", iprn(1:1)) .ne. 0) then
    do while (.true.)
!
!!     xsat, dtsat ---- position , velocity and clock correction
      call brdxyz('yyy', lstick, neph, ephem, iprn, jeph, jd0, tsend, xsat, v, dtsat, dtref)
      dtsat = dtsat - 2.d0*dot(3, xsat, xsat(4))/vlight**2 ! relativity effects
!
!! not found
      if (dabs(dtref) .gt. 2.d0) then
!       write(*,'(a,a3,f5.1)') '###WARNING(elevation): no ephemeris for PRN, &
!                              dt(hours)',iprn,dtref
        dist = -1.d0
        return
      end if
!
!! bad orbit
      if (any(dabs(xsat(1:3)) .gt. 1.d9)) then
        dist = -1.d0
        return
      end if

      dist = dsqrt((xsat(1) - x)**2 + (xsat(2) - y)**2 + (xsat(3) - z)**2)
      tsend1 = ti - dist/VLIGHT
      ite = ite + 1
      if (dabs(tsend - tsend1) .lt. 1.d-11 .or. ite .ge. 100) exit
      tsend = tsend1
    end do
!
! zenith-distance (approx.)
    r12 = x*(xsat(1) - x) + y*(xsat(2) - y) + z*(xsat(3) - z)
    r1 = dsqrt(x*x + y*y + z*z)
    cosz = r12/r1/dist
    elev = 90.d0 - dacos(cosz)/PI*180.d0
    dtsat = dtsat*VLIGHT

    return
  elseif (iprn(1:1) .eq. "R") then
    do while (.true.)
!
!! xsat, dtsat ---- position , velocity and clock correction
      call glsbrd('yyy', neph, ephem, iprn, jeph, jd0, tsend, xsat, v, dtsat, dtref)
      dtsat = dtsat - 2.d0*dot(3, xsat, xsat(4))/vlight**2      ! relativity effects
!
!! not found
      if (dabs(dtref) .gt. 0.5d0) then
!        write (*,'(a,i3,f5.1)') '###WARNING(elevation): no ephemeris for PRN, &
!                                dt(hours)',iprn,dtref
        dist = -1.d0
        return
      end if
!
!! bad orbit
      if (any(dabs(xsat(1:3)) .gt. 1.d9)) then
        dist = -1.d0
        return
      end if

      dist = dsqrt((xsat(1) - x)**2 + (xsat(2) - y)**2 + (xsat(3) - z)**2)
      tsend1 = ti - dist/VLIGHT
      if (dabs(tsend - tsend1) .lt. 1.d-9) exit
      tsend = tsend1
    end do
    !
    ! zenith-distance (approx.)
    r12 = x*(xsat(1) - x) + y*(xsat(2) - y) + z*(xsat(3) - z)
    r1 = dsqrt(x*x + y*y + z*z)
    cosz = r12/r1/dist
    elev = 90.d0 - dacos(cosz)/PI*180.d0
    dtsat = dtsat*VLIGHT

    return
  end if
end subroutine
