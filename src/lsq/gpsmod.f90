!
!! gpsmod.f90
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
!! Contributor: Maorong Ge, Jianghui Geng, Songfeng Yang, Jihang Lin, Jing Zeng
!!
!!
!!
!! purpose   : modelling of GNSS observations, observation equaitons are storaged in
!!             OB%amat and OB%omc etc.
!! parameters:
!!             jd,sod -- epoch time of the observations (GPS time)
!!             LCF    -- LSQ control parameters
!!             SITE   -- SITE information
!!             OB     -- observations
!!             SAT    -- satellite information
!!
!
subroutine gpsmod(jd, sod, LCF, SITE, OB, SAT, IM)
  implicit none
  include '../header/const.h'
  include '../header/station.h'
  include '../header/rnxobs.h'
  include '../header/satellite.h'
  include 'lsqcfg.h'
  include '../header/ionex.h'

! common
  integer*4     idxfrq(MAXSYS, 2)
  common        idxfrq
! parameter
  integer*4     jd
  real*8        sod
  type(lsqcfg)  LCF
  type(station) SITE
  type(rnxobr)  OB
  type(satellite) SAT(MAXSAT)
  type(ionex)   IM
! local
  integer*4     i0, i, j, k, isat, ite, ierr
  integer*4     iy, doy, jdc, jd_send, jd_recv
  real*8        fjd, sodc, sod_send, sod_recv, date
  real*8        drecclk, drecclk_tmp, dsatclk, dsatclk0, dsatclkrat, drate
  real*8        r1leng, r2leng, phase1, phase2, range1, range2, ddelay, reldel, dphwp
  real*8        gast, xpole, ypole
  real*8        sitrad, satrad, nadir
  real*8        freq(2), delay(2), pcv(2)
  real*8        dump(3), dx(3), r1(3), r2(3), rf(3), dloudx(3)
  real*8        xsun(6), xlun(6)
  real*8        trpdel, trpart(3), stec(4)
  real*8        xant_f(6, 2), xant_j(6, 2), xsat_f(3), xsat_j(6, 2), xpco_j(3, 2)
  real*8        rot_f2j(3, 3), rot_rat(3, 3), rot_l2j(3, 3)
  logical*1     flag(MAXSAT)
! function called
  real*8        dot

!
!! initial
  dump = 0.d0
  xant_j  = 0.d0
  xant_f  = 0.d0
  rot_f2j = 0.d0
  rot_rat = 0.d0
  xsun = 0.d0
  xlun = 0.d0
  xsat_j  = 0.d0
  xsat_f  = 0.d0
  dx = 0.d0
  r1 = 0.d0
  r2 = 0.d0
  rf = 0.d0
  delay   = 0.d0
  trpart  = 0.d0
  dloudx  = 0.d0
  pcv  = 0.d0
!
!! frequency
  i0 = index(GNSS_PRIO, 'G')
  freq(1) = FREQ_SYS(idxfrq(i0, 1), i0)
  freq(2) = FREQ_SYS(idxfrq(i0, 2), i0)
!
!! convert antenna offset from east-north-up to x-y-z in earth fixed system
  do j = 1, 2
    dump(1:3) = SITE%enu0(1:3) + SITE%enu_sys(1:3, j, i0)
    call matmpy(SITE%rot_l2f, dump, dump, 3, 3, 1)
    if (SITE%skd(1:1) .eq. 'L')  dump(1:3) = 0.d0
    xant_f(1:3, j) = SITE%x(1:3) + dump(1:3)*1.d-3
    xant_f(4:6, j) = 0.d0
  end do
!
!! receiver clock correction
  drecclk = 0.d0
  drecclk = SITE%rclock_G/VLIGHT
!
!! if receiver clock is too bad we repeat model starting from here
  ite = 0
  flag = .true.
100 continue
  ite = ite + 1
  call timinc(OB%jd, OB%tsec, -drecclk, jd_recv, sod_recv)
!
!! compute transformation matrix from earth-fixed to inertial system
  call ef2int(LCF%flnerp, jd_recv, sod_recv, rot_f2j, rot_rat, gast, xpole, ypole)
  do i = 1, 2
    call matmpy(rot_f2j, xant_f(1, i), xant_j(1, i), 3, 3, 1)
    call matmpy(rot_f2j, xant_f(4, i), dx, 3, 3, 1)
    call matmpy(rot_rat, xant_f(1, i), xant_j(4, i), 3, 3, 1)
    xant_j(4:6, i) = xant_j(4:6, i) + dx(1:3)
  end do
!
!! calculate solar and lunar coordinates
  call  sunxyz(jd_recv, sod_recv - 0.075d0, xsun)
  call moonxyz(jd_recv, sod_recv - 0.075d0, xlun)
  xsun = xsun/1.d3
  xlun = xlun/1.d3
  call matmpy(rot_f2j, xsun, xsun, 3, 3, 1)
  call matmpy(rot_f2j, xlun, xlun, 3, 3, 1)
!
!! tide displacement
  call tide_displace(LCF%tide, jd_recv, sod_recv, xant_j(1, 1), xant_f(1, 1), xsun, &
                     xlun, rot_f2j, SITE%rot_l2f, SITE%geod(1), SITE%geod(2), gast, xpole, ypole, SITE%olc, dx)
  do j = 1, 2
    xant_j(1:3, j) = xant_j(1:3, j) + dx(1:3)
  end do
!
!! LEO attitude
  SITE%xsrf=0.d0
  SITE%ysrf=0.d0
  SITE%zsrf=0.d0
  if (SITE%skd(1:1) .eq. 'L' .and. LCF%latuse) then
    call read_lat(LCF%flnlat, jd, sod, rot_f2j, SITE%xsrf, SITE%ysrf, SITE%zsrf, rot_l2j)
    do j = 1, 2
        do i = 1, 3
          xant_j(i, j) = xant_j(i, j) + (SITE%enu_sys(1, j, i0)*SITE%xsrf(i) + &
                                         SITE%enu_sys(2, j, i0)*SITE%ysrf(i) + &
                                         SITE%enu_sys(3, j, i0)*SITE%zsrf(i))*1.d-3
        end do
    end do
  end if
!
!! loop over all satellites
  do isat = 1, LCF%nprn
    if (LCF%prn(isat)(1:1) .ne. 'G') cycle
    if (OB%obs(isat, 1) .eq. 0.d0 .or. OB%obs(isat, 3) .eq. 0.d0 .or. .not. flag(isat)) cycle
!
!! get satellite clock information
    call read_satclk(LCF%flnsck, OB%prn(isat), jd, sod, jdc, sodc, dsatclk0, dsatclkrat, ierr)
!
!### BLOCK 1
!### ITERATION OF SEND TIME and GET THE GEOMETRIC DISTANCE
    ddelay = 0.1d0
    do while (dabs(ddelay) .gt. 1.d-9)
      call timinc(OB%jd, OB%tsec, -drecclk - OB%delay(isat), jd_send, sod_send)
!
!! compute the satellite position and velocity coordinates
      call lagrange_interp_orbit(LCF%flnorb, .true., .true., OB%jd, jd_send, sod_send, LCF%prn(isat), xsat_j(1, 1), xsat_j(4, 1))
      if (all(xsat_j(1:6, 1) .eq. 1.d15)) exit
      xsat_j(1:6, 2) = xsat_j(1:6, 1)
!
!! the satellite unit vectors (should be corrected for yaw error)
      SAT(isat)%xscf = 0.d0
      SAT(isat)%yscf = 0.d0
      SAT(isat)%zscf = 0.d0
      if (LCF%attuse) then
        call read_att(LCF%flnatt, LCF%prn(isat), jd_send, sod_send, rot_f2j, SAT(isat)%xscf, SAT(isat)%yscf, SAT(isat)%zscf)
      end if
      if (all(SAT(isat)%xscf(1:3) .eq. 0.d0)) then
        call rot_scfix2j2000(SAT(isat)%typ, xsat_j(1, 1), xsun, SAT(isat)%xscf, SAT(isat)%yscf, SAT(isat)%zscf)
      end if
      xpco_j = 0.d0
      do j = 1, 2
        do i = 1, 3
          xpco_j(i, j) = xpco_j(i, j) + (SAT(isat)%xyz(1, j)*SAT(isat)%xscf(i) + &
                                         SAT(isat)%xyz(2, j)*SAT(isat)%yscf(i) + &
                                         SAT(isat)%xyz(3, j)*SAT(isat)%zscf(i))*1.d-3
        end do
        xsat_j(1:3 ,j) = xsat_j(1:3, j) + xpco_j(1:3, j)
      end do
!
!! geometric distance
      r1(1:3) = xsat_j(1:3, 1) - xant_j(1:3, 1)
      r1leng = dsqrt(dot(3, r1, r1))
      dloudx(1:3) = r1(1:3)/r1leng
      delay(1) = r1leng/VLIGHT * 1.d3
!
!! decide whether to iterate the delay calculation
      ddelay = delay(1) - OB%delay(isat)
      OB%delay(isat) = delay(1)
    end do

    if (all(xsat_j(1:6, 1) .eq. 1.d15)) then
      flag(isat) = .false.
      OB%omc(isat, 1:4) = 0.d0
      OB%var(isat, 1) = -10.d0
      cycle      ! satellite missing
    end if
    r2(1:3) = xsat_j(1:3, 2) - xant_j(1:3, 2)
    r2leng = dsqrt(dot(3, r2, r2))
    delay(2) = r2leng/VLIGHT * 1.d3
!### END OF BLOCK 1
!###
!### BLOCK 2
!### COMPUTE CORRECTON
!
!! compute the atmospheric corrections to the delay
!  compute the station-satellite elevation angle and azimuth from north
    call matmpy(r1, rot_f2j, r2, 1, 3, 3)
    call matmpy(r2, SITE%rot_l2f, dump, 1, 3, 3)
    OB%azim(isat) = datan2(dump(1), dump(2))
    OB%elev(isat) = datan(dump(3)/dsqrt(dump(1)**2 + dump(2)**2))
    if (ierr .ne. 0) then
      flag(isat) = .false.
      OB%omc(isat, 1:4) = 0.d0
      OB%var(isat, 1) = -20.d0
      cycle
    end if
!--x
!! magnitude of magnetic field & theta
    stec(1:4) = 0.d0
    if (LCF%lioh) then
      !! ionosphere pierce point
      call ipp(SITE%x, OB%elev(isat), OB%azim(isat), IM%bradius, &
               IM%height, OB%ilat(isat), OB%ilon(isat), OB%zeni(isat))
      call mjd2doy(jd, iy, doy)
      date = iy + (doy + sod/86400.d0)/365.25d0
      call matmpy(xsat_j(1:3, 1), rot_f2j, xsat_f, 1, 3, 3)
      rf(1:3) = xant_f(1:3, 1) - xsat_f(1:3)
      call calB0Theta(rf, date, IM%bradius + IM%height, OB%ilat(isat)/pi*180.d0, OB%ilon(isat)/pi*180.d0, &
                      OB%ib0(isat), OB%itheta(isat))
      call calstech(jd, sod, isat, SAT, IM, OB, stec)
    end if
!
!! nadir angle for GNSS satellite
    nadir = dot(3, xsat_j(1, 1), r1)/dsqrt(dot(3, xsat_j(1, 1), xsat_j(1, 1)))/dsqrt(dot(3, r1, r1))
    nadir = dacos(nadir)
    SAT(isat)%nadir = nadir
!
!! get the nominal zenith delay, mapping function, and compute the path delay
    trpdel = 0.d0; trpart = 0.d0
    call troposphere_map(jd, sod, OB%elev(isat), SITE, OB%dmap(isat), OB%wmap(isat))
    trpdel = (OB%dmap(isat)*OB%zdd + OB%wmap(isat)*OB%zwd)/VLIGHT
    trpart(1) = OB%wmap(isat)
    if (LCF%htgmod(1:3) .ne. 'NON') then
      trpart(2) = OB%wmap(isat)/dtan(max(OB%elev(isat), 1.d-3))*dcos(OB%azim(isat))
      trpart(3) = OB%wmap(isat)/dtan(max(OB%elev(isat), 1.d-3))*dsin(OB%azim(isat))
      trpdel = trpdel + OB%nhtg*trpart(2)/VLIGHT + OB%ehtg*trpart(3)/VLIGHT
    end if
!
!! compute antenna / transmitter orientation dependent phase (cycle)
!! corrections for Right Circularly Polarized electro magnetic waves
    dphwp = 0.d0
    if (SITE%skd(1:1) .ne. 'L') then
      call phase_windup(SITE%first(isat), rot_f2j, SITE%rot_l2f, SAT(isat)%xscf, &
                        SAT(isat)%yscf, SAT(isat)%zscf, r1, SITE%prephi(isat), dphwp)
    else
      call lphase_windup(SITE%first(isat), rot_l2j, SAT(isat)%xscf, &
                         SAT(isat)%yscf, SAT(isat)%zscf, r1, SITE%prephi(isat), dphwp)
    end if
!
!! general relativistic time delay due to the Earth gravity (second)
    reldel = 2.d0*dot(3, xsat_j(1, 1), xsat_j(4, 1))/VLIGHT**2*1.d6
!
!! gravitional change effects on the satellite oscillators.
    sitrad = dsqrt(dot(3, xant_j(1, 1), xant_j(1, 1)))
    satrad = dsqrt(dot(3, xsat_j(1, 1), xsat_j(1, 1)))
    reldel = reldel + 2.d0*GM/(VLIGHT/1.d3)**3*log((sitrad + satrad + r1leng)/(sitrad + satrad - r1leng))
!
!! pcv correction for satellite and receiver antenna
    call get_ant_pcv(SITE%iptatx, SAT(isat)%iptatx, PI/2.d0 - OB%elev(isat), OB%azim(isat), nadir, pcv, LCF%prn(isat)(1:1))
!
!! antenna phase center correction (PCO + PCV)
    do j = 1, 2
      SAT(isat)%pcc(j) = -dot(3, xpco_j(1, j), r1)/sqrt(dot(3, r1, r1))*1.d3 - pcv(j)
    end do
!
!! delay
    delay(1) = delay(1) + trpdel + pcv(1)/VLIGHT + reldel + dphwp/freq(1)
    delay(2) = delay(2) + trpdel + pcv(2)/VLIGHT + reldel + dphwp/freq(2)
!
!! satellite clock correction, and save it for output clock estimates
    dsatclk = (jd - jdc)*86400.d0 + (sod - sodc)
    dsatclk = dsatclk0 + dsatclkrat*dsatclk
    phase1 = -VLIGHT*(dsatclk - drecclk - delay(1))
    phase2 = -VLIGHT*(dsatclk - drecclk - delay(2))
    range1 = phase1 - dphwp/freq(1)*VLIGHT - pcv(1)
    range2 = phase2 - dphwp/freq(2)*VLIGHT - pcv(2)
!
!! finally, form observed minus calculated, in meters
    OB%omc(isat, 1) = OB%obs(isat, 1)/freq(1)*VLIGHT - phase1 + stec(1)
    OB%omc(isat, 2) = OB%obs(isat, 2)/freq(2)*VLIGHT - phase2 + stec(2)
    OB%omc(isat, 3) = OB%obs(isat, 3) - range1 - stec(3)
    OB%omc(isat, 4) = OB%obs(isat, 4) - range2 - stec(4)
!
!! weight of observations
    OB%var(isat, 1) = (SITE%sigp*SAT(isat)%var*(vlight/freq(1)))**2
    OB%var(isat, 2) = (SITE%sigp*SAT(isat)%var*(vlight/freq(1)))**2
    OB%var(isat, 3) = (SITE%sigr*SAT(isat)%var)**2
    OB%var(isat, 4) = (SITE%sigr*SAT(isat)%var)**2
    if (OB%elev(isat)/PI*180.d0 .le. 30.d0) then
      OB%var(isat, 1:4) = OB%var(isat, 1:4)/(2*dsin(OB%elev(isat)))**2
    end if
!
!! compute the delay rate and the predicted delay
    dump(1:3) = xsat_j(4:6, 1) - xant_j(4:6, 1)
    drate = dot(3, dump, r1)/(VLIGHT/1.d3*r1leng)
    OB%delay(isat) = delay(1) + drate*LCF%dintv
!
!! observation equation
    call partial_gnss(OB%npar, dloudx, drate, OB%pname, OB%ltog(1, isat), rot_f2j, trpart, OB%amat(1, isat), 'G')
!
!! next satellite
  end do
!
!! receiver clock offset
  if (ite .lt. 4) then
    k = 0
    drecclk_tmp = 0.d0
    do isat = 1, LCF%nprn
      if (LCF%prn(isat)(1:1) .ne. 'G') cycle
      if (OB%omc(isat, 3) .ne. 0.d0) k = k + 1
      drecclk_tmp = drecclk_tmp + OB%omc(isat, 3)*freq(1)/VLIGHT
    end do
    drecclk = drecclk_tmp
    if (k .ge. 1) drecclk = drecclk/(k*freq(1))
    if (dabs(drecclk) .gt. 1.d-6 .or. (k .ge. 1 .and. SITE%rclock_G .eq. 0.d0)) then
      SITE%rclock_G = SITE%rclock_G + drecclk*VLIGHT
      drecclk = SITE%rclock_G/VLIGHT
      if (dabs(drecclk) .lt. 1.d-1) goto 100
      write (*, '(a,i7,f9.2,e15.4)') '***ERROR(gpsmod): abnormal drecclk at ', jd, sod, drecclk
      call exit(1)
    end if
  end if

  return
end subroutine
