!
!! lsq_add_obs.f90
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
!! purpose   : add observation equations to estimator LSQ
!! paraemters:
!!             lfncid,lfnobs -- tmp files for recovery
!!             jd,sod        -- float julian day
!!             LCF           -- LSQ control parameters
!!             OB            -- Observation struct
!!             PM,NM         -- normal matrix & PAR table
!
subroutine lsq_add_obs(lfncid, lfnobs, lfnrem, jd, sod, LCF, OB, PM, NM, SAT, SITE)
  implicit none
  include '../header/const.h'
  include '../header/rnxobs.h'
  include 'lsq.h'
  include 'lsqcfg.h'
  include '../header/station.h'
  include '../header/satellite.h'

! common
  integer*4     idxfrq(MAXSYS, 2)
  common        idxfrq
! parameter
  integer*4     lfncid, lfnobs, lfnrem
  integer*4     jd
  real*8        sod
  type(lsqcfg)  LCF
  type(rnxobr)  OB
  type(prmt)    PM(1:*)
  type(norm)    NM
  type(satellite) SAT(MAXSAT)
  type(station) SITE
! local
  logical*1     lfirst
  integer*4     i0, ind, isat, ipar, i, k, ir, ic, nelem, ipt(0:MAXPAR), prn_int, az, el
  real*8        rwl, phase, range, wphs, wrng, wele, amat(MAXPAR)
  real*8        f1(MAXSYS), f2(MAXSYS), rg(MAXSYS), r2(MAXSYS), g1(MAXSYS), lamw(MAXSYS)
  real*8        FREQ1_R(-50:50), FREQ2_R(-50:50)
  real*8        pco_sit_obs(2), pco_f1, pco_f2

  data lfirst/.true./
  save lfirst, f1, f2, rg, r2, g1, lamw

  if (lfirst) then
    lfirst = .false.
    do i0 = 1, MAXSYS
      f1(i0) = FREQ_SYS(idxfrq(i0, 1), i0)
      if (f1(i0) .eq. 0.d0) goto 100
      f2(i0) = FREQ_SYS(idxfrq(i0, 2), i0)
      if (f2(i0) .eq. 0.d0) goto 100
      rg(i0) = f1(i0)/f2(i0)
      r2(i0) = rg(i0)*rg(i0)
      g1(i0) = rg(i0)/(r2(i0) - 1)
      lamw(i0) = VLIGHT/(f1(i0) - f2(i0))
    end do
  end if

  call frequency_glonass(FREQ1_R, FREQ2_R)

  do isat = 1, LCF%nprn
    if (OB%omc(isat, 1) .eq. 0.d0 .or. OB%omc(isat, 3) .eq. 0.d0) cycle
    i0 = index(GNSS_PRIO, LCF%prn(isat)(1:1))
    range = 0.d0; wrng = 0.d0
    phase = 0.d0; wphs = 0.d0
!
!! get az and el	
	if(OB%azim(isat) .le. 0.d0) then
	  az=floor(OB%azim(isat)*180.d0/PI)+1
	else
	  az=floor(OB%azim(isat)*180.d0/PI)+360+1
	end if
	el=floor(OB%elev(isat)*180.d0/PI)+1
!
!! right hand side
    NM%nobs = NM%nobs + 1
    wrng = 1.d0/(OB%var(isat, 3) + OB%var(isat, 4))
    range = (OB%omc(isat, 3)*r2(i0) - OB%omc(isat, 4))/(r2(i0) - 1.d0)-SITE%mhm(az,el,1,i0)
    NM%nobs = NM%nobs + 1
    wphs = 1.d0/(OB%var(isat, 1) + OB%var(isat, 2))
    phase = (OB%omc(isat, 1)*r2(i0) - OB%omc(isat, 2))/(r2(i0) - 1.d0)-SITE%mhm(az,el,2,i0)
!
!! observation equations
    nelem = 0
    ipt(0) = 0
    do ipar = 1, OB%npar
      if (OB%ltog(ipar, isat) .eq. 0) cycle
!
!! extract non-zero ones
      nelem = nelem + 1
      ipt(nelem) = OB%ltog(ipar, isat)
      amat(nelem) = OB%amat(ipar, isat)
      ind = OB%ltog(ipar, isat)
      PM(ind)%iobs = PM(ind)%iobs + 1
      if (LCF%prn(isat)(1:1) .eq. 'G') PM(ind)%iobs_G = PM(ind)%iobs_G + 1
      if (LCF%prn(isat)(1:1) .eq. 'R') then
        PM(ind)%iobs_R = PM(ind)%iobs_R + 1
        read (LCF%prn(isat)(1:3), '(1x,i2)') prn_int
        f1(i0) = FREQ1_R(OB%glschn(prn_int))
        f2(i0) = FREQ2_R(OB%glschn(prn_int))
        lamw(i0) = VLIGHT/(f1(i0) - f2(i0))
      end if
      if (LCF%prn(isat)(1:1) .eq. 'E') PM(ind)%iobs_E = PM(ind)%iobs_E + 1
      if (LCF%prn(isat)(1:1) .eq. 'C') then
        read (LCF%prn(isat) (2:3), '(i2)') prn_int
        if (prn_int .le. 17) then
          PM(ind)%iobs_C = PM(ind)%iobs_C + 1
        else
          PM(ind)%iobs_3 = PM(ind)%iobs_3 + 1
        end if
      end if
      if (LCF%prn(isat)(1:1) .eq. 'J') PM(ind)%iobs_J = PM(ind)%iobs_J + 1
!
!! initial value for the new ambiguity parameters
      pco_f1 = 0.d0
      pco_f2 = 0.d0
      if (OB%pname(ipar)(1:4) .eq. 'AMBC') then
        if (LCF%pcowl) then
            pco_sit_obs(1) = SITE%enu_sys(3, 1, i0)*dsin(OB%elev(isat)) + &
                             SITE%enu_sys(1, 1, i0)*dcos(OB%elev(isat))*dsin(OB%azim(isat)) + &
                             SITE%enu_sys(2, 1, i0)*dcos(OB%elev(isat))*dcos(OB%azim(isat))
            pco_sit_obs(2) = SITE%enu_sys(3, 2, i0)*dsin(OB%elev(isat)) + &
                             SITE%enu_sys(1, 2, i0)*dcos(OB%elev(isat))*dsin(OB%azim(isat)) + &
                             SITE%enu_sys(2, 2, i0)*dcos(OB%elev(isat))*dcos(OB%azim(isat))
            pco_f1 = SAT(isat)%pcc(1) + pco_sit_obs(1)
            pco_f2 = SAT(isat)%pcc(2) + pco_sit_obs(2)
        end if
        rwl = (OB%obs(isat, 1) + pco_f1*f1(i0)/VLIGHT) - (OB%obs(isat, 2) + pco_f2*f2(i0)/VLIGHT) &
              - (rg(i0)*(OB%obs(isat, 3) + pco_f1) + (OB%obs(isat, 4) + pco_f2))/(1.d0 + rg(i0))/lamw(i0)
        wele = 1.d0
        if (OB%elev(isat)*180.d0/PI .le. 30.d0) wele = 2.d0*dsin(OB%elev(isat))
        ipt(0) = nelem
        if (OB%flag(isat) .eq. 1) then
          PM(ind)%xini = phase - range
          PM(ind)%zw = rwl
        end if
        phase = phase - amat(nelem)*PM(ind)%xini
!
!! processing of wide-lane observation
        rwl = rwl - PM(ind)%zw
        PM(ind)%xrwl = PM(ind)%xrwl + rwl*wele
        PM(ind)%rw = PM(ind)%rw + wele
        PM(ind)%xswl = PM(ind)%xswl + wele*rwl**2
        PM(ind)%mele = PM(ind)%mele + OB%elev(isat)
      end if
!
!! next parameters
    end do
!
!! save to recover residuals
    if (lfnobs .ne. 0) then
      write (lfncid) 'ob'
      write (lfnobs) jd, sod, isat, nelem, ipt(0), (ipt(i), amat(i), i=1, nelem), &
        phase, range, wphs, wrng, OB%flag(isat), OB%elev(isat), OB%azim(isat), OB%dmap(isat), OB%wmap(isat), &
        OB%typuse(isat, 1:4)
    end if
!
!! transform to index of normal matrix
    do i = 1, nelem
      ipt(i) = PM(ipt(i))%ipt
    end do
!
!! upper normal matrix is preferred
    do i = 1, nelem
      do k = i, nelem
        ir = min(ipt(i), ipt(k))
        ic = max(ipt(i), ipt(k))
        NM%norx(ir, ic) = NM%norx(ir, ic) + amat(i)*amat(k)*wphs
        if (i .ne. ipt(0) .and. k .ne. ipt(0)) NM%norx(ir, ic) = NM%norx(ir, ic) + amat(i)*amat(k)*wrng
      end do
      ir = ipt(i)
      NM%norx(ir, NM%imtx + 1) = NM%norx(ir, NM%imtx + 1) + amat(i)*wphs*phase
      if (i .ne. ipt(0)) NM%norx(ir, NM%imtx + 1) = NM%norx(ir, NM%imtx + 1) + amat(i)*wrng*range
    end do
    NM%ltpl = NM%ltpl + wphs*phase**2
    NM%ltpl = NM%ltpl + wrng*range**2
!
!! next satellite
  end do

  return
100 continue
  write (*, '(a,5(1x,a,2i1))') '***ERROR(lsq_add_obs): invalid frequency number:', &
    'G', idxfrq(index(GNSS_PRIO, 'G'), 1:2), &
    'R', idxfrq(index(GNSS_PRIO, 'R'), 1:2), &
    'E', idxfrq(index(GNSS_PRIO, 'E'), 1:2), &
    'C', idxfrq(index(GNSS_PRIO, 'C'), 1:2), &
    'J', idxfrq(index(GNSS_PRIO, 'J'), 1:2)
  call exit(1)
end subroutine
