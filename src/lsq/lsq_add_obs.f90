!
!! lsq_add_obs.f90
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

  integer*4 lfnrem, lfnobs, lfncid, jd
  real*8 sod
  type(lsqcfg) LCF
  type(rnxobr) OB
  type(prmt) PM(1:*)
  type(norm) NM
  type(satellite) SAT(MAXSAT)
  type(station) SITE
!
!! local
  logical*1 lfirst
  integer*4 ind, isat, ipar, i, k, ir, ic, nelem, ipt(0:MAXPAR)
  real*8 rwl, phase, range, wphs, wrng, wele, amat(MAXPAR)
  ! G
  real*8 g_G,g1_G,lamdw_G,coef_G
  ! R
  integer*4 frequency_glo_nu,prn_int
  real*8 :: FREQ1_R(-50:50),FREQ2_R(-50:50)
  real*8 g_R,g1_R,lamdw_R(-50:50),coef_R
  ! E
  real*8 g_E,g1_E,lamdw_E,coef_E
  ! C
  real*8 g_C,g1_C,lamdw_C,coef_C
  ! J
  real*8 g_J,g1_J,lamdw_J,coef_J

  real*8 pco_sit_obs(2),pco_sat_obs(2),pco_f1,pco_f2
  data lfirst/.true./
  save lfirst, g_G,g1_G,lamdw_G, g_R,g1_R,lamdw_R,g_E,g1_E,lamdw_E,g_C,g1_C,lamdw_C,g_J,g1_J,lamdw_J,&
       coef_G,coef_R,coef_E,coef_C,coef_J

  if (lfirst) then
    lfirst = .false.
    ! G
    g_G=FREQ1_G/FREQ2_G
    g1_G=g_G/(g_G*g_G-1.d0)
    lamdw_G=VLIGHT/FREQ1_G/(1.d0-1.d0/g_G)
    coef_G=g_G**2
    ! R
    g_R=9.0d0/7.0d0
    g1_R=g_R/(g_R*g_R-1.d0)
    coef_R=g_R**2
    ! E
    g_E=FREQ1_E/FREQ2_E
    g1_E=g_E/(g_E*g_E-1.d0)
    lamdw_E=VLIGHT/FREQ1_E/(1.d0-1.d0/g_E)
    coef_E=g_E**2
    ! C
    g_C=FREQ1_C/FREQ2_C
    g1_C=g_C/(g_C*g_C-1.d0)
    lamdw_C=VLIGHT/FREQ1_C/(1.d0-1.d0/g_C)
    coef_C=g_C**2
    ! J
    g_J=FREQ1_J/FREQ2_J
    g1_J=g_J/(g_J*g_J-1.d0)
    lamdw_J=VLIGHT/FREQ1_J/(1.d0-1.d0/g_J)
    coef_J=g_J**2
  endif
  call frequency_glonass(FREQ1_R,FREQ2_R)

  do isat = 1, LCF%nprn
    if (OB%omc(isat, 1) .eq. 0.d0 .or. OB%omc(isat, 3) .eq. 0.d0) cycle
    range = 0.d0; wrng = 0.d0
    phase = 0.d0; wphs = 0.d0
!
!! right hand side
    NM%nobs = NM%nobs + 1
    wrng = 1.d0/(OB%var(isat, 3) + OB%var(isat, 4))
    if(LCF%prn(isat)(1:1).eq.'G')then
      range=(OB%omc(isat,3)*coef_G-OB%omc(isat,4))/(coef_G-1.d0)
    elseif(LCF%prn(isat)(1:1).eq.'R')then
      range=(OB%omc(isat,3)*coef_R-OB%omc(isat,4))/(coef_R-1.d0)
    elseif(LCF%prn(isat)(1:1).eq.'E')then
      range=(OB%omc(isat,3)*coef_E-OB%omc(isat,4))/(coef_E-1.d0)
    elseif(LCF%prn(isat)(1:1).eq.'C')then
      range=(OB%omc(isat,3)*coef_C-OB%omc(isat,4))/(coef_C-1.d0)
    elseif(LCF%prn(isat)(1:1).eq.'J')then
      range=(OB%omc(isat,3)*coef_J-OB%omc(isat,4))/(coef_J-1.d0)
    end if
    NM%nobs = NM%nobs + 1
    wphs = 1.d0/(OB%var(isat, 1) + OB%var(isat, 2))
    if(LCF%prn(isat)(1:1).eq.'G')then
      phase=(OB%omc(isat,1)*coef_G-OB%omc(isat,2))/(coef_G-1.d0)
    elseif(LCF%prn(isat)(1:1).eq.'R')then
      phase=(OB%omc(isat,1)*coef_R-OB%omc(isat,2))/(coef_R-1.d0)
    elseif(LCF%prn(isat)(1:1).eq.'E')then
      phase=(OB%omc(isat,1)*coef_E-OB%omc(isat,2))/(coef_E-1.d0)
    elseif(LCF%prn(isat)(1:1).eq.'C')then
      phase=(OB%omc(isat,1)*coef_C-OB%omc(isat,2))/(coef_C-1.d0)
    elseif(LCF%prn(isat)(1:1).eq.'J')then
      phase=(OB%omc(isat,1)*coef_J-OB%omc(isat,2))/(coef_J-1.d0)
    end if
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
      if(LCF%prn(isat)(1:1).eq.'G') PM(ind)%iobs_G = PM(ind)%iobs_G + 1
      if(LCF%prn(isat)(1:1).eq.'R') PM(ind)%iobs_R = PM(ind)%iobs_R + 1
      if(LCF%prn(isat)(1:1).eq.'E') PM(ind)%iobs_E = PM(ind)%iobs_E + 1
      if(LCF%prn(isat)(1:1).eq.'C') then
        read(LCF%prn(isat)(2:3),'(i2)') prn_int
        if(prn_int .le. 17) then
          PM(ind)%iobs_C = PM(ind)%iobs_C + 1
        else
          PM(ind)%iobs_3 = PM(ind)%iobs_3 + 1
        endif
      endif
      if(LCF%prn(isat)(1:1).eq.'J') PM(ind)%iobs_J = PM(ind)%iobs_J + 1
!
!! initial value for the new ambiguity parameters
      pco_f1=0.d0
      pco_f2=0.d0
      if (OB%pname(ipar) (1:4) .eq. 'AMBC') then
        if(LCF%prn(isat)(1:1).eq.'G')then
          if(LCF%pcowl) then
            pco_sit_obs(1)=SITE%enu_G(3,1)*dsin(OB%elev(isat)) !f1 on meter
            pco_sit_obs(2)=SITE%enu_G(3,2)*dsin(OB%elev(isat)) !f2 on meter
            pco_sat_obs(1)=SAT(isat)%xyz(3,1)*dcos(SAT(isat)%nadir) !f1 on meter
            pco_sat_obs(2)=SAT(isat)%xyz(3,2)*dcos(SAT(isat)%nadir) !f2 on meter
            pco_f1=pco_sat_obs(1)+pco_sit_obs(1)
            pco_f2=pco_sat_obs(2)+pco_sit_obs(2)
          endif
          rwl=(OB%obs(isat,1)+pco_f1*FREQ1_G/VLIGHT) - (OB%obs(isat,2)+pco_f2*FREQ2_G/VLIGHT) &
              -(g_G*(OB%obs(isat,3)+pco_f1) + (OB%obs(isat,4)+pco_f2))/(1.d0+g_G)/lamdw_G
        elseif(LCF%prn(isat)(1:1).eq.'R')then
          if(LCF%pcowl) then
            pco_sit_obs(1)=SITE%enu_R(3,1)*dsin(OB%elev(isat)) !f1 on meter
            pco_sit_obs(2)=SITE%enu_R(3,2)*dsin(OB%elev(isat)) !f2 on meter
            pco_sat_obs(1)=SAT(isat)%xyz(3,1)*dcos(SAT(isat)%nadir) !f1 on meter
            pco_sat_obs(2)=SAT(isat)%xyz(3,2)*dcos(SAT(isat)%nadir) !f2 on meter
            pco_f1=pco_sat_obs(1)+pco_sit_obs(1)
            pco_f2=pco_sat_obs(2)+pco_sit_obs(2)
          endif
          read(LCF%prn(isat),'(1x,i2)') prn_int
          frequency_glo_nu=OB%glschn(prn_int)
          lamdw_R(frequency_glo_nu)=VLIGHT/FREQ1_R(frequency_glo_nu)/(1.d0-1.d0/g_R)
          rwl=(OB%obs(isat,1)+pco_f1*FREQ1_R(frequency_glo_nu)/VLIGHT) &
              - (OB%obs(isat,2)+pco_f2*FREQ2_R(frequency_glo_nu)/VLIGHT) &
              - (g_R*(OB%obs(isat,3)+pco_f1) + (OB%obs(isat,4)+pco_f2))/(1.d0+g_R)/lamdw_R(frequency_glo_nu)
        elseif(LCF%prn(isat)(1:1).eq.'E')then
          if(LCF%pcowl) then
            pco_sit_obs(1)=SITE%enu_E(3,1)*dsin(OB%elev(isat)) !f1 on meter
            pco_sit_obs(2)=SITE%enu_E(3,2)*dsin(OB%elev(isat)) !f2 on meter
            pco_sat_obs(1)=SAT(isat)%xyz(3,1)*dcos(SAT(isat)%nadir) !f1 on meter
            pco_sat_obs(2)=SAT(isat)%xyz(3,2)*dcos(SAT(isat)%nadir) !f2 on meter
            pco_f1=pco_sat_obs(1)+pco_sit_obs(1)
            pco_f2=pco_sat_obs(2)+pco_sit_obs(2)
          endif
          rwl=(OB%obs(isat,1)+pco_f1*FREQ1_E/VLIGHT) - (OB%obs(isat,2)+pco_f2*FREQ2_E/VLIGHT) &
              -(g_E*(OB%obs(isat,3)+pco_f1) + (OB%obs(isat,4)+pco_f2))/(1.d0+g_E)/lamdw_E
        elseif(LCF%prn(isat)(1:1).eq.'C')then
          if(LCF%pcowl) then
            pco_sit_obs(1)=SITE%enu_C(3,1)*dsin(OB%elev(isat)) !f1 on meter
            pco_sit_obs(2)=SITE%enu_C(3,2)*dsin(OB%elev(isat)) !f2 on meter
            pco_sat_obs(1)=SAT(isat)%xyz(3,1)*dcos(SAT(isat)%nadir) !f1 on meter
            pco_sat_obs(2)=SAT(isat)%xyz(3,2)*dcos(SAT(isat)%nadir) !f2 on meter
            pco_f1=pco_sat_obs(1)+pco_sit_obs(1)
            pco_f2=pco_sat_obs(2)+pco_sit_obs(2)
          endif
          rwl=(OB%obs(isat,1)+pco_f1*FREQ1_C/VLIGHT) - (OB%obs(isat,2)+pco_f2*FREQ2_C/VLIGHT) &
              -(g_C*(OB%obs(isat,3)+pco_f1) + (OB%obs(isat,4)+pco_f2))/(1.d0+g_C)/lamdw_C
        elseif(LCF%prn(isat)(1:1).eq.'J')then
          if(LCF%pcowl) then
            pco_sit_obs(1)=SITE%enu_J(3,1)*dsin(OB%elev(isat)) !f1 on meter
            pco_sit_obs(2)=SITE%enu_J(3,2)*dsin(OB%elev(isat)) !f2 on meter
            pco_sat_obs(1)=SAT(isat)%xyz(3,1)*dcos(SAT(isat)%nadir) !f1 on meter
            pco_sat_obs(2)=SAT(isat)%xyz(3,2)*dcos(SAT(isat)%nadir) !f2 on meter
            pco_f1=pco_sat_obs(1)+pco_sit_obs(1)
            pco_f2=pco_sat_obs(2)+pco_sit_obs(2)
          endif
          rwl=(OB%obs(isat,1)+pco_f1*FREQ1_J/VLIGHT) - (OB%obs(isat,2)+pco_f2*FREQ2_J/VLIGHT) &
              -(g_J*(OB%obs(isat,3)+pco_f1) + (OB%obs(isat,4)+pco_f2))/(1.d0+g_J)/lamdw_J
        end if
        wele = 1.d0
        if (OB%elev(isat)*180.d0/PI .le. 30.d0) wele = 2.d0*dsin(OB%elev(isat))
        ipt(0) = nelem
        if (OB%flag(isat) .eq. 1) then
          PM(ind)%xini = phase - range
          PM(ind)%zw = rwl
        endif
        phase = phase - amat(nelem)*PM(ind)%xini
!
!! processing of wide-lane observation
        rwl = rwl - PM(ind)%zw
        PM(ind)%xrwl = PM(ind)%xrwl + rwl*wele
        PM(ind)%rw = PM(ind)%rw + wele
        PM(ind)%xswl = PM(ind)%xswl + wele*rwl**2
        PM(ind)%mele = PM(ind)%mele + OB%elev(isat)
      endif
!
!! next parameters
    enddo
!
!! save to recover residuals
    if (lfnobs .ne. 0) then
      write (lfncid) 'ob'
      write (lfnobs) jd, sod, isat, nelem, ipt(0), (ipt(i), amat(i), i=1, nelem), &
        phase, range, wphs, wrng, OB%flag(isat), OB%elev(isat), OB%azim(isat), OB%dmap(isat), OB%wmap(isat),&
        OB%typuse(isat,1:4)
    endif
!
!! transform to index of normal matrix
    do i = 1, nelem
      ipt(i) = PM(ipt(i))%ipt
    enddo
!
!! upper normal matrix is preferred
    do i = 1, nelem
      do k = i, nelem
        ir = min(ipt(i), ipt(k))
        ic = max(ipt(i), ipt(k))
        NM%norx(ir, ic) = NM%norx(ir, ic) + amat(i)*amat(k)*wphs
        if (i .ne. ipt(0) .and. k .ne. ipt(0)) NM%norx(ir, ic) = NM%norx(ir, ic) + amat(i)*amat(k)*wrng
      enddo
      ir = ipt(i)
      NM%norx(ir, NM%imtx + 1) = NM%norx(ir, NM%imtx + 1) + amat(i)*wphs*phase
      if (i .ne. ipt(0)) NM%norx(ir, NM%imtx + 1) = NM%norx(ir, NM%imtx + 1) + amat(i)*wrng*range
    enddo
    NM%ltpl = NM%ltpl + wphs*phase**2
    NM%ltpl = NM%ltpl + wrng*range**2
!
!! next satellite
  enddo

  return
end
