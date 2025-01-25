!
!! lsq.f90
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
!! Least Squares Estimator
!
program lsq
  implicit none
  include '../header/const.h'
  include '../header/ionex.h'
  include '../header/orbit.h'
  include '../header/absbia.h'
  include '../header/rnxobs.h'
  include '../header/station.h'
  include '../header/satellite.h'
  include '../redig/data_flag.h'
  include 'lsqcfg.h'
  include 'lsq.h'

! common
  integer*4     idxfrq(MAXSYS, 2)
  common        idxfrq
!
!! lsq configure
  type(lsqcfg)  LCF
!! absolute bias
  type(absbia)  BIAS(MAXSAT, MAXTYP)
!! station
  type(station) SITE
  type(rnxhdr)  HD
  type(rnxobr)  OB
!! orbit
  type(orbhdr)  OH
  type(satellite) SAT(MAXSAT)
!! normal matrix
  type(norm)    NM
!! parameter definition
  type(prmt), pointer :: PM(:)
!! ionosphere
  type(ionex)   IM
!
!! local
  logical*1     lopen
  logical*1     arsig_used(MAXSAT)
  integer*4     i0, i, j, k, jd, jd_sav, jdc, isat, iepo, ipar, iamb, id, ih
  integer*4     openid, ierr
  integer*4     lfncid, lfnobs, lfnrem, lfnres, lfnpos, lfnamb, lfnneq, iunit_next
  integer*4     nbias_used(MAXSAT)
  integer*4     nbias, ibias(MAXPAR)
  integer*4     index_g(2), index_r(2), index_e(2), index_c(2), index_3(2), index_j(2)
  real*8        freq, dwnd, sod, sec, sodc, deltax(3), pdop, pdop_sum, nbb(4, 4) 
  character*3   temprn
  character*20  antnum
!! GLONASS
  integer*4     frequency_glo_nu, prn_int
  real*8        FREQ1_R(-50:50), FREQ2_R(-50:50)
!
!! function called
  integer*4     get_valid_unit
  integer*4     pointer_string
  real*8        timdif
  character*10  run_tim
!
!! initial
  do i = 1, MAXSAT
    LCF%prn(i) = ''
  end do
  LCF%attuse = .false.
  LCF%otluse = .false.
  LCF%mhmuse = .false.
  arsig_used = .false.
  call frequency_glonass(FREQ1_R, FREQ2_R)
!
!! get arguments
  call get_lsq_args(LCF, SITE, OB, SAT, IM)
  index_g = 0
  index_r = 0
  index_e = 0
  index_c = 0
  index_3 = 0
  index_j = 0
  do i = 1, LCF%nprn
    if (LCF%prn(i)(1:1) .eq. 'G') then
      if (index_g(1) .eq. 0) index_g(1) = i
      index_g(2) = i
    end if
    if (LCF%prn(i)(1:1) .eq. 'R') then
      if (index_r(1) .eq. 0) index_r(1) = i
      index_r(2) = i
    end if
    if (LCF%prn(i)(1:1) .eq. 'E') then
      if (index_e(1) .eq. 0) index_e(1) = i
      index_e(2) = i
    end if
    if (LCF%prn(i)(1:1) .eq. 'C') then
      read (LCF%prn(i) (2:3), '(i2)') j
      if (j .le. 17) then
        if (index_c(1) .eq. 0) index_c(1) = i
        index_c(2) = i
      end if
      if (j .gt. 17) then
        if (index_3(1) .eq. 0) index_3(1) = i
        index_3(2) = i
      end if
    end if
    if (LCF%prn(i)(1:1) .eq. 'J') then
      if (index_j(1) .eq. 0) index_j(1) = i
      index_j(2) = i
    end if
  end do
!
!! write removing info for recovering
  lfncid = get_valid_unit(10)
  open (lfncid, file='tmp_cid', form='unformatted')
  write (lfncid) '00'
  lfnobs = get_valid_unit(10)
  open (lfnobs, file='tmp_obs', form='unformatted')
  lfnrem = get_valid_unit(10)
  open (lfnrem, file='tmp_rem', form='unformatted')
!
!! satellite antenna: absolute phase centers
  LCF%flnatx_real = ' '
  do isat = 1, LCF%nprn
    write (antnum, '(a3)') SAT(isat)%prn
    call get_ant_ipt(LCF%jd0 + LCF%sod0/864.d2, LCF%jd1 + LCF%sod1/864.d2,   &
                     SAT(isat)%typ, antnum, SAT(isat)%iptatx, SAT(isat)%xyz, &
                     SAT(isat)%prn(1:1), LCF%flnatx_real)
  end do
!
!! get multipath model
  SITE%mhm = 0.0d0
  if (LCF%mhmuse) then
    print *, 1
    call read_mhm(SITE%mhm,SITE%name)  
  end if

!! site antenna
!! observation
  SITE%iunit = get_valid_unit(10)
  open (SITE%iunit, file=SITE%obsfil, status='old', iostat=ierr)
  if (ierr .ne. 0) then
    write (*, '(a)') '***ERROR(lsq): open file ', trim(SITE%obsfil)
    call exit(1)
  end if
  call rdrnxoh(SITE%iunit, HD, ierr)
  if (ierr .ne. 0) then
    write (*, '(2a)') '***ERROR(lsq): read header ', trim(SITE%obsfil)
    call exit(1)
  end if
  SITE%anttyp = HD%anttyp
  SITE%rectyp = HD%rectyp
  if (SITE%skd(1:1) .eq. 'L') then
    SITE%anttyp = SITE%name
    SITE%rectyp = ''
  end if
  SITE%enu0(1) = HD%e
  SITE%enu0(2) = HD%n
  SITE%enu0(3) = HD%h
  antnum = ' '
  do i0 = 1, MAXSYS
    call get_ant_ipt(LCF%jd0 + LCF%sod0/864.d2, LCF%jd1 + LCF%sod1/864.d2, &
                     SITE%anttyp, antnum, SITE%iptatx, SITE%enu_sys(:, :, i0), &
                     GNSS_PRIO(i0:i0), LCF%flnatx_real)
  end do
!
!! count number of parameters
  call lsq_cnt_prmt(LCF, SITE, NM)
  NM%npm = 0
  NM%npm = NM%npm + OB%amb_tot
  NM%npm = NM%npm + NM%imtx
  if (NM%npm .gt. MAXPAR) then
    write (*, '(a,i8)') '***ERROR(lsq): too many parameters ', NM%npm
    call exit(1)
  end if
  allocate (PM(NM%npm), stat=ierr)
  if (ierr .ne. 0) then
    write (*, '(a,i8)') '***ERROR(lsq): parameter array allocation ', NM%npm
    call exit(1)
  end if
!
!! normal matrix size
  if (LCF%lrmbias) then
    NM%nmtx = 0
    NM%nmtx = NM%nmtx + OB%amb_epo
    NM%nmtx = NM%nmtx + NM%imtx + 1
  else
    NM%nmtx = NM%npm + 1
  end if
!! allocate memory
  allocate (NM%norx(1:NM%nmtx, 1:NM%nmtx), stat=ierr)
  allocate (NM%iptp(NM%nmtx), stat=ierr)
  allocate (NM%iptx(NM%nmtx), stat=ierr)
  do i = 1, NM%nmtx
    NM%iptp(i) = 0
    NM%iptx(i) = NM%nmtx * (i - 1)
  end do
  write (*, '(/,3(a,i8,/))') ' size of Normal Matrix :', NM%nmtx, &
                             ' size of NC Parameter  :', NM%nc,   &
                             ' size of NP Parameter  :', NM%np
!
!! initiate normal matrix
  call lsq_init(LCF, SITE, OB, NM, PM)
!
!!++++++++++++++++++++EPOCH LOOP+++++++++++++++++++++++++++++
  jd = LCF%jd0
  sod = LCF%sod0
  iepo = 0
  jd_sav = 0
  NM%ltpl = 0.d0
  NM%nuk = NM%nc
  NM%nobs = 0
 ! dwnd = min(LCF%dintv/10.d0, 0.01d0)
  dwnd=LCF%dwnd
  pdop_sum = 0.d0
  do while (timdif(jd, sod, LCF%jd1, LCF%sod1) .lt. MAXWND)
    OB%omc(1:LCF%nprn, 1:4) = 0.d0
    if (jd .gt. jd_sav) then
      !
      !! reset site antenna info
      SITE%anttyp = HD%anttyp
      SITE%rectyp = HD%rectyp
      if (SITE%skd(1:1) .eq. 'L') then
        SITE%anttyp = SITE%name
        SITE%rectyp = ''
      end if
      SITE%enu0(1) = HD%e
      SITE%enu0(2) = HD%n
      SITE%enu0(3) = HD%h
      antnum = ' '
      do i0 = 1, MAXSYS
        call get_ant_ipt(LCF%jd0 + LCF%sod0/864.d2, LCF%jd1 + LCF%sod1/864.d2, &
                         SITE%anttyp, antnum, SITE%iptatx, SITE%enu_sys(:, :, i0), &
                         GNSS_PRIO(i0:i0), LCF%flnatx_real)
      end do
      !
      !! read PPP-AR biases
      call read_bias(LCF%flnfcb, LCF%nprn, LCF%prn, BIAS, jd*1.d0, (jd + 1)*1.d0)
      jd_sav = jd
    end if
!
!! read one epoch of observations
50 continue
    OB%typuse = ' '
    if (SITE%iunit .eq. 0) cycle
    if (HD%ver .ge. 300) then
      call rdrnxoi3(SITE%iunit, jd, sod, dwnd, LCF%nprn, LCF%prn, HD, OB, BIAS, nbias_used, ierr)
    else
      call rdrnxoi2(SITE%iunit, jd, sod, dwnd, LCF%nprn, LCF%prn, HD, OB, BIAS, nbias_used, ierr)
    end if
    if (ierr .eq. 2) then
      call next_rinex(SITE%iunit, iunit_next, LCF%jd1)
      if (iunit_next .eq. 0) goto 45
      SITE%iunit = iunit_next
      goto 50
    end if
45 continue
    if (ierr .ne. 0) then
      write (*, '(a)') '###WARNING(lsq): no more data to be processed'
      exit
    end if
    call read_obsrhd(jd, sod, LCF%nprn, LCF%prn, OB)
    if (SITE%lfnjmp .ne. 0) call read_clkjmp(jd, sod, SITE%lfnjmp, LCF%nprn, OB)
!
!! +++++++++++++++++++++++++++SITE LOOP+++++++++++++++++++++++++++++++++++
!! add observation equations
!
!! read kinematic position: if not found, use value at last epoch
    if ((count(OB%obs(1:LCF%nprn, 3) .ne. 0.d0) .gt. 0) .and. &
        (SITE%skd(1:1) .eq. 'K' .or. SITE%skd(1:1) .eq. 'L')) then
      ipar = pointer_string(OB%npar, OB%pname, 'STAPX')
      ipar = OB%ltog(ipar, 1)
      if (PM(ipar)%iepo .eq. 0) then
        call read_kinpos(SITE, PM(ipar)%ptime(1), PM(ipar)%ptime(2), deltax(1), deltax(2), deltax(3), dwnd)
        if (all(deltax(1:3) .eq. 1.d0) .or. &
            all(deltax(1:3) .eq. 0.d0)) then
          goto 40
        end if
        do j = 0, 2 
          PM(ipar + j)%xini = deltax(j + 1)
        end do
        if (SITE%skd(1:1) .ne. 'L') then
          call xyzblh(deltax(1:3), 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, deltax)
          if (deltax(3) .lt. -4.d2 .or. deltax(3) .gt. 20.d3) goto 40
        end if
      end if
    end if
!
!! add new ambiguities
    call lsq_add_newamb(jd, sod, SITE%name, LCF, OB, NM, PM)
!
!! iteration for more precise position
    if (count(OB%obs(1:LCF%nprn, 3) .ne. 0.d0) .eq. 0) goto 40
    iepo = iepo + 1
    if (dmod(iepo*1.d0, dnint(300.d0/LCF%dintv)) .le. MAXWND) write (*, '(2a,i9,i7,f9.1)') run_tim(), ' Epoch ', iepo, jd, sod
!
!! prepare a priori information for kinematic station
    if (SITE%skd(1:1) .eq. 'K' .or. SITE%skd(1:1) .eq. 'L') then
      ipar = pointer_string(OB%npar, OB%pname, 'STAPX')
      do i = 0, 2
        SITE%x(i + 1) = PM(OB%ltog(ipar + i, 1))%xini * 1.d-3
      end do
      call xyzblh(SITE%x(1:3)*1.d3, 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, SITE%geod)
      SITE%geod(3) = SITE%geod(3)*1.d-3
      call rot_enu2xyz(SITE%geod(1), SITE%geod(2), SITE%rot_l2f)
    end if
!
!! prepare a priori receiver clock correction (unit: m)
    ipar = pointer_string(OB%npar, OB%pname, 'RECCLK_G'//trim(LCF%rckmod))
    if (ipar .ne. 0) SITE%rclock_G = PM(OB%ltog(ipar, index_g(1)))%xini
    ipar = pointer_string(OB%npar, OB%pname, 'RECCLK_R'//trim(LCF%rckmod))
    if (ipar .ne. 0) SITE%rclock_R = PM(OB%ltog(ipar, index_r(1)))%xini
    ipar = pointer_string(OB%npar, OB%pname, 'RECCLK_E'//trim(LCF%rckmod))
    if (ipar .ne. 0) SITE%rclock_E = PM(OB%ltog(ipar, index_e(1)))%xini
    ipar = pointer_string(OB%npar, OB%pname, 'RECCLK_C'//trim(LCF%rckmod))
    if (ipar .ne. 0) SITE%rclock_C = PM(OB%ltog(ipar, index_c(1)))%xini
    ipar = pointer_string(OB%npar, OB%pname, 'RECCLK_3'//trim(LCF%rckmod))
    if (ipar .ne. 0) SITE%rclock_3 = PM(OB%ltog(ipar, index_3(1)))%xini
    ipar = pointer_string(OB%npar, OB%pname, 'RECCLK_J'//trim(LCF%rckmod))
    if (ipar .ne. 0) SITE%rclock_J = PM(OB%ltog(ipar, index_j(1)))%xini
!
!! one way model
    call read_meteo(jd, sod, SITE%imet, SITE%map, SITE%geod, SITE%p0, SITE%t0, SITE%hr0, SITE%undu)
    OB%zdd = 0.d0
    OB%zwd = 0.d0
    OB%nhtg = 0.d0
    OB%ehtg = 0.d0
    if (SITE%skd(1:1) .ne. 'L')  call troposphere_delay(jd, sod, SITE, OB%zdd, OB%zwd)
    if (index(LCF%sys, 'G') .ne. 0 .and. count(OB%obs(index_g(1):index_g(2), 3) .ne. 0.d0) .ne. 0) &
      call gpsmod(jd, sod, LCF, SITE, OB, SAT, IM)
    if (index(LCF%sys, 'R') .ne. 0 .and. count(OB%obs(index_r(1):index_r(2), 3) .ne. 0.d0) .ne. 0) &
      call glomod(jd, sod, LCF, SITE, OB, SAT, IM)
    if (index(LCF%sys, 'E') .ne. 0 .and. count(OB%obs(index_e(1):index_e(2), 3) .ne. 0.d0) .ne. 0) &
      call galmod(jd, sod, LCF, SITE, OB, SAT, IM)
    if (index(LCF%sys, 'C') .ne. 0 .and. count(OB%obs(index_c(1):index_c(2), 3) .ne. 0.d0) .ne. 0) &
      call bd2mod(jd, sod, LCF, SITE, OB, SAT, IM)
    if (index(LCF%sys, '3') .ne. 0 .and. count(OB%obs(index_3(1):index_3(2), 3) .ne. 0.d0) .ne. 0) &
      call bd3mod(jd, sod, LCF, SITE, OB, SAT, IM)
    if (index(LCF%sys, 'J') .ne. 0 .and. count(OB%obs(index_j(1):index_j(2), 3) .ne. 0.d0) .ne. 0) &
      call qzsmod(jd, sod, LCF, SITE, OB, SAT, IM)
!
!! check elevation & cutoff angle
    do isat = 1, LCF%nprn
      if (all(OB%obs(isat, 1:4) .eq. 0.d0)) cycle
      if (any(OB%omc(isat, 1:4) .eq. 0.d0)) then
        if (OB%var(isat, 1) .eq. -10.d0) then
          write (*, '(a,i5,f10.2,a,a3)') '###WARNING(lsq): no orbit at ', jd, sod, ' for ', OB%prn(isat)
          write (lfncid) 'de'
          write (lfnrem) jd, sod, isat, DELORB
        else if (OB%var(isat, 1) .eq. -20.d0) then
          write (*, '(a,i5,f10.2,a,a3)') '###WARNING(lsq): no clock at ', jd, sod, ' for ', OB%prn(isat)
          write (lfncid) 'de'
          write (lfnrem) jd, sod, isat, DELCLK
        end if
      else if (OB%elev(isat) .lt. SITE%cutoff) then
        OB%omc(isat, 1:4) = 0.d0
        write (*, '(a,i5,f10.2,a,a3,f6.2)') '###WARNING(lsq): low elev at ', jd, sod, ' for ', OB%prn(isat), OB%elev(isat)*180.d0/PI
        write (lfncid) 'de'
        write (lfnrem) jd, sod, isat, DELLOW
      end if
      if (OB%prn(isat)(1:1) .ne. 'R') then
        if (LCF%del_nobia) then
          if (nbias_used(isat) .eq. -1) cycle
          if (nbias_used(isat) .lt.  4) then
            OB%omc(isat, 1:4) = 0.d0
            OB%var(isat, 1:4) = 0.d0
            write (*, '(a,i5,f10.2,a,a3,i5)') '###WARNING(lsq): no bias at ', jd, sod, ' for ', OB%prn(isat), nbias_used(isat)
            write (lfncid) 'de'
            write (lfnrem) jd, sod, isat, DELBIA
          end if
        end if
        if (.not. arsig_used(isat)) then
          arsig_used(isat) = (nbias_used(isat) .ge. 4)
        end if
      end if
    end do
!
!! save a priori receiver clock correction
    ipar = pointer_string(OB%npar, OB%pname, 'RECCLK_G')
    if (ipar .ne. 0) PM(OB%ltog(ipar, index_g(1)))%xini = SITE%rclock_G
    ipar = pointer_string(OB%npar, OB%pname, 'RECCLK_R')
    if (ipar .ne. 0) PM(OB%ltog(ipar, index_r(1)))%xini = SITE%rclock_R
    ipar = pointer_string(OB%npar, OB%pname, 'RECCLK_E')
    if (ipar .ne. 0) PM(OB%ltog(ipar, index_e(1)))%xini = SITE%rclock_E
    ipar = pointer_string(OB%npar, OB%pname, 'RECCLK_C')
    if (ipar .ne. 0) PM(OB%ltog(ipar, index_c(1)))%xini = SITE%rclock_C
    ipar = pointer_string(OB%npar, OB%pname, 'RECCLK_3')
    if (ipar .ne. 0) PM(OB%ltog(ipar, index_3(1)))%xini = SITE%rclock_3
    ipar = pointer_string(OB%npar, OB%pname, 'RECCLK_J')
    if (ipar .ne. 0) PM(OB%ltog(ipar, index_j(1)))%xini = SITE%rclock_J
!
!! save a priori horizontal troposphere gradients
    if (LCF%htgmod(1:3) .ne. 'NON') then
      ipar = pointer_string(OB%npar, OB%pname, 'HTGC')
      PM(OB%ltog(ipar, 1))%xini = OB%nhtg
      ipar = pointer_string(OB%npar, OB%pname, 'HTGS')
      PM(OB%ltog(ipar, 1))%xini = OB%ehtg
    end if
!
!! add to normal equation
    call lsq_add_obs(lfncid, lfnobs, lfnrem, jd, sod, LCF, OB, PM, NM, SAT, SITE)
!
!! count epoch number
    do ipar = 1, OB%npar
      j = OB%ltog(ipar, 1)
      if (j .eq. 0 .or. j .gt. OB%npar) cycle
      if (PM(j)%iobs .eq. 0) cycle
      PM(j)%iepo = PM(j)%iepo + 1
    end do
!
!! save a priori zenith troposphere delay for each epoch
    if (LCF%ztdmod(1:3) .ne. 'NON') then
      ipar = pointer_string(OB%npar, OB%pname, 'ZTD')
      ipar = OB%ltog(ipar, 1)
      if (PM(ipar)%iobs .gt. 0) then
        write (lfncid) 'pz'
        write (lfnrem) ipar, OB%zdd, OB%zwd, jd + sod/864.d2
      end if
    end if
!
!! reset ambiguity flag
40  do isat = 1, LCF%nprn
      if (OB%omc(isat, 1) .eq. 0.d0) cycle
      iamb = pointer_string(OB%npar, OB%pname, 'AMBC')
      ipar = OB%ltog(iamb, isat)
      if (ipar .eq. 0) cycle
      if (OB%flag(isat) .eq. 1) then
        PM(ipar)%ptime(1) = jd + sod/864.d2
        OB%flag(isat) = 0
        OB%lifamb(isat, 1:2) = 0.d0
        NM%nuk = NM%nuk + 1
      end if
      PM(ipar)%ptime(2) = jd + sod/864.d2
    end do
!
!! calculate PDOP
    call lsq_dop(jd, sod, LCF, OB, pdop)
    pdop_sum = pdop_sum + pdop
!
!! stochastic process
    call lsq_process(lfncid, lfnrem, jd, sod, LCF, NM, PM, pdop_sum)
!
!! ambiguity
    if (LCF%lrmbias) then
!
!! find bias to be removed
      nbias = 0
      do ipar = NM%nc + NM%np + 1, NM%ipm
        if (PM(ipar)%pname(1:3) .ne. 'AMB' .or. PM(ipar)%ipt .eq. 0) cycle
        if ((jd - PM(ipar)%ptend)*864.d2 + sod .lt. MAXWND - LCF%dintv + 1.d-3) cycle
        isat = PM(ipar)%psat
        iamb = pointer_string(OB%npar, OB%pname, 'AMBC')
        nbias = nbias + 1
        ibias(nbias) = ipar
        OB%flag(isat) = 1
        OB%ltog(iamb, isat) = 0
        OB%lifamb(isat, 1:2) = 0.d0
      end do
!
!! remove bias
      if (nbias .ne. 0) then
        call lsq_add_ambcon(lfncid, lfnobs, jd, sod, LCF, SITE, NM, PM, OB)
        call lsq_rmv_normal(lfncid, lfnrem, nbias, ibias, NM, PM)
      end if
    end if
!
!! next epoch
    call timinc(jd, sod, LCF%dintv, jd, sod)
  end do
!
!! apply all remaining ambiguity constraints
  call lsq_add_ambcon(lfncid, lfnobs, jd, sod, LCF, SITE, NM, PM, OB)
!
!! remove all 'P' and 'S' parameters
  nbias = 0
  do ipar = 1, NM%ipm
    if (PM(ipar)%ptype .eq. 'C' .or. PM(ipar)%ipt .eq. 0) cycle
    if (PM(ipar)%ptype .eq. 'S' .and.  .not. LCF%lrmbias) cycle
    nbias = nbias + 1
    ibias(nbias) = ipar
  end do
  if (nbias .ne. 0) call lsq_rmv_normal(lfncid, lfnrem, nbias, ibias, NM, PM)
!
!! solve normal equation
  write (*, '(a,3i6)') 'Resolving : nc, np, ns ', NM%nc, NM%np, NM%ns
  call lsq_slv_prmt(LCF, NM, PM)
!
!! Recover parameters and compute residual
  call lsq_rcv_prmt(lfncid, lfnobs, lfnrem, lfnres, LCF, SITE, OB, NM, PM)

  close (lfncid, status='delete')
  close (lfnobs, status='delete')
  close (lfnrem, status='delete')
  inquire (file=LCF%flnpos, opened=lopen, number=openid)
  if (lopen) close (openid)
  if (SITE%skd(1:1) .eq. 'S' .or. SITE%skd(1:1) .eq. 'F') then
    lfnpos = get_valid_unit(10)
    open (lfnpos, file=LCF%flnpos)
    call lsq_wrt_header(lfnpos, LCF, SITE, OB, 'pos', .true., .true., .true., .true.)
    ipar = pointer_string(OB%npar, OB%pname, 'STAPX')
    ipar = OB%ltog(ipar, 1)
    if (PM(ipar)%iobs .ne. 0) then
      write (lfnpos, '(1x,a4,1x,f11.4,3f16.5,7e25.14,i15)') &
        SITE%name, LCF%jd0 + (LCF%sod0 + timdif(LCF%jd1, LCF%sod1, LCF%jd0, LCF%sod0)/2.d0)/864.d2, &
        PM(ipar)%xini + PM(ipar)%xcor, PM(ipar + 1)%xini + PM(ipar + 1)%xcor, PM(ipar + 2)%xini + PM(ipar + 2)%xcor, &
        NM%norx(1, 1), NM%norx(2, 2), NM%norx(3, 3), &
        NM%norx(1, 2), NM%norx(1, 3), NM%norx(2, 3), &
        NM%sig0, PM(ipar)%iobs
    end if
    close (lfnpos)
  end if
!
!! resolvable ambiguities
  LCF%fcbnprn = 0
  LCF%fcbprn = ' '
  do ipar = 1, NM%ipm
    if (PM(ipar)%ptype .eq. 'S' .and. PM(ipar)%iobs .gt. 0) then
      k = pointer_string(LCF%nprn, LCF%prn, LCF%prn(PM(ipar)%psat))
      if (k .eq. 0) cycle
      if (all(BIAS(k,   1:9)%length .eq. 0) .or. all(BIAS(k, 10:18)%length .eq. 0) .or. &
          all(BIAS(k, 19:27)%length .eq. 0) .or. all(BIAS(k, 28:36)%length .eq. 0)) cycle
      if (arsig_used(k) .and. pointer_string(LCF%fcbnprn, LCF%fcbprn, LCF%prn(PM(ipar)%psat)) .eq. 0) then
        LCF%fcbnprn = LCF%fcbnprn + 1
        LCF%fcbprn(LCF%fcbnprn) = LCF%prn(PM(ipar)%psat)
      end if
    end if
  end do
  do i = 1, LCF%fcbnprn - 1
    do j = i + 1, LCF%fcbnprn
      id = pointer_string(LCF%nprn, LCF%prn, LCF%fcbprn(i))
      ih = pointer_string(LCF%nprn, LCF%prn, LCF%fcbprn(j))
      if (id .gt. ih) then
        temprn = LCF%fcbprn(i)
        LCF%fcbprn(i) = LCF%fcbprn(j)
        LCF%fcbprn(j) = temprn
      end if
    end do
  end do
!
!! output ambiguities
  k = 0
  inquire (file=LCF%flnamb, opened=lopen, number=openid)
  if (lopen) close (openid)
  lfnamb = get_valid_unit(10)
  open (lfnamb, file=LCF%flnamb)
  call lsq_wrt_header(lfnamb, LCF, SITE, OB, 'amb', .true., .true., .true., .true.)
  do ipar = 1, NM%ipm
    if (PM(ipar)%ptype .eq. 'S') then
      if (PM(ipar)%iobs .gt. 0) then
        PM(ipar)%xrms = dsqrt(PM(ipar)%xrms/PM(ipar)%iobs)
        PM(ipar)%mele = PM(ipar)%mele/PM(ipar)%iobs*180.d0/PI
        if (PM(ipar)%iobs .gt. 1) then
          PM(ipar)%xswl = dsqrt(PM(ipar)%xswl/PM(ipar)%rw - &
                               (PM(ipar)%xrwl/PM(ipar)%rw)**2)/dsqrt(PM(ipar)%iobs - 1.d0)
        else
          PM(ipar)%xswl = 999.9999d0
        end if
        PM(ipar)%xrwl = PM(ipar)%xrwl/PM(ipar)%rw + PM(ipar)%zw
        k = k + 1
      end if
      if (LCF%prn(PM(ipar)%psat)(1:1) .ne. 'R') then
        i0 = index(GNSS_PRIO, LCF%prn(PM(ipar)%psat)(1:1))
        freq = FREQ_SYS(idxfrq(i0, 1), i0)
      else
        read (LCF%prn(PM(ipar)%psat), '(1x,i2)') prn_int
        frequency_glo_nu = OB%glschn(prn_int)
        freq = FREQ1_R(frequency_glo_nu)
      end if
      write (lfnamb, '(a4,2f22.6,2f18.10,2f9.4,f6.1)') &
        LCF%prn(PM(ipar)%psat), (PM(ipar)%xini + PM(ipar)%xcor)*freq/VLIGHT, PM(ipar)%xrwl, &
        PM(ipar)%ptime(1:2), PM(ipar)%xrms*freq/VLIGHT, PM(ipar)%xswl, PM(ipar)%mele
    end if
  end do
  close (lfnamb)
!
!! output inversed normal matrix
  if (.not. LCF%lrmbias) then
    lfnneq = get_valid_unit(10)
    open (lfnneq, file=LCF%flnneq, form='unformatted')
    write (lfnneq) SITE%name
    write (lfnneq) LCF%fcbnprn, (LCF%fcbprn(i), i=1, LCF%fcbnprn)
    write (lfnneq) LCF%nprn, (LCF%prn(i), i=1, LCF%nprn)
    write (lfnneq) NM%imtx, NM%ltpl, NM%nobs - NM%nuk
    do i = 1, NM%imtx
      j = NM%iptp(i)
      write (lfnneq) PM(j)%pname, PM(j)%psat, PM(j)%xrwl, &
        PM(j)%xini + PM(j)%xcor, PM(j)%ptime(1:2), PM(j)%xswl, PM(j)%mele
    end do
    write (lfnneq) ((NM%norx(i, j), i=j, NM%imtx), j=1, NM%imtx)
    close (lfnneq)
  end if
!
!! clean
  if (SITE%lfnjmp .ne. 0) close (SITE%lfnjmp)
  deallocate (PM)
  deallocate (NM%norx)
  deallocate (NM%iptx)
  deallocate (NM%iptp)
  call clean_gpt3_1()
  if (LCF%lioh) call clean_rdionex(IM)
end program
