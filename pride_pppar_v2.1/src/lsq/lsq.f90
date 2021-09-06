!
!! lsq.f90
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
!! Least Squares Estimator
!
program lsq
  implicit none
  include '../header/const.h'
  include '../header/orbit.h'
  include '../header/rnxobs.h'
  include '../header/station.h'
  include '../header/satellite.h'
  include '../header/ionex.h'
  include 'lsqcfg.h'
  include 'lsq.h'
!
!! lsq configure
  type(lsqcfg) LCF
!! station
  type(station) SITE
  type(rnxhdr) HD
  type(rnxobr) OB
!! orbit
  type(orbhdr) OH
  type(satellite) SAT(MAXSAT)
!! normal matrix
  type(norm) NM
!! parameter definition
  type(prmt), pointer :: PM(:)
!
!! ionosphere
  type(ionex) IM
!!
  logical*1 lopen 
  integer*4 i, j, k, jd, jd_sav, jdc, isat, iepo, ipar, iamb, id, ih, openid, ierr
  integer*4 lfncid, lfnobs, lfnrem, lfnres, lfnpos, lfnamb, lfnneq, nbias, ibias(MAXPAR)
  real*8 dwnd, sod, sec, sodc, deltax(3), xsat(6), bias(MAXSAT, 44)
  character*1 cin
  character*3 temprn,prn_mat(MAXSAT)
  character*20 antnum
!
!! function called
  integer*4 get_valid_unit, pointer_string
  real*8 timdif
  character*10 run_tim
! 
  real*8 PDOP,TDOP,GDOP,nbb(4,4)
! R
  integer*4 frequency_glo_nu,prn_int
  real*8 :: FREQ1_R(-50:50),FREQ2_R(-50:50)

  call prn_matbld(prn_mat)
  call frequency_glonass(FREQ1_R,FREQ2_R)
!
!! initial
  do i = 1, MAXSAT
    LCF%prn(i) = ''
  enddo
  LCF%attuse = .false.
  LCF%otluse = .false.
!
!! get arguments
  call get_lsq_args(LCF, SITE, OB, SAT, IM)
  
  dwnd = min(LCF%dintv/10.d0, 0.3d0)
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
    write(antnum,'(a3)') SAT(isat)%prn
    call get_ant_ipt(LCF%jd0 + LCF%sod0/86400.d0, LCF%jd1 + LCF%sod1/86400.d0, &
                     SAT(isat)%typ, antnum, SAT(isat)%iptatx, SAT(isat)%xyz, SAT(isat)%prn(1:1), LCF%flnatx_real)
  enddo
!
!! site antenna
!! observation
  SITE%iunit = get_valid_unit(10)
  open (SITE%iunit, file=SITE%obsfil, status='old', iostat=ierr)
  if (ierr .ne. 0) then
    write (*, '(a)') '***ERROR(lsq): open file ', trim(SITE%obsfil)
    call exit(1)
  endif
  call rdrnxoh(SITE%iunit, HD, ierr)
  if (ierr .ne. 0) then
    write (*, '(2a)') '***ERROR(lsq): read header ', trim(SITE%obsfil)
    call exit(1)
  endif
  SITE%anttyp = HD%anttyp
  SITE%enu0(1) = HD%e
  SITE%enu0(2) = HD%n
  SITE%enu0(3) = HD%h
  antnum = ' '
  call get_ant_ipt(LCF%jd0 + LCF%sod0/86400.d0, LCF%jd1 + LCF%sod1/86400.d0, &
                   SITE%anttyp, antnum, SITE%iptatx, SITE%enu_G, 'G', LCF%flnatx_real)
  call get_ant_ipt(LCF%jd0 + LCF%sod0/86400.d0, LCF%jd1 + LCF%sod1/86400.d0, &
                   SITE%anttyp, antnum, SITE%iptatx, SITE%enu_R, 'R', LCF%flnatx_real)
  call get_ant_ipt(LCF%jd0 + LCF%sod0/86400.d0, LCF%jd1 + LCF%sod1/86400.d0, &
                   SITE%anttyp, antnum, SITE%iptatx, SITE%enu_E, 'E', LCF%flnatx_real)
  call get_ant_ipt(LCF%jd0 + LCF%sod0/86400.d0, LCF%jd1 + LCF%sod1/86400.d0, &
                   SITE%anttyp, antnum, SITE%iptatx, SITE%enu_C, 'C', LCF%flnatx_real)
  call get_ant_ipt(LCF%jd0 + LCF%sod0/86400.d0, LCF%jd1 + LCF%sod1/86400.d0, &
                   SITE%anttyp, antnum, SITE%iptatx, SITE%enu_J, 'J', LCF%flnatx_real)
!
!! count number of parameters
  call lsq_cnt_prmt(LCF, SITE, NM)
  NM%npm = 0
  NM%npm = NM%npm + OB%amb_tot
  NM%npm = NM%npm + NM%imtx
  if (NM%npm .gt. MAXPAR) then
    write (*, '(a,i8)') '***ERROR(lsq): too many parameters ', NM%npm
    call exit(1)
  endif
  allocate (PM(NM%npm), stat=ierr)
  if (ierr .ne. 0) then
    write (*, '(a,i8)') '***ERROR(lsq): parameter array allocation ', NM%npm
    call exit(1)
  endif
!
!! normal matrix size
  if (LCF%lrmbias) then
    NM%nmtx = 0
    NM%nmtx = NM%nmtx + OB%amb_epo
    NM%nmtx = NM%nmtx + NM%imtx + 1
  else
    NM%nmtx = NM%npm + 1
  endif
!! allocate memory
  allocate (NM%norx(1:NM%nmtx, 1:NM%nmtx), stat=ierr)
  if (ierr .ne. 0) then
    write (*, '(a,i8)') '***ERROR(lsq): normal matrix allocation ', NM%nmtx
    call exit(1)
  endif
  allocate (NM%iptp(NM%nmtx), stat=ierr)
  if (ierr .ne. 0) then
    write (*, '(a,i8)') '***ERROR(lsq): iptp allocation ', NM%nmtx
    call exit(1)
  endif
  allocate (NM%iptx(NM%nmtx), stat=ierr)
  if (ierr .ne. 0) then
    write (*, '(a,i8)') '***ERROR(lsq): iptx allocation ', NM%nmtx
    call exit(1)
  endif
  do i = 1, NM%nmtx
    NM%iptp(i) = 0
    NM%iptx(i) = (i - 1)*NM%nmtx
  enddo
  write (*, '(/,3(a,i8,/))') ' Size of Normal Matrix :', NM%nmtx, &
    ' Size of NC Parameter  :', NM%nc, &
    ' Size of NP Parameter  :', NM%np
!
!! initiate normal matrix
  call lsq_init(LCF, SITE, SAT, OB, NM, PM)
!
!!++++++++++++++++++++EPOCH LOOP+++++++++++++++++++++++++++++
  jd = LCF%jd0
  sod = LCF%sod0
  iepo = 0
  NM%ltpl = 0.d0
  NM%nuk = NM%nc
  NM%nobs = 0
  jd_sav = 0
  OB%typuse = ' '
  do while (timdif(jd, sod, LCF%jd1, LCF%sod1) .lt. MAXWND)
!
!! read pppar biases
    if (jd .gt. jd_sav) then
      bias = 0.d0
      if(LCF%fcbuse) call read_bias(LCF%flnfcb, bias, jd*1.d0, (jd+1)*1.d0)
      jd_sav = jd
    endif
!
!! read one epoch of observations
    k = 0
    if (SITE%iunit .eq. 0) cycle
    if (HD%ver .eq. 3) then
      call rdrnxoi3(SITE%iunit, jd, sod, dwnd, LCF%nprn, LCF%prn, HD, OB, bias, ierr)
    else
      call rdrnxoi(SITE%iunit, jd, sod, dwnd, LCF%nprn, LCF%prn, HD, OB, bias, ierr)
    endif
    if (ierr .ne. 0) SITE%iunit = 0
    if (ierr .eq. 0) then
      k = k + 1
      call read_obsrhd(jd, sod, LCF%nprn, LCF%prn, OB)
      if(SITE%lfnjmp .ne. 0) call read_clkjmp(jd, sod, SITE%lfnjmp, LCF%nprn, OB)
    endif
    if (k .eq. 0) then
      write (*, '(a)') '###WARNING(lsq): no more data to be processed'
      exit
    endif
!
!! +++++++++++++++++++++++++++SITE LOOP+++++++++++++++++++++++++++++++++++
!! add observation equations
!
!! read kinematic position: if not found, use value at last epoch
    if (SITE%skd(1:1) .eq. 'K' .and. count(OB%obs(1:LCF%nprn, 3) .ne. 0.d0) .gt. 0) then
      call read_kinpos(SITE, jd, sod, deltax(1), deltax(2), deltax(3))
      if (.not. all(deltax(1:3) .eq. 1.d0)) then
        ipar = pointer_string(OB%npar, OB%pname, 'STAPX')
        ipar = OB%ltog(ipar, 1)
        PM(ipar)%xini     = deltax(1)
        PM(ipar + 1)%xini = deltax(2)
        PM(ipar + 2)%xini = deltax(3)
      else
        goto 40
      endif
    endif
!
!! add new ambiguities
    call lsq_add_newamb(jd, sod, SITE%name, LCF, OB, NM, PM)
!
!! iteration for more precise position
    do isat = 1, LCF%nprn
      OB%omc(isat, 1:4) = 0.d0
    enddo
    if (count(OB%obs(1:LCF%nprn, 3) .ne. 0.d0) .eq. 0) goto 40
    iepo = iepo + 1
    if (dmod(iepo*1.d0, dnint(300.d0/LCF%dintv)) .le. MAXWND) write (*, '(2a,i9,i7,f9.1)') run_tim(), ' Epoch ', iepo, jd, sod
!
!! prepare a priori information for kinematic station
    if (SITE%skd(1:1) .eq. 'K') then
      ipar = pointer_string(OB%npar, OB%pname, 'STAPX')
      do i = 0, 2
        SITE%x(i + 1) = PM(OB%ltog(ipar + i, 1))%xini*1.d-3
      enddo
      call xyzblh(SITE%x(1:3)*1.d3, 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, SITE%geod)
      SITE%geod(3) = SITE%geod(3)*1.d-3
      call rot_enu2xyz(SITE%geod(1), SITE%geod(2), SITE%rot_l2f)
    endif
!
!! prepare a priori receiver clock correction (unit: m)
    ipar = pointer_string(OB%npar, OB%pname, 'RECCLK_G')
    if (ipar .ne. 0) SITE%rclock_G = PM(OB%ltog(ipar, 1))%xini
    ipar = pointer_string(OB%npar, OB%pname, 'RECCLK_R')
    if (ipar .ne. 0) SITE%rclock_R = PM(OB%ltog(ipar, 1))%xini
    ipar = pointer_string(OB%npar, OB%pname, 'RECCLK_E')
    if (ipar .ne. 0) SITE%rclock_E = PM(OB%ltog(ipar, 1))%xini
    ipar = pointer_string(OB%npar, OB%pname, 'RECCLK_C')
    if (ipar .ne. 0) SITE%rclock_C = PM(OB%ltog(ipar, 1))%xini
    ipar = pointer_string(OB%npar, OB%pname, 'RECCLK_3')
    if (ipar .ne. 0) SITE%rclock_3 = PM(OB%ltog(ipar, 1))%xini
    ipar = pointer_string(OB%npar, OB%pname, 'RECCLK_J')
    if (ipar .ne. 0) SITE%rclock_J = PM(OB%ltog(ipar, 1))%xini
!
!! one way model
    call read_meteo(jd, sod, SITE%imet, SITE%map, SITE%geod, SITE%p0, SITE%t0, SITE%hr0, SITE%undu)
    OB%zdd  = 0.d0
    OB%zwd  = 0.d0
    OB%nhtg = 0.d0
    OB%ehtg = 0.d0
    call troposphere_delay(jd, sod, SITE, OB%zdd, OB%zwd)
    if(LCF%sys(1:1) .eq. 'G')call gpsmod(jd, sod, LCF, SITE, OB, SAT, IM)
    if(LCF%sys(2:2) .eq. 'R')call glomod(jd, sod, LCF, SITE, OB, SAT, IM)
    if(LCF%sys(3:3) .eq. 'E')call galmod(jd, sod, LCF, SITE, OB, SAT, IM)
    if(LCF%sys(4:4) .eq. 'C')call bd2mod(jd, sod, LCF, SITE, OB, SAT, IM)
    if(LCF%sys(5:5) .eq. '3')call bd3mod(jd, sod, LCF, SITE, OB, SAT, IM)
    if(LCF%sys(6:6) .eq. 'J')call qzsmod(jd, sod, LCF, SITE, OB, SAT, IM)
!
!! check elevation & cutoff angle
    do isat = 1, LCF%nprn
      if (OB%omc(isat, 1) .eq. 0.d0 .or. OB%omc(isat, 3) .eq. 0.d0) then
        if(OB%obs(isat, 1) .ne. 0.d0 .and. OB%obs(isat, 3) .ne. 0.d0) then
          if(OB%var(isat,1) .eq. -10.d0) then
            write (*, '(a,i5,f10.2,a,a3)') '###WARNING(lsq): no orbit at ', jd, sod, &
                                          ' for ', OB%prn(isat)
            write(lfncid) 'de'
            write(lfnrem) jd, sod, isat, 7
          else if(OB%var(isat, 1) .eq. -20.d0) then
            write (*, '(a,i5,f10.2,a,a3)') '###WARNING(lsq): no clock at ', jd, sod, &
                                          ' for ', OB%prn(isat)
            write(lfncid) 'de'
            write(lfnrem) jd, sod, isat, 8
          endif
        endif
      else if (OB%elev(isat) .lt. SITE%cutoff) then
        OB%omc(isat, 1:4) = 0.d0
        write (*, '(a,i5,f10.2,a,a3,f6.2)') '###WARNING(lsq): low elev at ', jd, sod, &
                                            ' for ', OB%prn(isat), OB%elev(isat)*180.d0/PI
        write (lfncid) 'de'
        write (lfnrem) jd, sod, isat, 4
      endif
    enddo
!
!! save a priori receiver clock correction
    ipar = pointer_string(OB%npar, OB%pname, 'RECCLK_G')
    if (ipar .ne. 0) PM(OB%ltog(ipar, 1))%xini = SITE%rclock_G
    ipar = pointer_string(OB%npar, OB%pname, 'RECCLK_R')
    if (ipar .ne. 0) PM(OB%ltog(ipar, 1))%xini = SITE%rclock_R
    ipar = pointer_string(OB%npar, OB%pname, 'RECCLK_E')
    if (ipar .ne. 0) PM(OB%ltog(ipar, 1))%xini = SITE%rclock_E
    ipar = pointer_string(OB%npar, OB%pname, 'RECCLK_C')
    if (ipar .ne. 0) PM(OB%ltog(ipar, 1))%xini = SITE%rclock_C
    ipar = pointer_string(OB%npar, OB%pname, 'RECCLK_3')
    if (ipar .ne. 0) PM(OB%ltog(ipar, 1))%xini = SITE%rclock_3
    ipar = pointer_string(OB%npar, OB%pname, 'RECCLK_J')
    if (ipar .ne. 0) PM(OB%ltog(ipar, 1))%xini = SITE%rclock_J
!
!! save a priori horizontal troposphere gradients
    if (LCF%htgmod(1:3) .ne. 'NON' .and. LCF%htgmod(1:3) .ne. 'FIX') then
      ipar = pointer_string(OB%npar, OB%pname, 'HTGC'//trim(LCF%htgmod))
      PM(OB%ltog(ipar, 1))%xini = OB%nhtg
      ipar = pointer_string(OB%npar, OB%pname, 'HTGS'//trim(LCF%htgmod))
      PM(OB%ltog(ipar, 1))%xini = OB%ehtg
    endif
!
!! add to normal equation
    call lsq_add_obs(lfncid, lfnobs, lfnrem, jd, sod, LCF, OB, PM, NM, SAT, SITE)
    !!add for DOP
    !
    call lsq_dop(jd,sod,LCF,OB,nbb)
    PDOP=dsqrt(nbb(1,1)+nbb(2,2)+nbb(3,3))
    TDOP=dsqrt(nbb(4,4))
    GDOP=dsqrt(nbb(1,1)+nbb(2,2)+nbb(3,3)+nbb(4,4))
    PM%pdop=PDOP
!
!! save a priori zenith troposphere delay for each epoch
    if (LCF%ztdmod(1:3) .ne. 'NON' .and. LCF%ztdmod(1:3) .ne. 'FIX') then
      ipar = pointer_string(OB%npar, OB%pname, 'ZTD'//trim(LCF%ztdmod))
      ipar = OB%ltog(ipar, 1)
      if (PM(ipar)%iobs .gt. 0) then
        write (lfncid) 'pz'
        write (lfnrem) ipar, OB%zdd, OB%zwd, jd + sod/86400.d0
      endif
    endif
!
!! reset ambiguity flag

40  do isat = 1, LCF%nprn
      if (OB%omc(isat, 1) .eq. 0.d0) cycle
      iamb = pointer_string(OB%npar, OB%pname, 'AMBC')
      ipar = OB%ltog(iamb, isat)
      if (ipar .eq. 0) cycle
      if (OB%flag(isat) .eq. 1) then
        PM(ipar)%ptime(1) = jd + sod/86400.d0
        OB%flag(isat) = 0
        OB%lifamb(isat, 1:2) = 0.d0
        NM%nuk = NM%nuk + 1
      endif
      PM(ipar)%ptime(2) = jd + sod/86400.d0
    enddo

!
!! stochastic process
    call lsq_process(lfncid, lfnrem, jd, sod, LCF%dintv, NM, PM)
!
!! ambiguity
    if (LCF%lrmbias) then
!
!! find bias to be removed
      nbias = 0
      do ipar = NM%nc + NM%np + 1, NM%ipm
        if (PM(ipar)%pname(1:3) .ne. 'AMB' .or. PM(ipar)%ipt .eq. 0) cycle
        if ((jd - PM(ipar)%ptend)*86400.d0 + sod .lt. MAXWND - LCF%dintv) cycle
        isat = PM(ipar)%psat
        iamb = pointer_string(OB%npar, OB%pname, 'AMBC')
        nbias = nbias + 1
        ibias(nbias) = ipar
        OB%flag(isat) = 1
        OB%ltog(iamb, isat) = 0
        OB%lifamb(isat, 1:2) = 0.d0
      enddo
!
!! remove bias
      if (nbias .ne. 0) then
        call lsq_add_ambcon(jd, sod, LCF, SITE, NM, PM, OB)
        call lsq_rmv_normal(lfncid, lfnrem, nbias, ibias, NM, PM)
      endif
    endif
!
!! next epoch
    call timinc(jd, sod, LCF%dintv, jd, sod)
  enddo
!
!! apply all remaining ambiguity constraints
  call lsq_add_ambcon(jd, sod, LCF, SITE, NM, PM, OB)
!
!! remove all 'P' and 'S' parameters
  nbias = 0
  do ipar = 1, NM%ipm
    if (PM(ipar)%ptype .eq. 'C' .or. PM(ipar)%ipt .eq. 0) cycle
    if (PM(ipar)%ptype .eq. 'S' .and. .not. LCF%lrmbias) cycle
    nbias = nbias + 1
    ibias(nbias) = ipar
  enddo
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
    cin = ''
    if (PM(ipar)%iobs .eq. 0) cin='*'
    write (lfnpos, '(1x,a4,a1,f11.4,3f15.4,7e25.14,i15)') &
      SITE%name, cin, LCF%jd0 + (LCF%sod0 + timdif(LCF%jd1, LCF%sod1, LCF%jd0, LCF%sod0)/2.d0)/86400.d0, &
      PM(ipar)%xini+PM(ipar)%xcor, PM(ipar + 1)%xini+PM(ipar + 1)%xcor, PM(ipar + 2)%xini+PM(ipar + 2)%xcor, &
      NM%norx(1, 1), NM%norx(2, 2), NM%norx(3, 3), &
      NM%norx(1, 2), NM%norx(1, 3), NM%norx(2, 3), &
      NM%sig0, PM(ipar)%iobs
    close (lfnpos)
  endif
!
!! resolvable ambiguities
  LCF%fcbnprn=0
  LCF%fcbprn=' '
  do ipar=1,NM%ipm
    if(PM(ipar)%ptype.eq.'S' .and. PM(ipar)%iobs.gt.0) then
      k=pointer_string(MAXSAT, prn_mat, LCF%prn(PM(ipar)%psat))
      if(k.ne.0 .and. any(bias(k,1:11).ne.0.d0) .and. any(bias(k,12:22).ne.0.d0)) then
        if(pointer_string(LCF%fcbnprn,LCF%fcbprn,LCF%prn(PM(ipar)%psat)).eq.0) then
          LCF%fcbnprn=LCF%fcbnprn+1
          LCF%fcbprn(LCF%fcbnprn)=LCF%prn(PM(ipar)%psat)
        endif
      endif
    endif
  enddo
  do i=1,LCF%fcbnprn-1
    do j=i+1,LCF%fcbnprn
      id=pointer_string(MAXSAT, prn_mat, LCF%fcbprn(i))
      ih=pointer_string(MAXSAT, prn_mat, LCF%fcbprn(j))
      if(id .gt. ih) then
        temprn=LCF%fcbprn(i)
        LCF%fcbprn(i)=LCF%fcbprn(j)
        LCF%fcbprn(j)=temprn
      endif
    enddo
  enddo
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
        endif
        PM(ipar)%xrwl = PM(ipar)%xrwl/PM(ipar)%rw + PM(ipar)%zw
        k = k + 1
      endif
      if(LCF%prn(PM(ipar)%psat)(1:1).eq.'G')then
        write (lfnamb, '(a4,2f22.6,2f18.10,2f9.4,f6.1)') &
          LCF%prn(PM(ipar)%psat), (PM(ipar)%xini + PM(ipar)%xcor)*FREQ1_G/VLIGHT, PM(ipar)%xrwl, &
          PM(ipar)%ptime(1:2), PM(ipar)%xrms*FREQ1_G/VLIGHT, PM(ipar)%xswl, PM(ipar)%mele
      elseif(LCF%prn(PM(ipar)%psat)(1:1).eq.'R')then
        read(LCF%prn(PM(ipar)%psat),'(1x,i2)') prn_int
        frequency_glo_nu=OB%glschn(prn_int)
        write (lfnamb, '(a4,2f22.6,2f18.10,2f9.4,f6.1)') &
          LCF%prn(PM(ipar)%psat), (PM(ipar)%xini + PM(ipar)%xcor)*FREQ1_R(frequency_glo_nu)/VLIGHT, PM(ipar)%xrwl, &
          PM(ipar)%ptime(1:2), PM(ipar)%xrms*FREQ1_R(frequency_glo_nu)/VLIGHT, PM(ipar)%xswl, PM(ipar)%mele
      elseif(LCF%prn(PM(ipar)%psat)(1:1).eq.'E')then
        write (lfnamb, '(a4,2f22.6,2f18.10,2f9.4,f6.1)') &
          LCF%prn(PM(ipar)%psat), (PM(ipar)%xini + PM(ipar)%xcor)*FREQ1_E/VLIGHT, PM(ipar)%xrwl, &
          PM(ipar)%ptime(1:2), PM(ipar)%xrms*FREQ1_E/VLIGHT, PM(ipar)%xswl, PM(ipar)%mele
      elseif(LCF%prn(PM(ipar)%psat)(1:1).eq.'C')then
        write (lfnamb, '(a4,2f22.6,2f18.10,2f9.4,f6.1)') &
          LCF%prn(PM(ipar)%psat), (PM(ipar)%xini + PM(ipar)%xcor)*FREQ1_C/VLIGHT, PM(ipar)%xrwl, &
          PM(ipar)%ptime(1:2), PM(ipar)%xrms*FREQ1_C/VLIGHT, PM(ipar)%xswl, PM(ipar)%mele
      elseif(LCF%prn(PM(ipar)%psat)(1:1).eq.'J')then
        write (lfnamb, '(a4,2f22.6,2f18.10,2f9.4,f6.1)') &
          LCF%prn(PM(ipar)%psat), (PM(ipar)%xini + PM(ipar)%xcor)*FREQ1_J/VLIGHT, PM(ipar)%xrwl, &
          PM(ipar)%ptime(1:2), PM(ipar)%xrms*FREQ1_J/VLIGHT, PM(ipar)%xswl, PM(ipar)%mele
      endif
    endif
  enddo
  close (lfnamb)
!
!! output inversed normal matrix
  if (.not. LCF%lrmbias) then
    lfnneq = get_valid_unit(10)
    open (lfnneq, file=LCF%flnneq, form='unformatted')
    write (lfnneq) LCF%fcbnprn, (LCF%fcbprn(i), i=1, LCF%fcbnprn)
    write (lfnneq) SITE%name
    write (lfnneq) LCF%nprn, (LCF%prn(i), i=1, LCF%nprn)
    write (lfnneq) NM%imtx, NM%ltpl, NM%nobs - NM%nuk
    do i = 1, NM%imtx
      j = NM%iptp(i)
      write (lfnneq) PM(j)%pname, PM(j)%psat, PM(j)%xrwl, &
        PM(j)%xini + PM(j)%xcor, PM(j)%ptime(1:2), PM(j)%xswl, PM(j)%mele
    enddo
    write (lfnneq) ((NM%norx(i, j), i=j, NM%imtx), j=1, NM%imtx)
    close (lfnneq)
  endif
!
!! clean
  if (SITE%lfnjmp .ne. 0) close(SITE%lfnjmp)
  deallocate (PM)
  deallocate (NM%norx)
  deallocate (NM%iptx)
  deallocate (NM%iptp)
  call clean_gpt3_1()
  if(LCF%lioh) call clean_rdionex(IM)
!
!! End of lsq
end program
