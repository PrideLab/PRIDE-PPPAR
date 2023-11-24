!
!! lsq_rcv_prmt.f90
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
!! Contributor: Maorong Ge, Jianghui Geng, Songfeng Yang, Jihang Lin, Wenyi Li
!!
!!
!!
!! purpose   : recover pre-eliminated parameters and residual
!! parameter :
!!    input  : lfnres               -- residual file
!!             LCF                  -- least squares control
!!             SITE                 -- station control
!!             NM,PM                -- normal equation & parameters
!
subroutine lsq_rcv_prmt(lfncid, lfnobs, lfnrem, lfnres, LCF, SITE, OB, NM, PM)
  implicit none
  include '../header/const.h'
  include '../header/station.h'
  include '../header/rnxobs.h'
  include 'lsqcfg.h'
  include 'lsq.h'

!
!! common
  integer*4 idxfrq(MAXSYS, 2)
  common idxfrq
!
!! parameter
  integer*4 lfncid, lfnobs, lfnrem, lfnres
  type(lsqcfg) LCF
  type(station) SITE
  type(rnxobr) OB
  type(norm) NM
  type(prmt) PM(1:*)
!
!! local
  character*2 cid
  integer*2 iflg
  integer*4 i, j, k, ig, ik, ip(2), ipar, isat, ntot, nprn, ipt(0:MAXPAR), jd, jdw, jdx, iunit, ierr
  integer*4 lfnrck, lfnztd, lfnhtg, jdr, iy, imon, id, ih, im, idoy
  real*8 cof(MAXPAR), x(3), xc(6), sodr, phase, range, wphs, wrng, mw, t0, t1, elev, azim, dmap, wmap, pdop
  real*8 zdd, zwd, xini, xcor, dummy, rckref, sodw, sodx
!
  character*60 line
  character*3 temprn, prn(MAXSAT), typuse(4)
  real*8 rck(1:6)
!
!! function called
  integer*4 get_valid_unit, pointer_string
  real*8 timdif

  lfnres = get_valid_unit(10)
  open (lfnres, file=LCF%flnres)
  lfnrck = get_valid_unit(10)
  open (lfnrck, file=LCF%flnrck)
  if (LCF%ztdmod(1:3) .ne. 'NON') then
    lfnztd = get_valid_unit(10)
    open (lfnztd, file=LCF%flnztd)
  end if
  if (LCF%htgmod(1:3) .ne. 'NON') then
    lfnhtg = get_valid_unit(10)
    open (lfnhtg, file=LCF%flnhtg)
  end if
  if (SITE%skd(1:1) .eq. 'P' .or. SITE%skd(1:1) .eq. 'K' .or. SITE%skd(1:1) .eq. 'L') then
    if (SITE%ikin .ne. 0) close (SITE%ikin)
    SITE%ikin = get_valid_unit(10)
    open (SITE%ikin, file=SITE%kinfil)
  end if
!
!! recover parameters and residuals
  nprn = 0
  pdop = 0.d0
  jdr = 0
  sodr = 0.d0
  rck = 0.d0
  rckref = 0.d0
  do while (.true.)
    backspace lfncid
    read (lfncid, err=100) cid
    if (cid .eq. '00') exit
    if (cid .eq. 'ob') then
      backspace lfnobs
      read (lfnobs) jd, mw, isat, ntot, ipt(0), (ipt(k), cof(k), k=1, ntot), &
        phase, range, wphs, wrng, iflg, elev, azim, dmap, wmap, typuse(1:4)
      dummy = 0.d0
      do i = 1, ntot
        dummy = dummy + PM(ipt(i))%xcor*cof(i)
      end do
      if (ipt(0) .eq. 0) then
        range = range - dummy
      else
        phase = phase - dummy
        range = range - dummy + PM(ipt(ipt(0)))%xcor*cof(ipt(0))
        PM(ipt(ipt(0)))%xrms = PM(ipt(ipt(0)))%xrms + phase**2
      end if
      if (jdr .eq. 0) then
        jdr = jd
        sodr = mw
        LCF%jd_end = jd
        LCF%sod_end = mw
      else if ((jdr - jd)*86400.d0 + sodr - mw .gt. MAXWND) then
        call mjd2date(jdr, sodr, iy, imon, id, ih, im, dummy)
        write (lfnres, '(a3,i5,4i3,f11.7,i7,f10.2)') 'TIM', iy, imon, id, ih, im, dummy, jdr, sodr
        jdr = jd
        sodr = mw
      end if
      elev = elev/PI*180.d0
      azim = azim/PI*180.d0
      write (lfnres, '(a3,2f10.4,2d16.8,i3,f8.3,f9.3,4(1x,a3))') &
        LCF%prn(isat), phase, range, wphs, wrng, iflg, elev, azim, typuse(1:4)
      if (pointer_string(nprn, prn, LCF%prn(isat)) .eq. 0) then
        nprn = nprn + 1
        prn(nprn) = LCF%prn(isat)
      end if
      backspace lfnobs
    else if (cid .eq. 'cn') then
      backspace lfnobs
      read (lfnobs) jdw, sodw, jdx, sodx, i, j, ip(1:2), phase, wphs
      if (jdr .ne. 0 .and. (jdr - jdx)*86400.d0 + sodr - sodx .gt. MAXWND) then
        call mjd2date(jdr, sodr, iy, imon, id, ih, im, dummy)
        write (lfnres, '(a3,i5,4i3,f11.7,i7,f10.2)') 'TIM', iy, imon, id, ih, im, dummy, jdr, sodr
        jdr = jdx
        sodr = sodx
      end if
      phase = phase - PM(ip(1))%xcor + PM(ip(2))%xcor
      write (lfnres, '(a3,f20.15,d16.8,2x,2(i7,f10.2),2(1x,a3))') 'CST', phase, wphs, jdw, sodw, jdx, sodx, LCF%prn(i), LCF%prn(j)
      backspace lfnobs
    else if (cid .eq. 'de') then
      backspace lfnrem
      read (lfnrem) jd, mw, isat, iflg
      if (jdr .eq. 0) then
        jdr = jd
        sodr = mw
        LCF%jd_end = jd
        LCF%sod_end = mw
      else if ((jdr - jd)*86400.d0 + sodr - mw .gt. MAXWND) then
        call mjd2date(jdr, sodr, iy, imon, id, ih, im, dummy)
        write (lfnres, '(a3,i5,4i3,f11.7,i7,f10.2)') 'TIM', iy, imon, id, ih, im, dummy, jdr, sodr
        jdr = jd
        sodr = mw
      end if
      write (lfnres, '(a3,2f10.4,2d16.8,i3,f8.3,f9.3)') LCF%prn(isat), 0.d0, 0.d0, 0.d0, 0.d0, iflg, 0.d0, 0.d0
      if (pointer_string(nprn, prn, LCF%prn(isat)) .eq. 0) then
        nprn = nprn + 1
        prn(nprn) = LCF%prn(isat)
      end if
      backspace lfnrem
    else if (cid .eq. 'pc') then
      backspace lfnrem
      read (lfnrem) ipar, ntot, ipt(0), (ipt(k), k=1, ntot - 1), (cof(k), k=1, ntot), &
        PM(ipar)%iepo, PM(ipar)%iobs, PM(ipar)%iobs_G, PM(ipar)%iobs_R, PM(ipar)%iobs_E, &
        PM(ipar)%iobs_C, PM(ipar)%iobs_3, PM(ipar)%iobs_J, &
        PM(ipar)%xini, PM(ipar)%ptime(1:2)
      dummy = cof(ipt(0))
      cof(ipt(0)) = -PM(ipar)%map*PM(ipar)%rw**2
      if (PM(ipar)%ipt .lt. 0) then
        PM(ipar)%ipt = 0
        cof(ipt(0)) = 0.d0
        NM%np = NM%np + 1
      end if
      do i = 1, ntot - 1
        cof(ntot) = cof(ntot) - cof(i)*PM(ipt(i))%xcor
      end do
      PM(ipar)%xcor = cof(ntot)/dummy
      if (PM(ipar)%iobs .gt. 0) then
        mw = PM(ipar)%ptime(1)
        if (PM(ipar)%pname(1:6) .eq. 'RECCLK') then
          if (rckref .ne. 0.d0 .and. (rckref - mw)*86400.d0 .gt. MAXWND) then
            call mjd2date(0, rckref, iy, imon, id, ih, im, dummy)
            if (dabs(dummy - 60.d0) .lt. MAXWND) then
              call carrysec(iy, imon, id, ih, im, dummy)
            end if
            write (lfnrck, '(5i6,f10.6,6f17.6)') iy, imon, id, ih, im, dummy, (rck(j), j=1, 6)
            rck = 0.d0
            rckref = mw
          else if (rckref .eq. 0.d0) then
            rckref = mw
          end if
          if (PM(ipar)%pname(1:8) .eq. 'RECCLK_G') then
            rck(1) = PM(ipar)%xini + PM(ipar)%xcor
          else if (PM(ipar)%pname(1:8) .eq. 'RECCLK_R') then
            rck(2) = PM(ipar)%xini + PM(ipar)%xcor
          else if (PM(ipar)%pname(1:8) .eq. 'RECCLK_E') then
            rck(3) = PM(ipar)%xini + PM(ipar)%xcor
          else if (PM(ipar)%pname(1:8) .eq. 'RECCLK_C') then
            rck(4) = PM(ipar)%xini + PM(ipar)%xcor
          else if (PM(ipar)%pname(1:8) .eq. 'RECCLK_3') then
            rck(5) = PM(ipar)%xini + PM(ipar)%xcor
          else if (PM(ipar)%pname(1:8) .eq. 'RECCLK_J') then
            rck(6) = PM(ipar)%xini + PM(ipar)%xcor
          end if
        else if (PM(ipar)%pname(1:4) .eq. 'HTGC') then
          if (LCF%htgmod(1:3) .ne. 'NON') then
            call mjd2date(0, PM(ipar)%ptime(1), iy, imon, id, ih, im, xini)
            if (dabs(xini - 60.d0) .lt. MAXWND) then
              call carrysec(iy, imon, id, ih, im, xini)
            end if
            write (lfnhtg, '(5i6,f10.6,$)') iy, imon, id, ih, im, xini
            call mjd2date(0, PM(ipar)%ptime(2), iy, imon, id, ih, im, xini)
            if (dabs(xini - 60.d0) .lt. MAXWND) then
              call carrysec(iy, imon, id, ih, im, xini)
            end if
            write (lfnhtg, '(5i6,f10.6,4f11.6)') iy, imon, id, ih, im, xini, &
              (PM(ipar + i)%xini, PM(ipar + i)%xcor, i=0, 1)
          end if
        else if (PM(ipar)%pname(1:5) .eq. 'STAPX') then
          cid = ' '
          if (PM(ipar)%iobs .lt. 5) cid = ' *'
          t0 = PM(ipar)%ptime(1)
          if (LCF%jd0 + LCF%sod0/864.d2 .gt. t0) t0 = LCF%jd0 + LCF%sod0/864.d2
          t1 = PM(ipar)%ptime(2)
          if (LCF%jd1 + LCF%sod1/864.d2 .lt. t1) t1 = LCF%jd1 + LCF%sod1/864.d2
          mw = (t0 + t1)/2.d0
          call xyzblh(PM(ipar:ipar + 2)%xini + PM(ipar:ipar + 2)%xcor, 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, SITE%geod)
          write (SITE%ikin, '(i5,f10.2,a2,3f14.4,2f17.11,f14.4,i7,6(1x,i2.2),f7.2)') int(mw), (mw - int(mw))*864.d2, cid, &
            (PM(ipar + i)%xini + PM(ipar + i)%xcor, i=0, 2), &
            SITE%geod(1)/PI*180.d0, SITE%geod(2)/PI*180.d0, SITE%geod(3), &
            PM(ipar)%iobs/PM(ipar)%iepo, PM(ipar)%iobs_G/PM(ipar)%iepo,   &
            PM(ipar)%iobs_R/PM(ipar)%iepo, PM(ipar)%iobs_E/PM(ipar)%iepo, &
            PM(ipar)%iobs_C/PM(ipar)%iepo, PM(ipar)%iobs_3/PM(ipar)%iepo, &
            PM(ipar)%iobs_J/PM(ipar)%iepo, pdop
        end if
      end if
      backspace lfnrem
    else if (cid .eq. 'pz') then
      backspace lfnrem
      read (lfnrem) ipar, zdd, zwd, dummy
      if ((dummy - PM(ipar)%ptime(1))*86400.d0 .gt. -MAXWND .and. &
          (PM(ipar)%ptime(2) - dummy)*86400.d0 .gt. -MAXWND .and. PM(ipar)%iobs .gt. 0) then
        call mjd2date(0, dummy, iy, imon, id, ih, im, mw)
        if (dabs(mw - 60.d0) .lt. maxwnd) then
          call carrysec(iy, imon, id, ih, im, mw)
        end if
        if (PM(ipar)%pname(1:3) .eq. 'ZTD') then
          if (im .ge. 60) then
            print *, dummy, iy, imon, id, ih, im, mw
            call exit(1)
          end if
          write (lfnztd, '(5i6,f10.6,3f11.6)') iy, imon, id, ih, im, mw, zdd, zwd, PM(ipar)%xcor
        end if
      end if
      backspace lfnrem
    else if (cid .eq. 'am') then
      backspace lfnrem
      read (lfnrem) ipar, ntot, ipt(0), (ipt(k), k=1, ntot - 1), (cof(k), k=1, ntot)
      NM%ns = NM%ns + 1
      do i = 1, ntot - 1
        if (ipt(i) .ne. ipar) then
          cof(ntot) = cof(ntot) - cof(i)*PM(ipt(i))%xcor
        end if
      end do
      PM(ipar)%xcor = cof(ntot)/cof(ipt(0))
      backspace lfnrem
    else if (cid .eq. 'dp') then
      backspace lfnrem
      read (lfnrem) pdop
      backspace lfnrem
    end if
    backspace lfncid
  end do
!
!! wirte residual file
  if (jdr .ne. 0) then
    call mjd2date(jdr, sodr, iy, imon, id, ih, im, dummy)
    write (lfnres, '(a3,i5,4i3,f11.7,i7,f10.2)') 'TIM', iy, imon, id, ih, im, dummy, jdr, sodr
    LCF%jd_beg = jdr
    LCF%sod_beg = sodr
  end if
  write (lfnres, '(60x,a)') 'END OF HEADER'
  write (lfnres, '(a15,45x,a)') LCF%amb_at_dbd, 'AMB AT DAT-BOUNDARY'
  do i = MAXSYS, 1, -1
    write (lfnres, '(a3,7x,2(a3,2x),40x,a)') GNSS_NAME_LEN3(i:i), &
      FREQ_NAME_SYS(idxfrq(i, 1), i), FREQ_NAME_SYS(idxfrq(i, 2), i), 'SYS / FREQUENCY BAND'
  end do
  call mjd2date(LCF%jd_beg, LCF%sod_beg, iy, imon, id, ih, im, dummy)
  write (lfnres, '(i5,4i3,f11.7,f12.2,20x,a)') iy, imon, id, ih, im, dummy, &
    timdif(LCF%jd_end, LCF%sod_end, LCF%jd_beg, LCF%sod_beg), 'RES TIME BEG/LEN'
  write (lfnres, '(f10.2,10x,a12,28x,a)') LCF%dintv, 'LCPC (meter)', 'INT / OBS TYPE'
  write (lfnres, '(f10.3,50x,a)') NM%sig0, 'WEIGHTED SIGMA (METER)'
  write (lfnres, '(2i10,40x,a)') NM%nuk, NM%nobs, '# OF UNKNOWN / OBS'
  ! sort prn
  do i = 1, nprn - 1
    do j = i + 1, nprn
      ig = pointer_string(LCF%nprn, LCF%prn, prn(i))
      ik = pointer_string(LCF%nprn, LCF%prn, prn(j))
      if (ig .gt. ik) then
        temprn = prn(i)
        prn(i) = prn(j)
        prn(j) = temprn
      end if
    end do
  end do
  k = nprn/15 + 1
  do j = k, 1, -1
    isat = (j - 1)*15 + 1
    if (isat .gt. nprn) cycle
    if (j .eq. k) then
      line = ' '
      write (line, '(15(a3,1x))') (prn(i), i=isat, nprn)
      write (lfnres, '(a)') line//'SATELLITE LIST'
    else
      write (lfnres, '(15(a3,1x),a)') (prn(i), i=isat, isat + 14), 'SATELLITE LIST'
    end if
  end do
  write (lfnres, '(i10,50x,a)') nprn, '# OF SAT'
  write (lfnres, '(f8.2,52x,a)') SITE%cutoff/PI*180.d0, 'OBS MASK ANGLE (deg)'
  if (SITE%skd(1:1) .eq. 'S') then
    write (lfnres, '(a6,5x,3f10.6,19x,a)') 'Static', (SITE%dx0(i), i=1, 3), 'POS MODE/PRIORI (meter)'
  else if (SITE%skd(1:1) .eq. 'P') then
    write (lfnres, '(a10,1x,3f10.6,19x,a)') 'Piece-wise', (SITE%dx0(i), i=1, 3), 'POS MODE/PRIORI (meter)'
  else if (SITE%skd(1:1) .eq. 'K') then
    write (lfnres, '(a9,2x,3f10.6,19x,a)') 'Kinematic', (SITE%dx0(i), i=1, 3), 'POS MODE/PRIORI (meter)'
  else if (SITE%skd(1:1) .eq. 'F') then
    write (lfnres, '(a11,3f10.6,19x,a)') 'Quasi-fixed', (SITE%dx0(i), i=1, 3), 'POS MODE/PRIORI (meter)'
  else if (SITE%skd(1:1) .eq. 'L') then
    write (lfnres,'(a3,8x,3f10.6,19x,a)') 'LEO', (SITE%dx0(i), i=1, 3), 'POS MODE/PRIORI (meter)'
  else
    write (lfnres, '(a7,4x,3f10.6,19x,a)') 'Unknown', (SITE%dx0(i), i=1, 3), 'POS MODE/PRIORI (meter)'
  end if
  write (lfnres, '(a4,56x,a)') SITE%name, 'STATION'
  write (lfnres, '(a9,51x,a)') 'Residuals', 'COMMENT'
  close (lfnres)
  iunit = get_valid_unit(10)
  call reverse_file(LCF%flnres, iunit)
!
!! write receiver clock file
  if (rckref .ne. 0.d0) then
    call mjd2date(0, rckref, iy, imon, id, ih, im, dummy)
    if (dabs(dummy - 60.d0) .lt. MAXWND) then
      call carrysec(iy, imon, id, ih, im, dummy)
    end if
    write (lfnrck, '(5i6,f10.6,6f17.6)') iy, imon, id, ih, im, dummy, (rck(j), j=1, 6)
  end if
  call lsq_wrt_header(lfnrck, LCF, SITE, OB, 'rck', .false., .true., .true., .true.)
  close (lfnrck)
  iunit = get_valid_unit(10)
  call reverse_file(LCF%flnrck, iunit)
!
!! write atm zenith delay file
  if (LCF%ztdmod(1:3) .ne. 'NON') then
    call lsq_wrt_header(lfnztd, LCF, SITE, OB, 'ztd', .false., .true., .true., .true.)
    close (lfnztd)
    iunit = get_valid_unit(10)
    call reverse_file(LCF%flnztd, iunit)
  end if
!
!! write horizontal troposphere gradient
  if (LCF%htgmod(1:3) .ne. 'NON') then
    call lsq_wrt_header(lfnhtg, LCF, SITE, OB, 'htg', .false., .true., .true., .true.)
    close (lfnhtg)
    iunit = get_valid_unit(10)
    call reverse_file(LCF%flnhtg, iunit)
  end if
!
!! write kinematic files
  if (SITE%skd(1:1) .eq. 'P' .or. SITE%skd(1:1) .eq. 'K' .or. SITE%skd(1:1) .eq. 'L') then
    call lsq_wrt_header(SITE%ikin, LCF, SITE, OB, 'kin', .false., .true., .true., .true.)
    close (SITE%ikin)
    iunit = get_valid_unit(10)
    call reverse_file(SITE%kinfil, iunit)
  end if

  return
100 write (*, '(a)') '***ERROR(lsq_rcv_prmt): back space or read file'
  call exit(1)

contains

  subroutine carrysec(iy, imon, iday, ih, im, fsec)
    implicit none
    integer*4 iy, imon, iday, ih, im
    real*8 fsec
    integer*4 imjd, idoy
    integer*4 modified_julday
    im = im + 1
    fsec = 0.d0
    if (im .ge. 60) then
      ih = ih + 1
      im = im - 60
    end if
    if (ih .ge. 24) then
      imjd = modified_julday(iday, imon, iy)
      call mjd2doy(imjd + 1, iy, idoy)
      call yeardoy2monthday(iy, idoy, imon, id)
      ih = ih - 24
    end if
  end subroutine

end subroutine
