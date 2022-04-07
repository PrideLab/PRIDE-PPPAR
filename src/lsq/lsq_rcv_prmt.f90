!
!! lsq_rcv_prmt.f90
!!
!!    Copyright (C) 2022 by Wuhan University
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
!! Contributor: Maorong Ge, Jianghui Geng, Songfeng Yang, Jihang Lin
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
  integer*4 i, j, k, ig, ik, ip(2), ipar, isat, ntot, nprn, ipt(0:MAXPAR), jd, jdw, jdx, ierr
  integer*4 lfnrck, lfnztd, lfnhtg, jdr, iy, imon, id, ih, im, idoy
  real*8 cof(MAXPAR), x(3), xc(6), sodr, phase, range, wphs, wrng, mw, elev, azim, dmap, wmap, pdop
  real*8 zdd, zwd, xini, xcor, dummy, rckref, sodw, sodx
!
  character*256 cflip
  character*60 line
  character*3 temprn,prn(MAXSAT),typuse(4)
  real*8 rck(1:6)
!
!! function called
  integer*4 get_valid_unit,pointer_string
  real*8 timdif

  cflip='newtac() { out=`uname`; if [ "$out" = "Darwin" ]; then tail -r $1 > tmp_scratch; &
                    else tac $1 > tmp_scratch; fi; mv -f tmp_scratch $1; }; newtac '
  lfnres = get_valid_unit(10)
  open (lfnres, file=LCF%flnres)
  lfnrck = get_valid_unit(10)
  open (lfnrck, file=LCF%flnrck)
  lfnztd = get_valid_unit(10)
  open (lfnztd, file=LCF%flnztd)
  if (LCF%htgmod(1:3) .ne. 'NON') then
    lfnhtg = get_valid_unit(10)
    open (lfnhtg, file=LCF%flnhtg)
  endif
  if (index(SITE%skd, 'K') .ne. 0) then ! kinematic
    if (SITE%ikin .ne. 0) close(SITE%ikin)
    SITE%ikin = get_valid_unit(10)
    open (SITE%ikin, file=SITE%kinfil)
  endif
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
                    phase, range, wphs, wrng, iflg, elev, azim, dmap, wmap,typuse(1:4)
      dummy = 0.d0
      do i = 1, ntot
        dummy = dummy + PM(ipt(i))%xcor*cof(i)
      enddo
      if (ipt(0) .eq. 0) then
        range = range - dummy
      else
        phase = phase - dummy
        range = range - dummy + PM(ipt(ipt(0)))%xcor*cof(ipt(0))
        PM(ipt(ipt(0)))%xrms = PM(ipt(ipt(0)))%xrms + phase**2
      endif
      if (jdr .eq. 0) then
        jdr = jd
        sodr = mw
        LCF%jd_end = jd
        LCF%sod_end = mw
      else if((jdr - jd)*86400.d0 + sodr - mw .gt. MAXWND) then
        call mjd2date(jdr, sodr, iy, imon, id, ih, im, dummy)
        write (lfnres, '(a3,i5,4i3,f11.7,i7,f10.2)') 'TIM', iy, imon, id, ih, im, dummy, jdr, sodr
        jdr = jd
        sodr = mw
      endif
      elev = elev/PI*180.d0
      azim = azim/PI*180.d0
      write (lfnres, '(a3,2f10.3,2d16.8,i3,f8.3,f9.3,4(1x,a3))') &
             LCF%prn(isat), phase, range, wphs, wrng, iflg, elev, azim, typuse(1:4)
      if(pointer_string(nprn,prn,LCF%prn(isat)).eq.0) then
        nprn=nprn+1
        prn(nprn)=LCF%prn(isat)
      endif
      backspace lfnobs
    else if (cid .eq. 'cn') then
      backspace lfnobs
      read(lfnobs) jdw,sodw,jdx,sodx,i,j,ip(1:2),phase,wphs
      if(jdr.ne.0.and.(jdr - jdx)*86400.d0 + sodr - sodx .gt. MAXWND) then
        call mjd2date(jdr, sodr, iy, imon, id, ih, im, dummy)
        write (lfnres, '(a3,i5,4i3,f11.7,i7,f10.2)') 'TIM', iy, imon, id, ih, im, dummy, jdr, sodr
        jdr = jdx
        sodr = sodx
      endif
      phase=phase-PM(ip(1))%xcor+PM(ip(2))%xcor
      write(lfnres, '(a3,f20.15,d16.8,2x,2(i7,f10.2),2(1x,a3))') 'CST',phase,wphs,jdw,sodw,jdx,sodx,LCF%prn(i),LCF%prn(j)
      backspace lfnobs
    else if (cid .eq. 'de') then
      backspace lfnrem
      read (lfnrem) jd, mw, isat, iflg
      if (jdr .eq. 0) then
        jdr = jd
        sodr = mw
        LCF%jd_end = jd
        LCF%sod_end = mw
      else if((jdr - jd)*86400.d0 + sodr - mw .gt. MAXWND) then
        call mjd2date(jdr, sodr, iy, imon, id, ih, im, dummy)
        write (lfnres, '(a3,i5,4i3,f11.7,i7,f10.2)') 'TIM', iy, imon, id, ih, im, dummy, jdr, sodr
        jdr = jd
        sodr = mw
      endif
      write (lfnres, '(a3,2f10.3,2d16.8,i3,f8.3,f9.3)') LCF%prn(isat), 0.d0, 0.d0, 0.d0, 0.d0, iflg, 0.d0, 0.d0
      if(pointer_string(nprn,prn,LCF%prn(isat)).eq.0) then
        nprn=nprn+1
        prn(nprn)=LCF%prn(isat)
      endif
      backspace lfnrem
    else if (cid .eq. 'pc') then
      backspace lfnrem
      read (lfnrem) ipar, ntot, ipt(0), (ipt(k), k=1, ntot - 1), (cof(k), k=1, ntot), &
                    PM(ipar)%iobs, PM(ipar)%iobs_G, PM(ipar)%iobs_R, PM(ipar)%iobs_E, &
                    PM(ipar)%iobs_C, PM(ipar)%iobs_3, PM(ipar)%iobs_J, &
                    PM(ipar)%xini, PM(ipar)%ptime(1:2)
      dummy = cof(ipt(0))
      cof(ipt(0)) = -PM(ipar)%map*PM(ipar)%rw**2
      if (PM(ipar)%ipt .lt. 0) then
        PM(ipar)%ipt = 0
        cof(ipt(0)) = 0.d0
        NM%np = NM%np + 1
      endif
      do i = 1, ntot - 1
        cof(ntot) = cof(ntot) - cof(i)*PM(ipt(i))%xcor
      enddo
      PM(ipar)%xcor = cof(ntot)/dummy
      if (PM(ipar)%iobs .gt. 0) then
        mw = PM(ipar)%ptime(1)
        if (PM(ipar)%pname(1:6) .eq. 'RECCLK') then
          if (rckref.ne.0.d0.and.(rckref-mw)*86400.d0.gt.MAXWND) then
            call mjd2date(0, rckref, iy, imon, id, ih, im, dummy)
            if (dabs(dummy-60.d0) .lt. MAXWND) then
              call carrysec(iy, imon, id, ih, im, dummy)
            endif
            write (lfnrck, '(5i6,f10.6,6f17.6)') iy, imon, id, ih, im, dummy, (rck(j), j=1, 6)
            rck = 0.d0
            rckref = mw
          else if (rckref.eq.0.d0) then
            rckref = mw
          endif
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
          endif
        else if (PM(ipar)%pname(1:4) .eq. 'HTGC') then
          if (LCF%htgmod(1:3) .ne. 'NON') then
            call mjd2date(0, PM(ipar)%ptime(1), iy, imon, id, ih, im, xini)
            write (lfnhtg, '(5i6,f10.6,$)') iy, imon, id, ih, im, xini
            call mjd2date(0, PM(ipar)%ptime(2), iy, imon, id, ih, im, xini)
            write (lfnhtg, '(5i6,f10.6,4f11.6)') iy, imon, id, ih, im, xini, &
                   (PM(ipar + i)%xini, PM(ipar + i)%xcor, i=0, 1)
          endif
        else if (PM(ipar)%pname(1:5) .eq. 'STAPX') then
          cid = ' '
          if (PM(ipar)%iobs .lt. 5) cid = ' *'
          call xyzblh(PM(ipar : ipar+2)%xini + PM(ipar : ipar+2)%xcor, 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, SITE%geod)
          write (SITE%ikin, '(i5,f10.2,a2,3f13.3,2f16.10,f13.3,i7,6(1x,i2.2),f7.2)') int(mw), (mw - int(mw))*864.d2, cid, &
                 (PM(ipar + i)%xini + PM(ipar + i)%xcor, i=0, 2), &
                 SITE%geod(1)/PI*180.d0, SITE%geod(2)/PI*180.d0, SITE%geod(3), PM(ipar)%iobs, PM(ipar)%iobs_G, &
                 PM(ipar)%iobs_R, PM(ipar)%iobs_E, PM(ipar)%iobs_C, PM(ipar)%iobs_3, PM(ipar)%iobs_J, pdop
        endif
      endif
      backspace lfnrem
    else if (cid .eq. 'pz') then
      backspace lfnrem
      read (lfnrem) ipar, zdd, zwd, dummy
      if ((dummy - PM(ipar)%ptime(1))*86400.d0 .gt. -MAXWND .and. &
          (PM(ipar)%ptime(2) - dummy)*86400.d0 .gt. -MAXWND .and. PM(ipar)%iobs .gt. 0) then
        call mjd2date(0, dummy, iy, imon, id, ih, im, mw)
        if (dabs(mw-60.d0) .lt. maxwnd) then
          call carrysec(iy, imon, id, ih, im, mw)
        endif
        if (PM(ipar)%pname(1:3) .eq. 'ZTD') then
          if (im .ge. 60) then
            print *, dummy, iy, imon, id, ih, im, mw
            call exit(1)
          endif
          write (lfnztd, '(5i6,f10.6,3f11.6)') iy, imon, id, ih, im, mw, zdd, zwd, PM(ipar)%xcor
        endif
      endif
      backspace lfnrem
    else if (cid .eq. 'am') then
      backspace lfnrem
      read (lfnrem) ipar, ntot, ipt(0), (ipt(k), k=1, ntot - 1), (cof(k), k=1, ntot)
      NM%ns = NM%ns + 1
      do i = 1, ntot - 1
        if (ipt(i) .ne. ipar) then
          cof(ntot) = cof(ntot) - cof(i)*PM(ipt(i))%xcor
        endif
      enddo
      PM(ipar)%xcor = cof(ntot)/cof(ipt(0))
      backspace lfnrem
    else if(cid.eq.'dp') then
      backspace lfnrem
      read(lfnrem) pdop
      backspace lfnrem
    endif
    backspace lfncid
  enddo
!
!! wirte residual file
  if(jdr .ne. 0) then
    call mjd2date(jdr, sodr, iy, imon, id, ih, im, dummy)
    write (lfnres, '(a3,i5,4i3,f11.7,i7,f10.2)') 'TIM', iy, imon, id, ih, im, dummy, jdr, sodr
    LCF%jd_beg=jdr
    LCF%sod_beg=sodr
  endif
  write (lfnres, '(60x,a)') 'END OF HEADER'
  call mjd2date(LCF%jd_beg, LCF%sod_beg, iy, imon, id, ih, im, dummy)
  write (lfnres, '(i5,4i3,f11.7,f12.2,20x,a)') iy, imon, id, ih, im, dummy, &
         timdif(LCF%jd_end,LCF%sod_end,LCF%jd_beg,LCF%sod_beg), 'RES TIME BEG/LEN'
  write (lfnres, '(f10.2,10x,a11,29x,a)') LCF%dintv, 'LCPC(meter)', 'INT / OBS TYPE'
  write (lfnres, '(f10.3,50x,a)') NM%sig0, 'WEIGHTED SIGMA (METER)'
  write (lfnres, '(2i10,40x,a)') NM%nuk, NM%nobs, '# OF UNKNOWN / OBS'
  ! sort prn
  do i=1,nprn-1
    do j=i+1,nprn
      ig=pointer_string(LCF%nprn, LCF%prn, prn(i))
      ik=pointer_string(LCF%nprn, LCF%prn, prn(j))
      if(ig .gt. ik) then
        temprn=prn(i)
        prn(i)=prn(j)
        prn(j)=temprn
      endif
    enddo
  enddo
  k = nprn/15 + 1
  do j = k, 1, -1
    isat = (j-1)*15 + 1
    if (isat.gt.nprn) cycle
    if (j .eq. k) then
      line = ' '
      write (line, '(15(a3,1x))') (prn(i), i=isat, nprn)
      write (lfnres, '(a)') line//'SATELLITE LIST'
    else
      write (lfnres, '(15(a3,1x),a)') (prn(i), i=isat, isat + 14), 'SATELLITE LIST'
    endif
  enddo
  write (lfnres, '(i10,50x,a)') nprn, '# OF SAT'
  write (lfnres, '(a4,56x,a)') SITE%name, 'STATION'
  write (lfnres, '(a9,51x,a)') 'Residuals', 'COMMENT'
  call system(cflip//LCF%flnres)
  close (lfnres)
!
!! write receiver clock file
  if(rckref.ne.0.d0) then
    call mjd2date(0, rckref, iy, imon, id, ih, im, dummy)
    if (dabs(dummy-60.d0) .lt. MAXWND) then
      call carrysec(iy, imon, id, ih, im, dummy)
    endif
    write (lfnrck, '(5i6,f10.6,6f17.6)') iy, imon, id, ih, im, dummy, (rck(j), j=1, 6)
  endif
  call lsq_wrt_header(lfnrck, LCF, SITE, OB, 'rck', .false., .true., .true., .true.)
  call system(cflip//LCF%flnrck)
  close (lfnrck)
!
!! write atm zenith delay file
  call lsq_wrt_header(lfnztd, LCF, SITE, OB, 'ztd', .false., .true., .true., .true.)
  call system(cflip//LCF%flnztd)
  close (lfnztd)
!
!! write horizontal troposphere gradient
  if (LCF%htgmod(1:3) .ne. 'NON') then
    call lsq_wrt_header(lfnhtg, LCF, SITE, OB, 'htg', .false., .true., .true., .true.)
    call system(cflip//LCF%flnhtg)
    close (lfnhtg)
  endif
!
!! write kinematic files
  if (index(SITE%skd, 'K') .ne. 0) then
    call lsq_wrt_header(SITE%ikin, LCF, SITE, OB, 'kin', .false., .true., .true., .true.)
    call system(cflip//SITE%kinfil)
  endif
  if (index(SITE%skd, 'K') .ne. 0) close (SITE%ikin)

  return
100 write (*, '(a)') '***ERROR(lsq_rcv_prmt): back space or read file'
  call exit(1)

contains

  subroutine carrysec(iy, imon, iday, ih, im, fsec)
    implicit none
    integer*4 iy, imon, iday, ih, im
    real*8    fsec
    integer*4 imjd, idoy
    integer*4 modified_julday
    im = im + 1
    fsec = 0.d0
    if (im .ge. 60) then
      ih = ih + 1
      im = im - 60
    endif
    if (ih .ge. 24) then
      imjd = modified_julday(iday, imon, iy)
      call mjd2doy(imjd+1, iy, idoy)
      call yeardoy2monthday(iy, idoy, imon, id)
      ih = ih - 24
    endif
  end subroutine

end subroutine
