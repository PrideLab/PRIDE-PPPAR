!
!! read_rinex_file.f90
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
!! flnrnx ------------ rinex O_file
! tstart,sstart ----- start time: tstart(1):year, tstart(2):month, tstart(3):day, tstart(4):hour, tstart(5):minute, sstart:second.
! session_length ---- sesstion time length
! interval ---------- Epoch interval of processing. (rinex o_file Epoch interval: OH.intv)
! check_pc
! pclimit -----------
! cutoff_elevation -- cut off elevation
! use_brdeph -------- true or false
! neph -------------- the number of eph
! ephem -------------
! stanam ------------
! x,y,z -------------
! t_first_in_rinex --
! t_last_in_rinex ---
! v -----------------
subroutine read_rinex_file(flnrnx, tstart, sstart, session_length, interval, &
                           check_pc, pclimit, &
                           cutoff_elevation, use_brdeph, neph, ephem, lm_edit, ltighter, &
                           stanam, x, y, z, t_first_in_rinex, t_last_in_rinex, v, &
                           nepo, nsat, jd0, nobs, flagall, ti, ts, obs, itypuse, bias,dwnd,GNSS_SYS)
  implicit none
  include '../header/const.h'
  include '../header/absbia.h'
  include '../header/brdeph.h'
  include '../header/rnxobs.h'

! common
  integer*4     idxfrq(MAXSYS, 2)
  common        idxfrq
! parameter
  character*(*) flnrnx
  integer*4     tstart(1:*)
  real*8        sstart, session_length, interval
  logical*1     check_pc
  real*8        pclimit
  real*8        cutoff_elevation
  logical*1     use_brdeph, lm_edit, ltighter
  integer*4     neph
  type(brdeph)  ephem(1:*)
  character*(*) stanam
  integer*4     nepo, nsat, jd0
  real*4        v(nepo, MAXSAT)
  real*8        x(nepo), y(nepo), z(nepo)
  integer*4     nobs(MAXSAT), nbias_used(MAXSAT)
  integer*4     flagall(nepo, MAXSAT)
  real*8        ti(nepo), ts(nepo)
  real*8        obs(nepo, MAXSAT, 7)
  type(absbia)  bias(MAXSAT, MAXTYP)
  integer*4     itypuse(MAXSAT, 4)
! local
  type(rnxhdr)  HD
  type(rnxobr)  OB
  integer*4     nprn, nused
  integer*4     ichn, ieph, iepo, isat, ilast(MAXSAT)
  integer*4     ierr
  integer*4     iunit, iunit_next
  integer*4     prn_int
  integer*4     i0, ij, j, k
  integer*1     lli_thre
  real*8        lglimit, nwlimit, lglimit_t
  character*3   jj
  character*3   prn(MAXSAT), prn0(MAXSAT)
  real*8        tobs, nwdif, rgdif, lgdif, elev, range, dtsat, dt, dwnd
  real*8        t_first_in_rinex, t_last_in_rinex, vv
!! coefficient
  real*8        f1(MAXSYS), f2(MAXSYS), rg(MAXSYS), c1(MAXSYS), c2(MAXSYS)
  real*8        lam1(MAXSYS), lam2(MAXSYS), lamw(MAXSYS)
!! GLONASS
  real*8        FREQ1_R(-50:50), FREQ2_R(-50:50)
!! range check
  integer*4     flg(MAXSAT), it(MAXSAT)
  real*8        pmb(MAXSAT), wpmb(MAXSAT), rms, sig
!! receiver clock jump check
  real*8        jumpsum(2), jumpval(2), rangeval(2, MAXSAT), deltap, deltal
  real*8        tmp(2), ratio, sec, obsval(0:2, 4, MAXSAT), chkval(2, MAXSAT)
  integer*4     lfnjmp, nvalid, njump, cnt
  integer*4     jumpflag    !! 1:range, 2:phase
! funtion called
  logical*1     istrue
  integer*4     modified_julday, set_flag, get_valid_unit, pointer_string
  real*8        timdif
  character*8   GNSS_SYS 

! initialize
  call prn_matbld(prn0)
  call frequency_glonass(FREQ1_R, FREQ2_R)

  do i0 = 1, MAXSYS
  !! frequency
    f1(i0) = FREQ_SYS(idxfrq(i0, 1), i0)
    if (f1(i0) .eq. 0.d0) goto 100
    f2(i0) = FREQ_SYS(idxfrq(i0, 2), i0)
    if (f2(i0) .eq. 0.d0) goto 100
  !! ratio and coeff
    rg(i0) = f1(i0)/f2(i0)
    c1(i0) = 1.d0/(1.d0 - 1.d0/rg(i0)**2)
    c2(i0) = -1.d0/(rg(i0)**2 - 1.d0)
  !! wave length
    lam1(i0) = VLIGHT/f1(i0)
    lam2(i0) = VLIGHT/f2(i0)
    lamw(i0) = VLIGHT/(f1(i0) - f2(i0))
  end do

  lglimit = 50.d0
  nwlimit = 50.d0
  if (ltighter) nwlimit = 10.d0
  lfnjmp = 0
  jumpsum = 0.d0
  obsval = 0.d0
  lli_thre = 3
  itypuse = 0

  !dwnd = min(interval/10.d0, 0.01d0)

! open rinex observation file
  iunit = 10
  open (iunit, file=flnrnx, form='FORMATTED', status='OLD', iostat=ierr)
  if (ierr .ne. 0) then
    write (*, '(a)') '***ERROR(read_rinex_file): open file error, '//flnrnx
    call exit(1)
  end if

! read head of rinex file
  call rdrnxoh(iunit, HD, ierr)
  if (ierr .ne. 0) then
    call exit(1)
  end if

! start time of first session
  if (tstart(2) .eq. 0) then
    do j = 1, 5
      tstart(j) = HD%t0(j)
    end do
    sstart = HD%t0s
  end if
  call yr2year(tstart(1))

!! jd0:  julday
!! tobs: fractional julday
  jd0 = modified_julday(tstart(3), tstart(2), tstart(1))
  tobs = tstart(4)*3600.d0 + tstart(5)*60.d0 + sstart

! initialize
  nsat = 0
  do isat = 1, MAXSAT
    ilast(isat) = 0
  end do

! ti --- time of each processing epoch
! nepo - the number of epoch
  ti(1) = tobs
  do iepo = 1, nepo
    if (iepo .ge. 2) ti(iepo) = ti(iepo - 1) + interval
    ts(iepo) = ti(iepo)
    do isat = 1, MAXSAT
      flagall(iepo, isat) = set_flag(0, 'nodata')
      v(iepo, isat) = 0.d0
      do j = 1, 6
        obs(iepo, isat, j) = 0.d0
      end do
    end do
  end do
!
! loop over session and epoch
  iepo = 1
  do while (iepo .le. nepo)
    if (ierr .ne. 0 .and. ierr .ne. 2) cycle
    nprn = 0
    prn = ''
! ***************************************************************************************** !
! purpose  : read one epoch data from a RINEX o-file
! parameter: jd0,tobs ----------- start time
!            dwnd --------------- window for time matching. If the obsersing time from the
!                                 rinex-file is close to the requested time within the
!                                 window, we take the data. Be careful with this parameter
!                                 when you are working sampling rate larger than 1Hz.
!            nprn,prn ----------- satellite prn number
!            HD ----------------- rinex header structure
!            OB ----------------- o_file body structure
! **************************************************************************************** !
50 continue
    call read_glschn(jd0, tobs, OB%glschn)
    if (HD%ver .ge. 300) then
      call rdrnxoi3(iunit, jd0, tobs, dwnd, nprn, prn, HD, OB, bias, nbias_used, ierr)
    else
      call rdrnxoi2(iunit, jd0, tobs, dwnd, nprn, prn, HD, OB, bias, nbias_used, ierr)
    end if
    ! judge satellite system
    do ichn = 1, OB%nprn
      read (OB%prn(ichn) (2:3), '(i2)') j
      if (index(GNSS_SYS, OB%prn(ichn)(1:1)) > 0) then
        continue  
      else if (OB%prn(ichn)(1:1)  .eq. 'C') then
        if  (index(GNSS_SYS, '2') > 0 .and. index(GNSS_SYS, '3') > 0 ) then
          continue
        else if  (index(GNSS_SYS, '2') > 0) then
          if (j .gt. 17) then
            OB%prn(ichn)=" "
            continue
          end if 
        else if  (index(GNSS_SYS, '3') > 0) then
          if (j .le. 17) then
            OB%prn(ichn)=" "
            continue
          end if 
        else
          OB%prn(ichn)=" "  
          continue
        end if
      else
        OB%prn(ichn)=" "
      endif
    end do
    if (ierr .eq. 2) then
      call next_rinex(iunit, iunit_next, &
        int(jd0 + (tstart(4)*3600.d0 + tstart(5)*60.d0 + sstart + session_length)/864.d2))
      if (iunit_next .eq. 0) exit
      close (iunit)
      iunit = iunit_next
      goto 50
    end if
    if (ierr .ne. 0) continue
    obsval(0, :, :) = obsval(1, :, :)
    obsval(1, :, :) = obsval(2, :, :)
    obsval(2, :, :) = 0.d0
    ti(iepo) = tobs
    if (OB%nprn .ne. 0) then
      ti(iepo) = ti(iepo) - timdif(jd0, tobs, OB%jd, OB%tsec)
    end if
    do ichn = 1, OB%nprn
      jj = OB%prn(ichn)
      j = pointer_string(MAXSAT, prn0, jj)
      if (j .eq. 0) cycle
! nsat : largest PRN for a satellite
      if (nsat .lt. j) nsat = j
! itypuse: index of obs type used at the day boundary
      if ((abs(OB%tsec) .le. MAXWND .or. abs(OB%tsec - 864.d2) .le. MAXWND)) then
        itypuse(j, 1:4) = OB%itypuse(ichn, 1:4)
      end if
      if (use_brdeph) then
        ieph = 0
! range --------- the distance between satellite and reciver
! dtsat --------- clock correction, units: meter
! vv    --------- true anomaly
        call elevation(neph, ephem, jj, jd0, ti(iepo), x(iepo), y(iepo), z(iepo), elev, range, dtsat, vv, ieph, .false.)
        if (range .lt. 0.d0) then
          flagall(iepo, j) = set_flag(0, 'pcbad')
          cycle
        end if
      end if
      i0 = index(GNSS_PRIO, jj(1:1))
      if ('R' .eq. jj(1:1)) then
        read (OB%prn(ichn), '(1x,i2)') prn_int
        lam1(i0) = VLIGHT/FREQ1_R(OB%glschn(prn_int))
        lam2(i0) = VLIGHT/FREQ2_R(OB%glschn(prn_int))
        lamw(i0) = VLIGHT/(FREQ1_R(OB%glschn(prn_int)) - FREQ2_R(OB%glschn(prn_int)))
      end if
      if (all(dabs(OB%obs(ichn, 1:4)) .gt. 1.d-3)) then
        flagall(iepo, j) = 0
        rangeval(1, j) = rangeval(2, j)
        rangeval(2, j) = range
        obsval(2, 2, j) = OB%obs(ichn, 1)*lam1(i0)  !! L1
        obsval(2, 4, j) = OB%obs(ichn, 2)*lam2(i0)  !! L2
        obsval(2, 1, j) = OB%obs(ichn, 3)           !! P1
        obsval(2, 3, j) = OB%obs(ichn, 4)           !! P2
        !! geometry-free : ionosphere observations
        obs(iepo, j, 1) = (lam1(i0)*OB%obs(ichn, 1) - lam2(i0)*OB%obs(ichn, 2))/(lam2(i0) - lam1(i0))
        ! **************************************************************** !
        !                       check geometry-free(lg)                    !
        ! **************************************************************** !
        if (ilast(j) .ne. 0) then
          if (.not. ltighter) then
            lgdif = dabs(obs(iepo, j, 1) - obs(ilast(j), j, 1))/(ti(iepo) - ti(ilast(j)))
            if (lgdif .gt. lglimit) then
              flagall(iepo, j) = set_flag(flagall(iepo, j), 'lgjump')
            end if
          else
            lgdif = dabs(obs(iepo, j, 1) - obs(ilast(j), j, 1))*(lam2(i0) - lam1(i0))
            dt = ti(iepo) - ti(ilast(j))
            if (dt .gt. 0.d0 .and. dt .le. 1.d0) then
              lglimit_t = 5.d-2
            else if (dt .gt. 1.d0 .and. dt .le. 2.d1) then
              lglimit_t = 1.d-1/2.d1*dt+5.d-2
            else if (dt .gt. 2.d1 .and. dt .le. 6.d1) then
              lglimit_t = 1.5d-1
            else if (dt .gt. 6.d1 .and. dt .le. 1.d2) then
              lglimit_t = 2.5d-1
            else
              lglimit_t = 3.5d-1
            end if
            if (lgdif .gt. lglimit_t) then
              flagall(iepo, j) = set_flag(flagall(iepo, j), 'lgjump')
            end if
          end if
        end if
        ! **************************************************************** ! 
        !                check Lost of lock indicator (LLI)                !
        ! **************************************************************** !
        if (lm_edit .and. (iand(OB%lli(ichn, 1), lli_thre) .ne. 0 .or. iand(OB%lli(ichn, 2), lli_thre) .ne. 0)) then
          flagall(iepo, j) = set_flag(flagall(iepo, j), 'lli')
        end if
      else
        flagall(iepo, j) = set_flag(0, 'no4')
      end if
    end do

    ! **************************************************************** !
    !               check & recover receiver clock jump                !
    ! **************************************************************** !
    nvalid = 0 ! number of satellites with flag 'ok'
    njump = 0
    cnt = 0
    jumpval = 0.d0
    deltap = 0.d0
    deltal = 0.d0
    do ichn = 1, OB%nprn
      jj = OB%prn(ichn)
      j = pointer_string(MAXSAT, prn0, jj)
      if (j .eq. 0) cycle
      if (.not. istrue(flagall(iepo, j), 'ok')) cycle
      if (dabs(obsval(0, 1, j)) .lt. 1.d-3 .or. dabs(obsval(1, 1, j)) .lt. 1.d-3 .or. dabs(obsval(2, 1, j)) .lt. 1.d-3) cycle
      nvalid = nvalid + 1
      chkval(1, j) = obsval(2, 1, j) - obsval(1, 1, j) - (obsval(2, 2, j) - obsval(1, 2, j)) ! f1
      chkval(2, j) = obsval(2, 3, j) - obsval(1, 3, j) - (obsval(2, 4, j) - obsval(1, 4, j)) ! f2
      if (dabs(chkval(1, j)) .gt. 20.d0 .and. dabs(chkval(2, j)) .gt. 20.d0) then            ! 0.1 us jump: 30.d0
        njump = njump + 1
        if (rangeval(1, j) .gt. 1.d0 .and. rangeval(2, j) .gt. 1.d0) then
          cnt = cnt + 1
          rgdif = rangeval(2, j) - rangeval(1, j)
          tmp(1) = (obsval(2, 1, j) - obsval(1, 1, j) + obsval(2, 3, j) - obsval(1, 3, j))/2.d0 - rgdif
          tmp(2) = (obsval(2, 2, j) - obsval(1, 2, j) + obsval(2, 4, j) - obsval(1, 4, j))/2.d0 - rgdif
          deltap = deltap + tmp(1)
          deltal = deltal + tmp(2)
        end if
        jumpval(1) = jumpval(1) + chkval(1, j)
        jumpval(2) = jumpval(2) + chkval(2, j)
      end if
    end do

    if (nvalid .ne. 0 .and. nvalid .eq. njump) then
      jumpval = jumpval/(nvalid)
      deltap = deltap/cnt !! P: range jump value
      deltal = deltal/cnt !! L: phase jump value
      ratio = dabs(deltap/deltal)
      if (dabs(jumpval(1)) .gt. VLIGHT*1E-3 - 15) then
        sec = (jumpval(1) + jumpval(2))/2.d0/VLIGHT*1E3   ! ms
        if (dabs(sec - nint(sec)) .lt. 1E-5) then
          sec = nint(sec) * VLIGHT * 1E-3   ! m
        else
          sec = 0
        end if
      else if (dabs(jumpval(1)) .gt. 1E-7 * VLIGHT - 15 .and. dabs(jumpval(1)) .lt. 1E-5 * VLIGHT + 15) then
        sec = (jumpval(1) + jumpval(2))/2.d0/VLIGHT*1E6   ! us
        sec = sec * VLIGHT * 1E-6   ! m
      else
        sec = 0
      end if
      ! inverse fix
      if (ratio < 1.d0) then  ! phase jump
        jumpflag = 1
        jumpsum(1) = jumpsum(1) + sec !! range
      else
        jumpflag = 2
        jumpsum(2) = jumpsum(2) + sec !! phase
      end if
      if (lfnjmp .eq. 0) then
        lfnjmp = get_valid_unit(10)
        open (lfnjmp, file='.'//stanam//'.jmp', status='replace')
      end if
      write (lfnjmp, '(f8.1,2x,i1,x,f28.14)') ti(iepo), jumpflag, sec
    end if

    do ichn = 1, OB%nprn
      jj = OB%prn(ichn)
      j = pointer_string(MAXSAT, prn0, jj)
      if (j .eq. 0) cycle
      if (nsat .lt. j) nsat = j
      if (use_brdeph) then
        ieph = 0
        call elevation(neph, ephem, jj, jd0, ti(iepo), x(iepo), y(iepo), z(iepo), elev, range, dtsat, vv, ieph, .false.)
        if (range .lt. 0.d0) then
          flagall(iepo, j) = set_flag(0, 'pcbad')
          cycle
        end if
        obs(iepo, j, 7) = elev
      end if
      i0 = index(GNSS_PRIO, jj(1:1))
      if ('R' .eq. jj(1:1)) then
        read (OB%prn(ichn), '(1x,i2)') prn_int
        lam1(i0) = VLIGHT/FREQ1_R(OB%glschn(prn_int))
        lam2(i0) = VLIGHT/FREQ2_R(OB%glschn(prn_int))
        lamw(i0) = VLIGHT/(FREQ1_R(OB%glschn(prn_int)) - FREQ2_R(OB%glschn(prn_int)))
      end if
      if (all(dabs(OB%obs(ichn, 1:4)) .gt. 1.d-3)) then
        OB%obs(ichn, 1) = OB%obs(ichn, 1) + jumpsum(2)/lam1(i0)   !! L1
        OB%obs(ichn, 2) = OB%obs(ichn, 2) + jumpsum(2)/lam2(i0)   !! L2
        OB%obs(ichn, 3:4) = OB%obs(ichn, 3:4) + jumpsum(1)        !! P1 & P2
        !! Melbourne-Wubbena (N1 - N2)
        obs(iepo, j, 2) = (OB%obs(ichn, 1) - OB%obs(ichn, 2)) - &
                   (rg(i0)*OB%obs(ichn, 3) + OB%obs(ichn, 4))/(1.d0 + rg(i0))/lamw(i0)
        !! ionosphere-free (pp37 LC)
        obs(iepo, j, 3) = c1(i0)*OB%obs(ichn, 1) + rg(i0)*c2(i0)*OB%obs(ichn, 2)
        if (use_brdeph) then
          ieph = 0
          call elevation(neph, ephem, jj, jd0, ti(iepo), x(iepo), y(iepo), z(iepo), elev, range, dtsat, vv, ieph, .false.)
          v(iepo, j) = vv/PI*180.0
          !! (sit-sat distance from PC) - (sit-sat distance from broadcast) to check recv clock
          obs(iepo, j, 6) = c1(i0)*OB%obs(ichn, 3) + c2(i0)*OB%obs(ichn, 4) - (range - dtsat)
          ! ******************************************* !
          !              check elevation                !
          ! ******************************************* !
          if (elev .lt. cutoff_elevation) then
            flagall(iepo, j) = set_flag(flagall(iepo, j), 'lowele')
          end if
        end if
        ! **************************************************************** !
        !                        check Melbourne-Wubbenamw                 !
        ! **************************************************************** !
        if (ilast(j) .ne. 0) then
          nwdif = dabs(obs(iepo, j, 2) - obs(ilast(j), j, 2))
          if (nwdif .gt. nwlimit) then
            flagall(iepo, j) = set_flag(flagall(iepo, j), 'lwjump')
          end if
        end if
        ilast(j) = iepo
      else
        flagall(iepo, j) = set_flag(0, 'no4')
      end if
    end do
    iepo = iepo + 1
    tobs = tobs + interval
  end do
  close (iunit)
  if (lfnjmp .ne. 0) close (lfnjmp)
!
!! check receiver clock
  if (.not. use_brdeph .or. pclimit .eq. 0.d0) return
  iepo = 1
  do while (iepo .le. nepo)
    do i0 = 1, MAXSYS
      nused = 0
      pmb = 0.d0
      flg = 0
      wpmb = 0.d0
      it = 0
      ij = 0
      do j = 1, i0-1
        ij = ij + MAXSAT_SYS(j)
      end do
      do j = ij+1, ij+MAXSAT_SYS(i0)
        if ('C' .eq. prn0(j)(1:1)) then
          read (prn0(j) (2:3), '(i2)') prn_int
          if (prn_int .le. 5) cycle
        end if
        if (istrue(flagall(iepo, j), 'ok')) then
          nused = nused + 1
          pmb(nused) = obs(iepo, j, 6)
          flg(nused) = 0
          wpmb(nused) = 1.d0
          it(nused) = j
        end if
      end do
      if (nused .eq. 0) cycle
      call get_wgt_mean(.true., pmb, flg, wpmb, nused, 30.d0, k, dt, rms, sig)
      do while (nused - k .gt. 2 .and. rms .gt. 30.d0)
        j = k
        call sign_robust(nused, pmb, flg, 60.d0, k)
        if (k .eq. j) exit
        call get_wgt_mean(.false., pmb, flg, wpmb, nused, 30.d0, k, dt, rms, sig)
      end do
      if (k .gt. 0) then
        do j = 1, nused
          if (flg(j) .le. 1) then
            obs(iepo, it(j), 6) = pmb(j) - dt
          end if
          if (flg(j) .ge. 2 .or. dabs(obs(iepo, it(j), 6)) .gt. pclimit) then
            if (dabs(obs(iepo, it(j), 6)) .gt. 2*1.d5) then
              flagall(iepo, it(j)) = set_flag(flagall(iepo, it(j)), 'pc1ms')
            else
              flagall(iepo, it(j)) = set_flag(flagall(iepo, it(j)), 'pcbad')
            end if
          end if
        end do
      end if
    end do
    iepo = iepo + 1
  end do
!
!! return time tag in seconds
  j = modified_julday(HD%t0(3), HD%t0(2), HD%t0(1))
  t_first_in_rinex = HD%t0(4)*3600.d0 + HD%t0(5)*60.d0 + HD%t0s + (j - jd0)*86400.d0
  t_last_in_rinex = ti(nepo)
  return

100 continue
  write (*, '(a,5(1x,a,2i1))') '***ERROR(read_rinex_file): invalid frequency number:', &
    'G', idxfrq(index(GNSS_PRIO, 'G'), 1:2), &
    'R', idxfrq(index(GNSS_PRIO, 'R'), 1:2), &
    'E', idxfrq(index(GNSS_PRIO, 'E'), 1:2), &
    'C', idxfrq(index(GNSS_PRIO, 'C'), 1:2), &
    'J', idxfrq(index(GNSS_PRIO, 'J'), 1:2)
  call exit(1)
end subroutine
