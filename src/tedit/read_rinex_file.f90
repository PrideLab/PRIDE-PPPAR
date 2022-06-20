!
!! read_rinex_file.f90
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
!! flnrnx ------------ rinex O_file
! tstart,sstart ----- start time: tstart(1):year, tstart(2):month, tstart(3):day, tstart(4):hour, tstart(5):minute, sstart:second.
! sesstion_length---- sesstion time length
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
subroutine read_rinex_file(flnrnx, tstart, sstart, interval, &
                           check_pc, pclimit, cutoff_elevation, use_brdeph, &
                           neph, ephem, stanam, x, y, z, t_first_in_rinex, t_last_in_rinex, v, &
                           nepo, nsat, jd0, nobs, flagall, ti, ts, obs, bias)
  implicit none
  include '../header/const.h'
  include '../header/rnxobs.h'
  include '../header/brdeph.h'

  logical*1 check_pc, use_brdeph
  integer*4 neph, tstart(1:*), nepo
  real*4 v(nepo, MAXSAT)
  real*8 sstart, x(nepo), y(nepo), z(nepo), interval, lglimit, nwlimit, cutoff_elevation, pclimit
  character*(*) flnrnx, stanam
  type(brdeph) ephem(1:*)

  integer*4 nsat, flagall(nepo, MAXSAT), jd0, nobs(MAXSAT), iday
  real*8 ti(nepo), ts(nepo), obs(nepo, MAXSAT, 6), bias(MAXSAT, MAXTYP)

! local
  integer*4 nprn, iepo, isat, j, k, iunit, iunit_next, ichn, ieph, nused, ilast(MAXSAT)
  character*3 jj
  character*3 prn(MAXSAT), prn0(MAXSAT)
  integer*4 ierr
  real*8 tobs, nwdif, lgdif, pgdif, elev, range, dtsat, dt, dwnd
  ! G
  real*8 g_G, lambda1_G, lambda2_G, lambdaw_G, c1_G, c2_G
  ! R
  integer*4 frequency_glo_nu,prn_int,i
  character*3 prn_str
  real*8 :: FREQ1_R(-50:50),FREQ2_R(-50:50)
  real*8 g_R,lambda1_R(-50:50),lambda2_R(-50:50),lambdaw_R(-50:50),c1_R,c2_R

  ! E
  real*8 g_E,lambda1_E,lambda2_E,lambdaw_E,c1_E,c2_E
  
  ! C
  real*8 g_C,lambda1_C,lambda2_C,lambdaw_C,c1_C,c2_C
  
  ! J
  real*8 g_J,lambda1_J,lambda2_J,lambdaw_J,c1_J,c2_J

  real*8 t_first_in_rinex, t_last_in_rinex, vv
  type(rnxhdr) HD
  type(rnxobr) OB

! for range check
  integer*4 flg(MAXSAT), it(MAXSAT)
  real*8 pmb(MAXSAT), wpmb(MAXSAT), rms, sig

! receiver clock jump check
  real*8        jumpsum(2), jumpval(2), rangeval(2, MAXSAT), deltap, deltal
  real*8        tmp(2), ratio, sec, obsval(0:2, 4, MAXSAT), chkval(2, MAXSAT)
  integer*4     lfnjmp, nvalid, njump, cnt, jumpflag  ! 1:range, 2:phase

! funtion called
  logical*1 istrue
  integer*4 modified_julday, set_flag, get_valid_unit, pointer_string
  real*8 timdif
  
  call prn_matbld(prn0)
  call frequency_glonass(FREQ1_R,FREQ2_R)
  ! G
  g_G=FREQ1_G/FREQ2_G
  c1_G=g_G*g_G/(g_G*g_G-1.d0)
  c2_G=-1.d0/(g_G*g_G-1.d0)
  lambda1_G=VLIGHT/FREQ1_G
  lambda2_G=VLIGHT/FREQ2_G
  lambdaw_G=VLIGHT/(FREQ1_G-FREQ2_G)
  ! R
  g_R=9.0d0/7.0d0
  c1_R=g_R*g_R/(g_R*g_R-1.d0)
  c2_R=-1.d0/(g_R*g_R-1.d0)
  ! E
  g_E=FREQ1_E/FREQ2_E
  c1_E=g_E*g_E/(g_E*g_E-1.d0)
  c2_E=-1.d0/(g_E*g_E-1.d0)
  lambda1_E=VLIGHT/FREQ1_E
  lambda2_E=VLIGHT/FREQ2_E
  lambdaw_E=VLIGHT/(FREQ1_E-FREQ2_E)
  ! C
  g_C=FREQ1_C/FREQ2_C
  c1_C=g_C*g_C/(g_C*g_C-1.d0)
  c2_C=-1.d0/(g_C*g_C-1.d0)
  lambda1_C=VLIGHT/FREQ1_C
  lambda2_C=VLIGHT/FREQ2_C
  lambdaw_C=VLIGHT/(FREQ1_C-FREQ2_C)
  ! J
  g_J=FREQ1_J/FREQ2_J
  c1_J=g_J*g_J/(g_J*g_J-1.d0)
  c2_J=-1.d0/(g_J*g_J-1.d0)
  lambda1_J=VLIGHT/FREQ1_J
  lambda2_J=VLIGHT/FREQ2_J
  lambdaw_J=VLIGHT/(FREQ1_J-FREQ2_J)

  lglimit = 50.d0
  nwlimit = 50.d0

  lfnjmp = 0
  jumpsum = 0.d0
  obsval = 0.d0

  bias = 0.d0

  dwnd = min(interval/10.d0, 0.3d0)
! open rinex observation file
  iunit = 10
  open (iunit, file=flnrnx, form='FORMATTED', status='OLD', iostat=ierr)
  if (ierr .ne. 0) then
!  write(*,'(a)') '***ERROR(read_rinex_file): open file error, '//flnrnx
    call exit(1)
  endif

! read head of rinex file
  call rdrnxoh(iunit, HD, ierr)
  if (ierr .ne. 0) then
    call exit(1)
  endif

! start time of first session
  if (tstart(2) .eq. 0) then
    do j = 1, 5
      tstart(j) = HD%t0(j)
    enddo
    sstart = HD%t0s
  endif
  call yr2year(tstart(1))

!! jd0:  julday
!! tobs: fractional julday
  jd0 = modified_julday(tstart(3), tstart(2), tstart(1))
  tobs = tstart(4)*3600.d0 + tstart(5)*60.d0 + sstart
  call read_glschn(jd0, tobs, OB%glschn)

! initialize
  nsat = 0
  do isat = 1, MAXSAT
    ilast(isat) = 0
  enddo

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
      enddo
    enddo
  enddo
!
! loop over session and epoch
  iday = 1
  iepo = 1
  do while (iepo .le. nepo)
    if (ierr .ne. 0 .and. ierr .ne. 2) cycle
    nprn = 0
    prn = ''
    !
    !! try reading the next rnxo file
    if (ti(iepo) .ge. iday*86400.d0) then
      iday = iday + 1
      call next_rinex(iunit, iunit_next)
      if (iunit_next .ne. 0) then
        close(iunit)
        iunit = iunit_next
      endif
      ierr = 0
    endif
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
    if (ierr .ne. 2) then
      if (HD%ver .eq. 3) then
        call rdrnxoi3(iunit, jd0, tobs, dwnd, nprn, prn, HD, OB, bias, ierr)
      else
        call rdrnxoi2(iunit, jd0, tobs, dwnd, nprn, prn, HD, OB, bias, ierr)
      endif
    endif
    if (ierr .eq. 0) then
      obsval(0,:,:) = obsval(1,:,:)
      obsval(1,:,:) = obsval(2,:,:)
      obsval(2,:,:) = 0.d0
      ti(iepo) = tobs
      if (OB%nprn .ne. 0) then
        ti(iepo) = ti(iepo) - timdif(jd0, tobs, OB%jd, OB%tsec)
      endif
      do ichn = 1, OB%nprn
        jj = OB%prn(ichn)
        j=pointer_string(MAXSAT,prn0,jj)
        if (j .eq. 0) cycle
! nsat : largest PRN for a satellite
        if (nsat .lt. j) nsat = j
        if (use_brdeph) then
          ieph = 0
! range --------- the distance between satellite and reciver
! dtsat --------- clock correction, units: meter
! vv    --------- true anomaly
          call elevation(neph, ephem, jj, jd0, ti(iepo), x(iepo), y(iepo), z(iepo), elev, range, dtsat, vv, ieph, .false.)
          if (range .lt. 0.d0) then
            flagall(iepo, j) = set_flag(0, 'pcbad')
            cycle
          endif
        endif
        if (dabs(OB%obs(ichn, 1)) .gt. 1.d-3 .and. dabs(OB%obs(ichn, 2)) .gt. 1.d-3 .and. &
            dabs(OB%obs(ichn, 3)) .gt. 1.d-3 .and. dabs(OB%obs(ichn, 4)) .gt. 1.d-3) then
          flagall(iepo, j) = 0
          obsval(2,1,j) = OB%obs(ichn, 3) ! P1
          obsval(2,3,j) = OB%obs(ichn, 4) ! P2
          if(OB%prn(ichn)(1:1) .eq. 'G')then
            obsval(2,2,j) = OB%obs(ichn, 1)*lambda1_G ! L1
            obsval(2,4,j) = OB%obs(ichn, 2)*lambda2_G ! L2
          elseif(OB%prn(ichn)(1:1) .eq. 'R')then
            read(OB%prn(ichn),'(1x,i2)') prn_int
            frequency_glo_nu=OB%glschn(prn_int)
            lambda1_R(frequency_glo_nu)=VLIGHT/FREQ1_R(frequency_glo_nu)
            lambda2_R(frequency_glo_nu)=VLIGHT/FREQ2_R(frequency_glo_nu)
            obsval(2,2,j) = OB%obs(ichn, 1)*lambda1_R(frequency_glo_nu) ! L1
            obsval(2,4,j) = OB%obs(ichn, 2)*lambda2_R(frequency_glo_nu) ! L2
          elseif(OB%prn(ichn)(1:1) .eq. 'E')then
            obsval(2,2,j) = OB%obs(ichn, 1)*lambda1_E ! L1
            obsval(2,4,j) = OB%obs(ichn, 2)*lambda2_E ! L2
          elseif(OB%prn(ichn)(1:1) .eq. 'C')then
            obsval(2,2,j) = OB%obs(ichn, 1)*lambda1_C ! L1
            obsval(2,4,j) = OB%obs(ichn, 2)*lambda2_C ! L2
          elseif(OB%prn(ichn)(1:1) .eq. 'J')then
            obsval(2,2,j) = OB%obs(ichn, 1)*lambda1_J ! L1
            obsval(2,4,j) = OB%obs(ichn, 2)*lambda2_J ! L2
          endif
          rangeval(1,j) = rangeval(2,j)
          rangeval(2,j) = range
          if(OB%prn(ichn)(1:1) .eq. 'G')then
            ! geometry-free : ionosphere observations
            obs(iepo, j, 1) = (lambda1_G*OB%obs(ichn, 1) - lambda2_G*OB%obs(ichn, 2))/(lambda2_G - lambda1_G)
          elseif(OB%prn(ichn)(1:1) .eq. 'R')then
            read(OB%prn(ichn),'(1x,i2)') prn_int
            frequency_glo_nu=OB%glschn(prn_int)
            lambda1_R(frequency_glo_nu)=VLIGHT/FREQ1_R(frequency_glo_nu)
            lambda2_R(frequency_glo_nu)=VLIGHT/FREQ2_R(frequency_glo_nu)
            obs(iepo, j, 1) = (lambda1_R(frequency_glo_nu)*OB%obs(ichn, 1) - lambda2_R(frequency_glo_nu)*OB%obs(ichn, 2))&
                             /(lambda2_R(frequency_glo_nu) - lambda1_R(frequency_glo_nu))
          elseif(OB%prn(ichn)(1:1) .eq. 'E')then
            obs(iepo, j, 1) = (lambda1_E*OB%obs(ichn, 1) - lambda2_E*OB%obs(ichn, 2))/(lambda2_E - lambda1_E)
          elseif(OB%prn(ichn)(1:1) .eq. 'C')then
            obs(iepo, j, 1) = (lambda1_C*OB%obs(ichn, 1) - lambda2_C*OB%obs(ichn, 2))/(lambda2_C - lambda1_C)
          elseif(OB%prn(ichn)(1:1) .eq. 'J')then
            obs(iepo, j, 1) = (lambda1_J*OB%obs(ichn, 1) - lambda2_J*OB%obs(ichn, 2))/(lambda2_J - lambda1_J)
          endif
          ! **************************************************************** !
          !                       check geometry-free(lg)                    !
          ! **************************************************************** !
          if (ilast(j) .ne. 0) then
            lgdif = dabs(obs(iepo, j, 1) - obs(ilast(j), j, 1))/(ti(iepo) - ti(ilast(j)))
            if (lgdif .gt. lglimit) then
              flagall(iepo, j) = set_flag(flagall(iepo, j), 'lgjump')
            endif
          endif
        else
          flagall(iepo, j) = set_flag(0, 'no4')
        endif
      enddo

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
        j=pointer_string(MAXSAT,prn0,jj)
        if (j .eq. 0) cycle
        if (.not.istrue(flagall(iepo,j), 'ok')) cycle
        if (dabs(obsval(0, 1, j)).lt.1.d-3 .or. dabs(obsval(1, 1, j)).lt.1.d-3 .or. dabs(obsval(2, 1, j)).lt.1.d-3) cycle
        nvalid = nvalid + 1
        chkval(1,j) = obsval(2,1,j)-obsval(1,1,j) - (obsval(2,2,j)-obsval(1,2,j)) ! f1
        chkval(2,j) = obsval(2,3,j)-obsval(1,3,j) - (obsval(2,4,j)-obsval(1,4,j)) ! f2
        if (dabs(chkval(1, j)).gt.20.d0 .and. dabs(chkval(2, j)).gt.20.d0) then ! 0.1 us jump: 30.d0
          njump = njump + 1
          if (rangeval(1, j).gt.1.d0 .and. rangeval(2, j).gt.1.d0) then
            cnt = cnt + 1
            tmp(1) = (obsval(2,1,j)-obsval(1,1,j) + obsval(2,3,j)-obsval(1,3,j))/2.d0 - (rangeval(2,j) - rangeval(1,j))
            tmp(2) = (obsval(2,2,j)-obsval(1,2,j) + obsval(2,4,j)-obsval(1,4,j))/2.d0 - (rangeval(2,j) - rangeval(1,j))
            deltap = deltap + tmp(1)
            deltal = deltal + tmp(2)
          endif
          jumpval(1) = jumpval(1) + chkval(1,j)
          jumpval(2) = jumpval(2) + chkval(2,j)
        endif
      enddo

      if (nvalid.ne.0 .and. nvalid.eq.njump) then
        jumpval = jumpval/(nvalid)
        deltap = deltap/cnt  ! P: range jump value
        deltal = deltal/cnt  ! L: phase jump value
        ratio = dabs(deltap/deltal)
        if (dabs(jumpval(1)) .gt. vlight*1E-3-15) then
          sec = (jumpval(1) + jumpval(2))/2.d0/vlight*1E3 ! ms
          if (dabs(sec - nint(sec)) .lt. 1E-5) then
            sec = nint(sec)*1E-3*vlight ! m
          else
            sec = 0
          endif
        else if (dabs(jumpval(1)) .gt. 1E-7*vlight-15 .and. dabs(jumpval(1)) .lt. 1E-5*vlight+15) then
          sec = (jumpval(1)+jumpval(2))/2.d0/vlight*1E6 ! us
          sec = sec*1E-6*vlight ! m
        else
          sec = 0
        endif
        ! inverse fix
        if(ratio < 1.d0) then ! phase jump
          jumpflag = 1
          jumpsum(1) = jumpsum(1) + sec!/1.d3*vlight ! range
        else
          jumpflag = 2
          jumpsum(2) = jumpsum(2) + sec!/1.d3*vlight ! phase
        end if
        if (lfnjmp.eq.0) then
          lfnjmp = get_valid_unit(10)
          open(lfnjmp,file='.'//stanam//'.jmp',status='replace')
        end if
        write(lfnjmp,'(f8.1,2x,i1,x,f28.14)') ti(iepo), jumpflag, sec
      end if

      do ichn = 1, OB%nprn
        jj=OB%prn(ichn)
        j=pointer_string(MAXSAT,prn0,jj)
        if (j .eq. 0) cycle
        if (nsat .lt. j) nsat = j
        if (use_brdeph) then
          ieph = 0
          call elevation(neph, ephem, jj, jd0, ti(iepo), x(iepo), y(iepo), z(iepo), elev, range, dtsat, vv, ieph, .false.)
          if (range .lt. 0.d0) then
            flagall(iepo, j) = set_flag(0, 'pcbad')
            cycle
          endif
        endif
        if (dabs(OB%obs(ichn, 1)) .gt. 1.d-3 .and. dabs(OB%obs(ichn, 2)) .gt. 1.d-3 .and. &
            dabs(OB%obs(ichn, 3)) .gt. 1.d-3 .and. dabs(OB%obs(ichn, 4)) .gt. 1.d-3) then
          OB%obs(ichn, 3:4) = OB%obs(ichn, 3:4) + jumpsum(1) ! P1 & P2
          if(OB%prn(ichn)(1:1) .eq. 'G')then
            OB%obs(ichn, 1) = OB%obs(ichn, 1) + jumpsum(2)/lambda1_G
            OB%obs(ichn, 2) = OB%obs(ichn, 2) + jumpsum(2)/lambda2_G
          elseif(OB%prn(ichn)(1:1) .eq. 'R')then
            read(OB%prn(ichn),'(1x,i2)') prn_int
            frequency_glo_nu=OB%glschn(prn_int)
            lambda1_R(frequency_glo_nu)=VLIGHT/FREQ1_R(frequency_glo_nu)
            lambda2_R(frequency_glo_nu)=VLIGHT/FREQ2_R(frequency_glo_nu)
            OB%obs(ichn, 1) = OB%obs(ichn, 1) + jumpsum(2)/lambda1_R(frequency_glo_nu)
            OB%obs(ichn, 2) = OB%obs(ichn, 2) + jumpsum(2)/lambda2_R(frequency_glo_nu)
          elseif(OB%prn(ichn)(1:1) .eq. 'E')then
            OB%obs(ichn, 1) = OB%obs(ichn, 1) + jumpsum(2)/lambda1_E
            OB%obs(ichn, 2) = OB%obs(ichn, 2) + jumpsum(2)/lambda2_E
          elseif(OB%prn(ichn)(1:1) .eq. 'C')then
            OB%obs(ichn, 1) = OB%obs(ichn, 1) + jumpsum(2)/lambda1_C
            OB%obs(ichn, 2) = OB%obs(ichn, 2) + jumpsum(2)/lambda2_C
          elseif(OB%prn(ichn)(1:1) .eq. 'J')then
            OB%obs(ichn, 1) = OB%obs(ichn, 1) + jumpsum(2)/lambda1_J
            OB%obs(ichn, 2) = OB%obs(ichn, 2) + jumpsum(2)/lambda2_J
          endif
          if(OB%prn(ichn)(1:1) .eq. 'G')then
            ! Melbourne-Wubbena (N1 - N2)
            obs(iepo, j, 2) = OB%obs(ichn, 1) - OB%obs(ichn, 2) - (g_G*OB%obs(ichn, 3) + OB%obs(ichn, 4))/(1.d0 + g_G)/lambdaw_G
            ! ionosphere-free (pp37 LC)
            obs(iepo, j, 3) = c1_G*OB%obs(ichn, 1) + g_G*c2_G*OB%obs(ichn, 2)
          elseif(OB%prn(ichn)(1:1) .eq. 'R')then
            read(OB%prn(ichn),'(1x,i2)') prn_int
            frequency_glo_nu=OB%glschn(prn_int)
            lambda1_R(frequency_glo_nu)=VLIGHT/FREQ1_R(frequency_glo_nu)
            lambda2_R(frequency_glo_nu)=VLIGHT/FREQ2_R(frequency_glo_nu)
            lambdaw_R(frequency_glo_nu)=VLIGHT/(FREQ1_R(frequency_glo_nu)-FREQ2_R(frequency_glo_nu))
            obs(iepo, j, 2) = OB%obs(ichn, 1) - OB%obs(ichn, 2) - (g_R*OB%obs(ichn, 3) + OB%obs(ichn, 4))/(1.d0 + g_R)&
                              /lambdaw_R(frequency_glo_nu)
            obs(iepo, j, 3) = c1_R*OB%obs(ichn, 1) + g_R*c2_R*OB%obs(ichn, 2)
          elseif(OB%prn(ichn)(1:1) .eq. 'E')then
            obs(iepo, j, 2) = OB%obs(ichn, 1) - OB%obs(ichn, 2) - (g_E*OB%obs(ichn, 3) + OB%obs(ichn, 4))/(1.d0 + g_E)/lambdaw_E
            obs(iepo, j, 3) = c1_E*OB%obs(ichn, 1) + g_E*c2_E*OB%obs(ichn, 2)
          elseif(OB%prn(ichn)(1:1) .eq. 'C')then
            obs(iepo, j, 2) = OB%obs(ichn, 1) - OB%obs(ichn, 2) - (g_C*OB%obs(ichn, 3) + OB%obs(ichn, 4))/(1.d0 + g_C)/lambdaw_C
            obs(iepo, j, 3) = c1_C*OB%obs(ichn, 1) + g_C*c2_C*OB%obs(ichn, 2)
          elseif(OB%prn(ichn)(1:1) .eq. 'J')then
            obs(iepo, j, 2) = OB%obs(ichn, 1) - OB%obs(ichn, 2) - (g_J*OB%obs(ichn, 3) + OB%obs(ichn, 4))/(1.d0 + g_J)/lambdaw_J
            obs(iepo, j, 3) = c1_J*OB%obs(ichn, 1) + g_J*c2_J*OB%obs(ichn, 2)
          end if
          if (use_brdeph) then
            ieph = 0
            call elevation(neph, ephem, jj, jd0, ti(iepo), x(iepo), y(iepo), z(iepo), elev, range, dtsat, vv, ieph, .false.)
            v(iepo, j) = vv/PI*180.0
            ! (sit-sat distance from PC)-(sit-sat distance from broadcast) to check recv clock
            if(OB%prn(ichn)(1:1) .eq. 'G')then
              obs(iepo, j, 6) = c1_G*OB%obs(ichn, 3) + c2_G*OB%obs(ichn, 4) - (range - dtsat)
            elseif(OB%prn(ichn)(1:1) .eq. 'R')then
              obs(iepo, j, 6) = c1_R*OB%obs(ichn, 3) + c2_R*OB%obs(ichn, 4) - (range - dtsat)
            elseif(OB%prn(ichn)(1:1) .eq. 'E')then
              obs(iepo, j, 6) = c1_E*OB%obs(ichn, 3) + c2_E*OB%obs(ichn, 4) - (range - dtsat)
            elseif(OB%prn(ichn)(1:1) .eq. 'C')then
              obs(iepo, j, 6) = c1_C*OB%obs(ichn, 3) + c2_C*OB%obs(ichn, 4) - (range - dtsat)
            elseif(OB%prn(ichn)(1:1) .eq. 'J')then
              obs(iepo, j, 6) = c1_J*OB%obs(ichn, 3) + c2_J*OB%obs(ichn, 4) - (range - dtsat)
            end if
            ! ******************************************* !
            !              check elevation                !
            ! ******************************************* !
            if (elev .lt. cutoff_elevation) then
              flagall(iepo, j) = set_flag(flagall(iepo, j), 'lowele')
            endif
          endif
          ! **************************************************************** !
          !                        check Melbourne-Wubbenamw                 !
          ! **************************************************************** !
          if (ilast(j) .ne. 0) then
            nwdif = dabs(obs(iepo, j, 2) - obs(ilast(j), j, 2))
            if (nwdif .gt. nwlimit) then
              flagall(iepo, j) = set_flag(flagall(iepo, j), 'lwjump')
            endif
          endif
          ilast(j) = iepo
        else
          flagall(iepo, j) = set_flag(0, 'no4')
        endif
      enddo
      iepo = iepo + 1
      tobs = tobs + interval
    endif
    if (ierr .eq. 2) then
      iepo = iepo + 1
      tobs = tobs + interval
      flagall(iepo, :) = set_flag(0, 'no4')
    endif
  enddo
  close (10)
  if (lfnjmp .ne. 0) close(lfnjmp)
!
!! check receiver clock
  if(.not.use_brdeph.or.pclimit.eq.0.d0) return
  iepo = 1
  do while (iepo .le. nepo)
!
!! GPS
    nused=0
    do j=1,MAXSAT_G
      if(istrue(flagall(iepo,j),'ok')) then
        nused=nused+1
        pmb(nused) =obs(iepo,j,6)
        flg(nused) =0
        wpmb(nused)=1.d0
        it(nused)  =j
      endif
    enddo
    if(nused.ne.0) then
      call get_wgt_mean(.true.,pmb,flg,wpmb,nused,30.d0,k,dt,rms,sig)
      do while(nused-k.gt.2.and.rms.gt.30.d0)
        j=k
        call sign_robust(nused,pmb,flg,60.d0,k)
        if(k.eq.j) exit
        call get_wgt_mean(.false.,pmb,flg,wpmb,nused,30.d0,k,dt,rms,sig)
      enddo
      if(k.gt.0) then
        do j=1,nused
          if(flg(j).le.1) then
            obs(iepo,it(j),6)=pmb(j)-dt
          endif
          if(flg(j).ge.2.or.dabs(obs(iepo,it(j),6)).gt.pclimit) then
            if(dabs(obs(iepo,it(j),6)).gt.2*1.d5) then
              flagall(iepo,it(j))=set_flag(flagall(iepo,it(j)),'pc1ms')
            else
              flagall(iepo,it(j))=set_flag(flagall(iepo,it(j)),'pcbad')
            endif
            !write(*,'(a15,i6,7x,a1,i2.2,f15.1,1x,f11.1)') ' bad range ',iepo,'G',it(j),pmb(j)-dt,dt
          endif
        enddo
      endif
    endif
!    
!! GLONASS
    nused=0
    do j=MAXSAT_G+1,MAXSAT_G+MAXSAT_R
      if(istrue(flagall(iepo,j),'ok')) then
        nused=nused+1
        pmb(nused) =obs(iepo,j,6)
        flg(nused) =0
        wpmb(nused)=1.d0
        it(nused)  =j
      endif
    enddo
    if(nused.ne.0) then
      call get_wgt_mean(.true.,pmb,flg,wpmb,nused,30.d0,k,dt,rms,sig)
      do while(nused-k.gt.2.and.rms.gt.30.d0)
        j=k
        call sign_robust(nused,pmb,flg,60.d0,k)
        if(k.eq.j) exit
        call get_wgt_mean(.false.,pmb,flg,wpmb,nused,30.d0,k,dt,rms,sig)
      enddo
      if(k.gt.0) then
        do j=1,nused
          if(flg(j).le.1) then
            obs(iepo,it(j),6)=pmb(j)-dt
          endif
          if(flg(j).ge.2.or.dabs(obs(iepo,it(j),6)).gt.pclimit) then
            if(dabs(obs(iepo,it(j),6)).gt.2*1.d5) then
              flagall(iepo,it(j))=set_flag(flagall(iepo,it(j)),'pc1ms')
            else
              flagall(iepo,it(j))=set_flag(flagall(iepo,it(j)),'pcbad')
            endif
            !write(*,'(a15,i6,7x,a1,i2.2,f15.1,1x,f11.1)') ' bad range ',iepo,'R',it(j)-MAXSAT_G,pmb(j)-dt,dt
          endif
        enddo
      endif
    endif
!
!! Galileo
    nused=0
    do j=MAXSAT_G+MAXSAT_R+1,MAXSAT_G+MAXSAT_R+MAXSAT_E
      if(istrue(flagall(iepo,j),'ok')) then
        nused=nused+1
        pmb(nused) =obs(iepo,j,6)
        flg(nused) =0
        wpmb(nused)=1.d0
        it(nused)  =j
      endif
    enddo
    if(nused.ne.0) then
      call get_wgt_mean(.true.,pmb,flg,wpmb,nused,30.d0,k,dt,rms,sig)
      do while(nused-k.gt.2.and.rms.gt.30.d0)
        j=k
        call sign_robust(nused,pmb,flg,60.d0,k)
        if(k.eq.j) exit
        call get_wgt_mean(.false.,pmb,flg,wpmb,nused,30.d0,k,dt,rms,sig)
      enddo
      if(k.gt.0) then
        do j=1,nused
          if(flg(j).le.1) then
            obs(iepo,it(j),6)=pmb(j)-dt
          endif
          if(flg(j).ge.2.or.dabs(obs(iepo,it(j),6)).gt.pclimit) then
            if(dabs(obs(iepo,it(j),6)).gt.2*1.d5) then
              flagall(iepo,it(j))=set_flag(flagall(iepo,it(j)),'pc1ms')
            else
              flagall(iepo,it(j))=set_flag(flagall(iepo,it(j)),'pcbad')
            endif
            !write(*,'(a15,i6,7x,a1,i2.2,f15.1,1x,f11.1)') ' bad range ',iepo,'E',it(j)-MAXSAT_G-MAXSAT_R,pmb(j)-dt,dt
          endif
        enddo
      endif
    endif
!
!! BDS
    nused=0
    do j=MAXSAT_G+MAXSAT_R+MAXSAT_E+1,MAXSAT_G+MAXSAT_R+MAXSAT_E+MAXSAT_C
      read(prn0(j)(2:3),'(i2)') prn_int
      if(prn_int .le. 5) cycle
      if(istrue(flagall(iepo,j),'ok')) then
        nused=nused+1
        pmb(nused) =obs(iepo,j,6)
        flg(nused) =0
        wpmb(nused)=1.d0
        it(nused)  =j
      endif
    enddo
    if(nused.ne.0) then
      call get_wgt_mean(.true.,pmb,flg,wpmb,nused,30.d0,k,dt,rms,sig)
      do while(nused-k.gt.2.and.rms.gt.30.d0)
        j=k
        call sign_robust(nused,pmb,flg,60.d0,k)
        if(k.eq.j) exit
        call get_wgt_mean(.false.,pmb,flg,wpmb,nused,30.d0,k,dt,rms,sig)
      enddo
      if(k.gt.0) then
        do j=1,nused
          if(flg(j).le.1) then
            obs(iepo,it(j),6)=pmb(j)-dt
          endif
          if(flg(j).ge.2.or.dabs(obs(iepo,it(j),6)).gt.pclimit) then
            if(dabs(obs(iepo,it(j),6)).gt.2*1.d5) then
              flagall(iepo,it(j))=set_flag(flagall(iepo,it(j)),'pc1ms')
            else
              flagall(iepo,it(j))=set_flag(flagall(iepo,it(j)),'pcbad')
            endif
            !write(*,'(a15,i6,7x,a1,i2.2,f15.1,1x,f11.1)') ' bad range ',iepo,'C', &
            !      it(j)-MAXSAT_G-MAXSAT_R-MAXSAT_E,pmb(j)-dt,dt
          endif
        enddo
      endif
    endif
!
!! QZSS
    nused=0
    do j=MAXSAT_G+MAXSAT_R+MAXSAT_E+MAXSAT_C+1,MAXSAT_G+MAXSAT_R+MAXSAT_E+MAXSAT_C+MAXSAT_J
      if(istrue(flagall(iepo,j),'ok')) then
        nused=nused+1
        pmb(nused) =obs(iepo,j,6)
        flg(nused) =0
        wpmb(nused)=1.d0
        it(nused)  =j
      endif
    enddo
    if(nused.ne.0) then
      call get_wgt_mean(.true.,pmb,flg,wpmb,nused,30.d0,k,dt,rms,sig)
      do while(nused-k.gt.2.and.rms.gt.30.d0)
        j=k
        call sign_robust(nused,pmb,flg,60.d0,k)
        if(k.eq.j) exit
        call get_wgt_mean(.false.,pmb,flg,wpmb,nused,30.d0,k,dt,rms,sig)
      enddo
      if(k.gt.0) then
        do j=1,nused
          if(flg(j).le.1) then
            obs(iepo,it(j),6)=pmb(j)-dt
          endif
          if(flg(j).ge.2.or.dabs(obs(iepo,it(j),6)).gt.pclimit) then
            if(dabs(obs(iepo,it(j),6)).gt.2*1.d5) then
              flagall(iepo,it(j))=set_flag(flagall(iepo,it(j)),'pc1ms')
            else
              flagall(iepo,it(j))=set_flag(flagall(iepo,it(j)),'pcbad')
            endif
            !write(*,'(a15,i6,7x,a1,i2.2,f15.1,1x,f11.1)') ' bad range ',iepo,'J', &
            !      it(j)-MAXSAT_G-MAXSAT_R-MAXSAT_E-MAXSAT_C,pmb(j)-dt,dt
          endif
        enddo
      endif
    endif
    iepo = iepo + 1
  enddo

! time tag in second
  j = modified_julday(HD%t0(3), HD%t0(2), HD%t0(1))
  t_first_in_rinex = HD%t0(4)*3600.d0 + HD%t0(5)*60.d0 + HD%t0s + (j - jd0)*86400.d0
  t_last_in_rinex = ti(nepo)

  return
end
