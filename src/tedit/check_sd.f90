!
!! check_sd.f90
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
!! purpose  : check single-difference LC
!! parameter:
!!    input : lfnsd --   0 : debug-sd = true
!!                      -1 : debug-sd = false
!!            ndgr  --   2
!!            niter --   3
!!            mepo  --   20 :arc length for checking
!!            nstep --   8  :step for next checking
!!          lclimit --
!!
!
subroutine check_sd(lfnsd, neph, ephem, ndgr, niter, mepo, nstep, x, y, z, interval, lclimit, &
                    nepo, nsat, jd0, nobs, flagall, ti, ts, obs)
  implicit none
  include '../header/const.h'
  include '../header/brdeph.h'
  include 'data_flag.h'

! common
  integer*4     idxfrq(MAXSYS, 2)
  common        idxfrq
! parameter
  integer*4     lfnsd
  integer*4     neph
  type(brdeph)  ephem(1:*)
  integer*4     ndgr, niter, mepo, nstep
  real*8        x(nepo), y(nepo), z(nepo)
  real*8        interval, lclimit
  integer*4     nepo, nsat
  integer*4     jd0
  integer*4     nobs(MAXSAT)
  integer*4     flagall(nepo, MAXSAT)
  real*8        ti(nepo), ts(nepo)
  real*8        obs(nepo, MAXSAT, 7)
! local
  logical*1     found, ref_ok, lwrite
  integer*4     i, i0, i1, j, k
  integer*4     isys, iG, istart, istop, iepo, kobs, ierr, jeph, jepo, isat, iref
  integer*4     ngood, nsd, njump
  integer*4     prn_int
  integer*4     ichecked(nepo, MAXSAT), jump(nepo*10, 4)
  integer*4     nobs_epoch(nepo), nobsprn(MAXSAT), nobsflg(MAXSAT)
  integer*4     flglg(nepo), flglc(nepo), flglgref(nepo)
  integer*4     ilc(nepo), ilg(nepo)
  real*8        lc(nepo), lg(nepo), range_ref(nepo), range(nepo)
  real*8        vlc(nepo), vlg(nepo), vlgref(nepo)
  real*8        rmslg(MAXSAT), rmslc(MAXSAT)
  real*8        dtsat, elev, dummy, a0
  real*8        lam1(MAXSYS)
  real*8        FREQ1_R(-50:50), FREQ2_R(-50:50)
  integer*4     glschn(MAXSAT_R)
  character*3   prn0(MAXSAT)
  logical*1     lonlyG, lfound
! function called
  logical*1     istrue
  integer*4     set_flag

!
!! initialization
  lwrite = .true.
  if (lfnsd .eq. 0) lfnsd = 6
  if (lfnsd .lt. 0) lwrite = .false.
  call read_glschn(jd0, ti(1), glschn)
  iG = index(GNSS_PRIO, "G")
  do i0 = 1, MAXSYS
    i1 = idxfrq(i0, 1)
    lam1(i0) = VLIGHT/FREQ_SYS(i1, i0)
  end do
  call frequency_glonass(FREQ1_R, FREQ2_R)
  call prn_matbld(prn0)
  ierr = 0
!! find out first and last epoch with enough data
!! nobs_epoch -- # of obs. per epoch, the number of available satellite
  istart = 0    ! start epoch
  istop = 0     ! end epoch
  ngood = 0     ! Effective number of epoch
  do iepo = 1, nepo
    nobs_epoch(iepo) = 0
    do isat = 1, nsat
      if (prn0(isat)(1:1) .eq. 'C') then
        cycle
      end if
      ichecked(iepo, isat) = 10
      if (istrue(flagall(iepo, isat), 'ok')) then
        nobs_epoch(iepo) = nobs_epoch(iepo) + 1    ! the number of available satellite
      end if
    end do
    if (nobs_epoch(iepo) .ge. 2) then
      if (istart .eq. 0) istart = iepo
      istop = iepo
      ngood = ngood + 1
    end if
  end do
  if (istop - istart + 1 .le. ndgr + 2 .or. ngood .le. ndgr + 2) then
! remove all
    return
  end if

  njump = 0
  iepo = istart
  do while (iepo .le. istop)
    ngood = 0

! extract an epoch section from iepo
! mepo = 20; default value, the constant length of check arc is 20 epochs
    i1 = min(iepo + mepo - 1, istop)   ! set the start epoch of check arc
    i0 = max(istart, i1 - mepo + 1)    ! set the end epoch of check arc
    do k = i0, i1
      if (nobs_epoch(k) .ge. 2) ngood = ngood + 1
    end do

! ndgr = 2; default value
    if (ngood .le. ndgr + 2) goto 100

! kobs ------ the number of epoch in check arc
    kobs = i1 - i0 + 1

! poly fit for each satellite in epoch section to choose reference satellite
!-nobsprn : # of observations for each satellite
!-nobsflg : # of flags (amb.) for each satellite
    do isat = 1, nsat
      if (prn0(isat)(1:1) .eq. 'C') then
        cycle
      end if
      nobsprn(isat) = 0
      nobsflg(isat) = 0
      do jepo = 1, kobs
        k = i0 + jepo - 1
        if (istrue(flagall(k, isat), 'ok')) then
          nobsprn(isat) = nobsprn(isat) + 1     ! the number of available epochs for each satellite
        end if
      end do
      if (nobsprn(isat) .gt. ndgr + 2) then
        do jepo = 1, kobs
          k = i0 + jepo - 1
          flglg(jepo) = 10
          if (istrue(flagall(k, isat), 'ok')) then
            flglg(jepo) = 1
            lg(jepo) = obs(k, isat, 1)          ! geometry-free obs value
          end if
          if (istrue(flagall(k, isat), 'good')) then
            flglg(jepo) = 0
          end if
        end do
        call check_for_jump(' LG ', lfnsd, kobs, ti(i0), lg, flglg, ndgr, niter, &
                            2.0d0, 0.40d0, a0, rmslg(isat), vlg, ilg, ierr, interval)
        if (ierr .ne. 0) then
        end if
        do jepo = 1, kobs
          k = i0 + jepo - 1
          if (flglg(jepo) .eq. 1) nobsflg(isat) = nobsflg(isat) + 1
        end do
        if (ierr .ne. 0) rmslg(isat) = 1.d3
      else
        nobsflg(isat) = ngood
        rmslg(isat) = 1.d3
      end if
    end do
!!
!! select satisfactory referece satellte
    iref = 1
    lonlyG = .true.
    lfound = .false.
50  continue
    do isat = 1, nsat
      if (lonlyG .and. prn0(isat)(1:1) .ne. 'G') then !add for referece satellte
        cycle
      end if
      if (nobsprn(isat) .gt. nobsprn(iref)) then
        iref = isat
        lfound = .true.
      else if (nobsprn(isat) .eq. nobsprn(iref)) then
        if (nobsflg(isat) .lt. nobsflg(iref)) then
          iref = isat
          lfound = .true.
        else if (nobsflg(isat) .eq. nobsflg(iref)) then
          if (rmslg(isat) .lt. rmslg(iref)) then
            iref = isat
            lfound = .true.
          end if
        end if
      end if
    end do
    if (iref .eq. 1 .and. nobsprn(iref) .ne. 0) lfound = .true.
    if (.not. lfound .and. lonlyG) then
      lonlyG = .false.
      goto 50
    end if
!!
!! check referece satellite
!!-range_ref : distance from ref. satellites to static station in each epoch
    jeph = 0
    do jepo = 1, kobs
      k = i0 + jepo - 1
      range_ref(jepo) = 0.d0
      flglg(jepo) = 10
      if (istrue(flagall(k, iref), 'ok')) then
        call elevation(neph, ephem, prn0(iref), jd0, ti(k), x(k), y(k), z(k), elev, range_ref(jepo), dtsat, dummy, jeph, .true.)
        range_ref(jepo) = range_ref(jepo) - dtsat
        flglg(jepo) = 1
        if (istrue(flagall(k, iref), 'good')) flglg(jepo) = 0
        lg(jepo) = obs(k, iref, 1)
      end if
    end do
    call check_for_jump(' LGR ', lfnsd, kobs, ti(i0), lg, flglg, ndgr, niter, &
                        2.0d0, 0.40d0, a0, rmslg(iref), vlg, ilg, ierr, interval)
    if (ierr .ne. 0) then
    end if
    do jepo = 1, kobs
      k = i0 + jepo - 1
      flglgref(jepo) = flglg(jepo)
      vlgref(jepo) = 0.d0
      if (ilg(jepo) .ne. 0) then
        vlgref(jepo) = vlgref(ilg(jepo)) + vlg(jepo)
      end if
    end do
!!
!! check other satellites, form single diff LC obs
    do isat = 1, nsat
      if ('C' .eq. prn0(isat)(1:1)) then
        cycle
      end if
      isys = index(GNSS_PRIO, prn0(isat)(1:1))
      if ('R' .eq. prn0(isat)(1:1)) then
        read (prn0(isat), '(1x,i2)') prn_int
        lam1(isys) = VLIGHT/FREQ1_R(glschn(prn_int))
      end if
      if (nobsprn(isat) .gt. ndgr + 2 .and. isat .ne. iref) then
        jeph = 0
        nsd = 0
        do jepo = 1, kobs
          range(jepo) = 0.d0
          k = i0 + jepo - 1
          flglc(jepo) = 10
          if (istrue(flagall(k, isat), 'ok') .and. istrue(flagall(k, iref), 'ok')) then
            call elevation(neph, ephem, prn0(isat), jd0, ti(k), x(k), y(k), z(k), elev, range(jepo), dtsat, dummy, jeph, .true.)
            range(jepo) = range(jepo) - dtsat
            !! LC single-difference = troposphere delay or ...
            lc(jepo) = (obs(k, isat, 3)*lam1(isys) - range(jepo)) - (obs(k, iref, 3)*lam1(iG) - range_ref(jepo))
            if (lwrite) write (lfnsd, '(i6,f15.3)') k, lc(jepo)
            flglc(jepo) = 1
            if (istrue(flagall(k, isat), 'good') .and. istrue(flagall(k, iref), 'good')) then
              nsd = nsd + 1
              flglc(jepo) = 0
            end if
          end if
        end do
        if (nsd .gt. ndgr + 2) then
          if (lwrite) write (lfnsd, *) ' Satellite ', prn0(isat), prn0(iref)
          call check_for_jump(' LC ', lfnsd, kobs, ti(i0), lc, flglc, ndgr, niter, &
                              0.3d0, 0.15d0, a0, rmslc(isat), vlc, ilc, ierr, interval)
          if (ierr .eq. 0) then
            if (rmslc(isat) .gt. lclimit .and. lwrite) &
              write (lfnsd, '(a,f10.6)') ' LC big rms ', rmslc(isat)
            do jepo = 1, kobs
              if (ilc(jepo) .ne. 0) then
                k = i0 + jepo - 1
                if (ichecked(k, isat) .ne. 0) ichecked(k, isat) = 0
                if (ichecked(k, iref) .ne. 0) ichecked(k, iref) = 0
                if (ichecked(i0 + ilc(jepo) - 1, isat) .eq. 10) &
                  ichecked(i0 + ilc(jepo) - 1, isat) = 2
                if (ichecked(i0 + ilc(jepo) - 1, iref) .eq. 10) &
                  ichecked(i0 + ilc(jepo) - 1, iref) = 2
                if (dabs(vlc(jepo)) .gt. lclimit) then
                  if (lwrite) write (lfnsd, *) ' Found ', k, jepo, isat, iref, vlc(jepo)
                  if (flglgref(jepo) .ne. 1 .and. rmslg(iref) .le. 0.10) ref_ok = .true.
                  found = .false.
                  do i = 1, njump
                    if (jump(i, 1) .eq. k .and. (jump(i, 2) .eq. isat .or. &
                                                 jump(i, 3) .eq. isat) .and. (jump(i, 2) .eq. iref .or. &
                                                                              jump(i, 3) .eq. iref)) then
                      found = .true.
                      cycle
                    end if
                  end do
                  if (.not. found) then
                    njump = njump + 1
                    jump(njump, 1) = k
                    jump(njump, 2) = isat
                    jump(njump, 3) = iref
                    jump(njump, 4) = 2
                    if (ref_ok .and. istrue(flagall(k - 1, iref), 'ok')) jump(njump, 4) = 1
                  end if
                end if
              end if
            end do
          else
            if (lwrite) write (lfnsd, *) ' fit error, can not check this piece ', prn0(isat), i0, i1, ierr, rmslc(isat)
          end if
        end if
      end if
! next satellite
    end do
100 continue
! next epoch
    iepo = iepo + nstep
  end do

!!
!! count and remove unchecked data
  j = 0
  do iepo = 1, nepo
    do isat = 1, nsat
      if (prn0(isat)(1:1) .eq. 'C') then
        cycle
      end if
      if (ichecked(iepo, isat) .eq. 2) then
        if (lwrite) write (lfnsd, *) ' First Epoch ', iepo, prn0(isat)
        flagall(iepo, isat) = set_flag(flagall(iepo, isat), 'bigsd')
      else if (istrue(flagall(iepo, isat), 'ok') .and. ichecked(iepo, isat) .ne. 0) then
! shift flag to the next epoch
        if (.not. istrue(flagall(iepo, isat), 'good')) then
          do k = iepo + 1, nepo
            if (istrue(flagall(k, isat), 'ok')) then
              if (istrue(flagall(k, isat), 'good')) flagall(k, isat) = flagall(iepo, isat)
              exit
            end if
          end do
          j = j + 1
          flagall(iepo, isat) = set_flag(flagall(iepo, isat), 'lccheck')
          if (lwrite) write (lfnsd, *) ' Not checked ', iepo, prn0(isat)
        end if
      end if
    end do
  end do
  if (lwrite) write (lfnsd, *) ' Not checked', j

! statistics of flags and set in 'flagall'

  j = 0
  do i = 1, njump
    if (istrue(flagall(jump(i, 1), jump(i, 2)), 'ok')) then
      if (istrue(flagall(jump(i, 1), jump(i, 2)), 'good')) then
        j = j + 1
      end if
      flagall(jump(i, 1), jump(i, 2)) = set_flag(flagall(jump(i, 1), jump(i, 2)), 'bigsd')
    end if
    if (jump(i, 4) .eq. 2 .and. istrue(flagall(jump(i, 1), jump(i, 3)), 'ok')) then
      if (istrue(flagall(jump(i, 1), jump(i, 3)), 'good')) then
        j = j + 1
      end if
      flagall(jump(i, 1), jump(i, 3)) = set_flag(flagall(jump(i, 1), jump(i, 3)), 'bigsd')
    end if
  end do
  if (lwrite) write (lfnsd, *) ' TOTAL NEW FLAGs', j

  return
end subroutine
