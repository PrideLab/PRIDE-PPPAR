!
!! check_sd.f90
!!
!!    Copyright (C) 2021 by Wuhan University
!!
!!    This program is an open source software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License (version 3) as
!!    published by the Free Software Foundation.
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
                    nepo, nsat, jd0, nobs, flagall, ti, ts, obs, suitepo)
  implicit none
  include '../header/const.h'
  include '../header/brdeph.h'
  include 'data_flag.h'
!
  integer*4 suitepo
  integer*4 nepo, nsat, flagall(suitepo, MAXSAT), jd0, nobs(MAXSAT)
  real*8 ti(suitepo), ts(suitepo), obs(suitepo, MAXSAT, 6)
!
  integer*4 neph, lfnsd, ndgr, niter, mepo, nstep
  type(brdeph) EPHEM(MAXEPH)
  real*8 interval, lclimit, x(MAXEPO), y(MAXEPO), z(MAXEPO)
!
!! local
  logical*1 found, ref_ok, lwrite
  integer*4 i, j, k, i0, i1, istart, istop, ngood, iepo, kobs, ierr, jeph, jepo, nsd, isat, iref, njump
  integer*4 ichecked(suitepo, MAXSAT), jump(suitepo*10, 4), nobs_epoch(suitepo), nobsprn(MAXSAT), nobsflg(MAXSAT)
  integer*4 flglg(suitepo), flglc(suitepo), flglgref(suitepo), ilc(suitepo), ilg(suitepo)
  real*8 lc(suitepo), lg(suitepo), range_ref(suitepo), range(suitepo)
  real*8 vlc(suitepo), vlg(suitepo), vlgref(suitepo)
  real*8 rmslg(MAXSAT), rmslc(MAXSAT)
  real*8 dtsat, elev, dummy, a0
  ! G
  real*8 lambda1_G
  ! R
  integer*4 frequency_glo_nu,prn_int,glschn(30)
  character*3 prn_str
  real*8 :: FREQ1_R(-50:50),FREQ2_R(-50:50),lambda1_R(-50:50)
  ! E
  real*8 lambda1_E
  ! C
  real*8 lambda1_C
  ! J
  real*8 lambda1_J
  
  character*3 prn0(MAXSAT)
!
!! function called
  logical*1 istrue
  integer*4 set_flag,rdglschn

  call frequency_glonass(FREQ1_R,FREQ2_R)
  call prn_matbld(prn0)
!
!! initialization
  lwrite = .true.
  if (lfnsd .eq. 0) lfnsd = 6
  if (lfnsd .lt. 0) lwrite = .false.
  ! R
  do i=1,30
    write(prn_str,'(a1,i2.2)') 'R',i
    glschn(i)=rdglschn(prn_str,jd0,ti(1))
  enddo
  ! G
  lambda1_G=VLIGHT/FREQ1_G
  ! E
  lambda1_E=VLIGHT/FREQ1_E
  ! C
  lambda1_C=VLIGHT/FREQ1_C
  ! J
  lambda1_J=VLIGHT/FREQ1_J
  ierr = 0
!! find out first and last epoch with enough data
!! nobs_epoch -- # of obs. per epoch, the number of available satellite
  istart = 0     ! start epoch
  istop = 0     ! end epoch
  ngood = 0     ! Effective number of epoch
  do iepo = 1, nepo
    nobs_epoch(iepo) = 0
    do isat = 1, nsat
      if(prn0(isat)(1:1) .eq. 'C')then
        cycle
      endif
      ichecked(iepo, isat) = 10
      if (istrue(flagall(iepo, isat), 'ok')) then
        nobs_epoch(iepo) = nobs_epoch(iepo) + 1    ! the number of available satellite
      endif
    enddo
    if (nobs_epoch(iepo) .ge. 2) then
      if (istart .eq. 0) istart = iepo
      istop = iepo
      ngood = ngood + 1
    endif
  enddo
  if (istop - istart + 1 .le. ndgr + 2 .or. ngood .le. ndgr + 2) then
! remove all
    return
  endif

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
    enddo

! ndgr = 2; default value
    if (ngood .le. ndgr + 2) goto 100

! kobs ------ the number of epoch in check arc
    kobs = i1 - i0 + 1

! poly fit for each satellite in epoch section to choose reference satellite
!-nobsprn : # of observations for each satellite
!-nobsflg : # of flags (amb.) for each satellite
    do isat = 1, nsat
      if(prn0(isat)(1:1) .eq. 'C')then
        cycle
      endif
      nobsprn(isat) = 0
      nobsflg(isat) = 0
      do jepo = 1, kobs
        k = i0 + jepo - 1
        if (istrue(flagall(k, isat), 'ok')) then
          nobsprn(isat) = nobsprn(isat) + 1     ! the number of available epochs for each satellite
        endif
      enddo
      if (nobsprn(isat) .gt. ndgr + 2) then
        do jepo = 1, kobs
          k = i0 + jepo - 1
          flglg(jepo) = 10
          if (istrue(flagall(k, isat), 'ok')) then
            flglg(jepo) = 1
            lg(jepo) = obs(k, isat, 1)          ! geometry-free obs value
          endif
          if (istrue(flagall(k, isat), 'good')) then
            flglg(jepo) = 0
          endif
        enddo
        call check_for_jump(' LG ', lfnsd, kobs, ti(i0), lg, flglg, ndgr, niter, &
                            2.0d0, 0.40d0, a0, rmslg(isat), vlg, ilg, ierr, interval)
        if (ierr .ne. 0) then
        endif
        do jepo = 1, kobs
          k = i0 + jepo - 1
          if (flglg(jepo) .eq. 1) nobsflg(isat) = nobsflg(isat) + 1
        enddo
        if (ierr .ne. 0) rmslg(isat) = 1.d3
      else
        nobsflg(isat) = ngood
        rmslg(isat) = 1.d3
      endif
    enddo
!!
!! select satisfactory referece satellte
    iref = 1
    do isat = 1, nsat
      if(prn0(isat)(1:1) .ne. 'G')then !add for referece satellte
        cycle
      endif
      if (nobsprn(isat) .gt. nobsprn(iref)) then
        iref = isat
      else if (nobsprn(isat) .eq. nobsprn(iref)) then
        if (nobsflg(isat) .lt. nobsflg(iref)) then
          iref = isat
        else if (nobsflg(isat) .eq. nobsflg(iref)) then
          if (rmslg(isat) .le. rmslg(iref)) iref = isat
        endif
      endif
    enddo
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
      endif
    enddo
    call check_for_jump(' LGR ', lfnsd, kobs, ti(i0), lg, flglg, ndgr, niter, &
                        2.0d0, 0.40d0, a0, rmslg(iref), vlg, ilg, ierr, interval)
    if (ierr .ne. 0) then
    endif
    do jepo = 1, kobs
      k = i0 + jepo - 1
      flglgref(jepo) = flglg(jepo)
      vlgref(jepo) = 0.d0
      if (ilg(jepo) .ne. 0) then
        vlgref(jepo) = vlgref(ilg(jepo)) + vlg(jepo)
      endif
    enddo
!!
!! check other satellites, form single diff LC obs
    do isat = 1, nsat
      if(prn0(isat)(1:1) .eq. 'C')then
        cycle
      endif
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
            ! LC single-difference = troposphere delay or ...
            if(prn0(isat)(1:1).eq.'G')then 
              lc(jepo) = (obs(k, isat, 3)-range(jepo)/lambda1_G)-(obs(k, iref, 3)-range_ref(jepo)/lambda1_G)
            elseif(prn0(isat)(1:1).eq.'R')then
              read(prn0(isat),'(1x,i2)') prn_int
              frequency_glo_nu=glschn(prn_int)
              lambda1_R(frequency_glo_nu)=VLIGHT/FREQ1_R(frequency_glo_nu)
              lc(jepo) = (obs(k, isat, 3)-range(jepo)/lambda1_R(frequency_glo_nu))-(obs(k, iref, 3)-range_ref(jepo)/lambda1_G)
            elseif(prn0(isat)(1:1).eq.'E')then
              lc(jepo) = (obs(k, isat, 3)-range(jepo)/lambda1_E)-(obs(k, iref, 3)-range_ref(jepo)/lambda1_G)
            elseif(prn0(isat)(1:1).eq.'C')then
              lc(jepo) = (obs(k, isat, 3)-range(jepo)/lambda1_C)-(obs(k, iref, 3)-range_ref(jepo)/lambda1_G)
            elseif(prn0(isat)(1:1).eq.'J')then
              lc(jepo) = (obs(k, isat, 3)-range(jepo)/lambda1_J)-(obs(k, iref, 3)-range_ref(jepo)/lambda1_G)
            end if
            if (lwrite) write (lfnsd, '(i6,f15.3)') k, lc(jepo)
            flglc(jepo) = 1
            if (istrue(flagall(k, isat), 'good') .and. istrue(flagall(k, iref), 'good')) then
              nsd = nsd + 1
              flglc(jepo) = 0
            endif
          endif
        enddo
        if (nsd .gt. ndgr + 2) then
          if (lwrite) write (lfnsd, *) ' Satellite ', prn0(isat),prn0(iref)
          call check_for_jump(' LC ', lfnsd, kobs, ti(i0), lc, flglc, ndgr, niter, &
                              0.2d0, 0.12d0, a0, rmslc(isat), vlc, ilc, ierr, interval)
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
                    endif
                  enddo
                  if (.not. found) then
                    njump = njump + 1
                    jump(njump, 1) = k
                    jump(njump, 2) = isat
                    jump(njump, 3) = iref
                    jump(njump, 4) = 2
                    if (ref_ok .and. istrue(flagall(k - 1, iref), 'ok')) jump(njump, 4) = 1
                  endif
                endif
              endif
            enddo
          else
            if (lwrite) write (lfnsd, *) ' fit error, can not check this piece ', prn0(isat), i0, i1, ierr, rmslc(isat)
          endif
        endif
      endif
! next satellite
    enddo
100 continue
! next epoch
    iepo = iepo + nstep
  enddo

!!
!! count and remove unchecked data
  j = 0
  do iepo = 1, nepo
    do isat = 1, nsat
      if(prn0(isat)(1:1) .eq. 'C')then
        cycle
      endif
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
            endif
          enddo
          j = j + 1
          flagall(iepo, isat) = set_flag(flagall(iepo, isat), 'lccheck')
          if (lwrite) write (lfnsd, *) ' Not checked ', iepo, prn0(isat)
        endif
      endif
    enddo
  enddo
  if (lwrite) write (lfnsd, *) ' Not checked', j

! statistics of flags and set in 'flagall'
  j = 0
  do i = 1, njump
    if (istrue(flagall(jump(i, 1), jump(i, 2)), 'ok')) then
      if (istrue(flagall(jump(i, 1), jump(i, 2)), 'good')) then
        j = j + 1
      endif
      flagall(jump(i, 1), jump(i, 2)) = set_flag(flagall(jump(i, 1), jump(i, 2)), 'bigsd')
    endif
    if (jump(i, 4) .eq. 2 .and. istrue(flagall(jump(i, 1), jump(i, 3)), 'ok')) then
      if (istrue(flagall(jump(i, 1), jump(i, 3)), 'good')) then
        j = j + 1
      endif
      flagall(jump(i, 1), jump(i, 3)) = set_flag(flagall(jump(i, 1), jump(i, 3)), 'bigsd')
    endif
  enddo
  if (lwrite) write (lfnsd, *) ' TOTAL NEW FLAGs', j

  return
end
