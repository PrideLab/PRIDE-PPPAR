!
!! rdrnxoi3.f90
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
!! Contributor: Jianghui Geng, Songfeng Yang, Jihang Lin, Jing Zeng
!!
!!
!!
!! purpose  : read one epoch data from a RINEX o-file
!!
!! parameter: lfn -- file unit
!!            jd0, sod0 --- julian day and second of day of the requested epoch
!!                          if they are zero, take the epoch the file pointer
!!                          points at.
!!            dwnd  --- window for time matching. If the obsersing time from the
!!                      rinex-file is close to the requested time within the
!!                      window, we take the data. Be careful with this parameter
!!                      when you are working sampling rate larger than 1Hz.
!!            nprn0,prn0 -- number of satellite and satellite PRNs are chosen
!!                          If nprn is zero, take all observation of the matched
!!                          epoch.
!!            HD -- rinex header structure
!!            OB -- observation structure
!!            bias, nbias_used -- observable-specific biases and the number of 
!!                                observations that are applied these biases
!!            ierr -- error code, end of file or read fil error
!!
!
subroutine rdrnxoi3(lfn, jd0, sod0, dwnd, nprn0, prn0, HD, OB, bias, nbias_used, ierr)
  implicit none
  include '../header/const.h'
  include '../header/absbia.h'
  include '../header/rnxobs.h'

! common
  integer*4     idxfrq(MAXSYS, 2)
  common        idxfrq
! parameter
  integer*4     lfn, jd0, nprn0, ierr
  character*3   prn0(1:*)
  real*8        sod0, dwnd
  type(rnxhdr)  HD
  type(rnxobr)  OB
  type(absbia)  bias(MAXSAT, MAXTYP)
  integer*4     nbias_used(MAXSAT)
! local
  integer*4     ioerr, iy, im, id, ih, imi, iepo, nprn
  integer*4     iflag, i, j, k, i0, isys, nline, nobstyp
  integer*1     lli(MAXTYP), lli_thre
  real*8        sec, ds, dt
  real*8        obs(MAXTYP), vobs
  character*3   prn(MAXSAT)
  character*1   sysid(MAXSAT)
  character*80  line, msg, name
!! RINEX-3 signal priority
  integer*4     obs_prio_index(4)
  integer*4     phs_prio_index(4)
  integer*4     obs_used_index(4)
  integer*4     best_index              !! GPS: W, GLONASS: P, Galileo: X, BDS: I, QZSS: L
  integer*4     prio_index
  integer*4     ityp, imes
  character*512 string
  character*3   hobstyp(MAXTYP)         !! HD%obstyp3
  character*2   mes_type_sys(4)         !! Measurement type (L/C) + Frequency number(1-9)
  character*3   ctyp
  real*8        lambda(4)
! function used
  integer*4     modified_julday
  integer*4     pointer_string
  integer*4     prn_int

!
!! initialize record
  ierr = 0
  line = ' '
  prn  = ' '
  nbias_used = -1
  lli_thre = 3
!
!! read the next record
10 continue
  read (lfn, '(a)', end=200) line
  msg = ' '
!
!! in case of multiple headers in file
  if (index(line, "RINEX VERSION / TYPE") .ne. 0) then
    backspace lfn
    call rdrnxoh(lfn, HD, ierr)
    goto 10
  end if
!
!! the start line
  if (line(1:1) .ne. '>') then
    goto 100
  end if
!
!! number of satellite
  read (line(33:35), '(i3)', iostat=ioerr) nprn
  if (ioerr .ne. 0) then
    msg = 'read satellite number error.'
  end if
  if (len_trim(msg) .ne. 0) goto 100
!
!! check the RINEX 3 event flag
  read (line(30:32), '(i3)', iostat=ioerr) iflag
  if (ioerr .ne. 0) then
    msg = 'read event flag error.'
    goto 100
  else if (iflag .gt. 1) then
    msg = 'read internal antenna information error'
    do i = 1, nprn
      read (lfn, '(a80)', iostat=ioerr, end=200) line
      if (line(61:80) .eq. 'ANTENNA: DELTA H/E/N') then
        read (line, '(3f14.4)', err=100) HD%h, HD%e, HD%n
      end if
    end do
    goto 10
  end if
!
!! initialize observations
  OB%obs = 0.d0
  OB%lli = 0.d0
!
!! format of the time tag line
  msg = 'read time & svn error'
  read (line, '(1x,i5,4i3,f11.7,i3,i3)', err=100) iy, im, id, ih, imi, sec, iflag, nprn
  read (line(42:56), '(f15.12)', iostat=ioerr) dt
  if (ioerr .ne. 0) dt = 0.d0
  if (nprn .gt. MAXSAT) then
    write (*, '(a,i3)') '***ERROR(rdrnxoi3): nprn > maxsat ', nprn
    call exit(1)
  end if
!
!! Check epoch time
  if (im .le. 0 .or. im .gt. 12 .or. id .le. 0 .or. id .gt. 31 .or. ih .lt. 0 .or. ih .ge. 24 &
      .or. imi .lt. 0 .or. imi .gt. 60 .or. sec .lt. 0.d0 .or. sec .gt. 60.d0) then
    msg = 'epoch time incorrect'
    goto 100
  end if
!
!! check on time tags. do not change the requested time if there is no data
  ds = 0.d0
  OB%jd = modified_julday(id, im, iy)
  OB%tsec = ih*3600.d0 + imi*60.d0 + sec
  if (jd0 .ne. 0) then
    ds = (jd0 - OB%jd)*864.d2 + (sod0 - OB%tsec)
    if (ds .lt. -dwnd) then
      OB%jd = jd0
      OB%tsec = sod0
      backspace lfn
      OB%nprn = 0
      return
    else if (ds .gt. dwnd) then
      i = nprn
      do j = 1, i
        read (lfn, '(a)') line
      end do
      line = ' '
      goto 10
    end if
  end if
!
!! read data
!! if more than 10 type 3 line should be merged to one
  do i = 1, nprn
    read (lfn, '(a)', err=100, end=200) string
    sysid(i) = string(1:1)
    read (string, '(1x,i2)') prn_int
    isys = index(GNSS_PRIO, sysid(i))
!
!! get the index and the priority of each GNSS and signal type
    lambda = 1.d0
    nobstyp = HD%nobstyp3_sys(isys)
    hobstyp = HD%cobstyp3_sys(1:MAXTYP, isys)
    do j = 1, 2
      k = idxfrq(isys, j)
      lambda(j) = VLIGHT/FREQ_SYS(k, isys)
      if (sysid(i) .eq. 'C' .and. HD%ver .eq. 302 .and. k .eq. 2) k = 1
      write (mes_type_sys(j + 0), '(a1,i1)') 'L', k
      write (mes_type_sys(j + 2), '(a1,i1)') 'C', k
    end do
    if (sysid(i) .eq. 'G') then
      best_index = index(OBS_PRIO_G, 'W')
    elseif (sysid(i) .eq. 'R') then
      best_index = index(OBS_PRIO_R, 'P')
    elseif (sysid(i) .eq. 'E') then
      best_index = index(OBS_PRIO_E, 'X')
    elseif (sysid(i) .eq. 'C') then
      best_index = index(OBS_PRIO_C, 'I')
    elseif (sysid(i) .eq. 'J') then
      best_index = index(OBS_PRIO_J, 'L')
    else
      cycle
    end if
    !
    !! check if the sallite is requested
    read (string, '(a3,50(f14.3,i1,1x))', err=100) prn(i), (obs(j), lli(j), j=1, min(nobstyp, MAXTYP))
    if ((prn(i)(1:1) .ne. ' ') .and. (prn(i)(2:2) .eq. ' ') .and. (prn(i)(3:3) .ne. ' ')) prn(i)(2:2) = '0'
    i0 = 0
    if (nprn0 .gt. 0) then
      do j = 1, nprn0
        if (prn0(j) .eq. prn(i)) i0 = j
      end do
    else
      i0 = i
    end if
!
!! fill in the obs. structure
    if (i0 .ne. 0) then
      nbias_used(i0) = 0
      obs_used_index = 0
      obs_prio_index = 0
      phs_prio_index = 0
      do j = 1, nobstyp
        prio_index = 0
        ctyp = hobstyp(j)(1:3)
        if (dabs(obs(j)) .lt. MAXWND .or. dabs(obs(j)) .gt. 1.d12) cycle
        !
        !! check measurement type (L/C)
        if ((ctyp(3:3) .eq. ' ') .and. &
           ((ctyp(1:1) .eq. 'L') .or. ((ctyp(1:1) .eq. 'C') .and. (sysid(i) .eq. 'E')))) then
          prio_index = best_index
        else
          prio_index = index(obs_prio_sys(isys), ctyp(3:3))
        end if
        if (prio_index .eq. 0) cycle
        do imes = 1, 4
          if (index(ctyp, mes_type_sys(imes)) .ne. 0) exit
          if (imes .eq. 4) goto 50
        end do
        !
        !! GLONASS
        if (sysid(i) .eq. 'R') then
          if (ctyp(1:1) .ne. 'L') goto 20
          !
          !! no absolute bias for GLONASS phase observation
          if (obs_prio_index(imes) .lt. prio_index) then
            obs_prio_index(imes) = prio_index
            OB%obs(i0, imes) = obs(j)
            if (imes .le. 2) OB%lli(i0, imes) = iand(lli(j), lli_thre)
            OB%typuse(i0, imes) = ctyp
          end if
          cycle
        end if
20      continue
        !
        !! calculate epoch index
        ityp = prio_index + (imes - 1) * 9
        if (bias(i0, ityp)%length .gt. 0) then
          iepo = int(OB%tsec/bias(i0, ityp)%period) + 1
          if (iepo .lt. 1) iepo = 1
          if (iepo .gt. bias(i0, ityp)%length) iepo = bias(i0, ityp)%length
          ds = 0.d0
          if (bias(i0, ityp)%length .gt. 1) then
            ds = OB%tsec - (iepo-0.5)*bias(i0, ityp)%period
          end if
        end if
        !
        !! correct raw observations
        if ((.not. allocated(bias(i0, ityp)%val)) .or. &
            (abs(bias(i0, ityp)%val(iepo) - 1.d9) .lt. 1.d-3)) then
          if (phs_prio_index(imes) .lt. prio_index) then
            obs_used_index(imes) = j
            phs_prio_index(imes) = prio_index
          end if
        else
          if (obs_prio_index(imes) .lt. prio_index) then
            obs_prio_index(imes) = prio_index
            vobs = bias(i0, ityp)%val(iepo)/lambda(imes)
            if (bias(i0, ityp)%length .gt. 1) then
              if (abs(bias(i0, ityp)%grd(iepo) - 1.d9) .gt. 1.d-3) then
                vobs = vobs + bias(i0, ityp)%grd(iepo)/lambda(imes) * ds
              end if
            end if
            OB%obs(i0, imes) = obs(j) - vobs
            if (imes .le. 2) OB%lli(i0, imes) = iand(lli(j), lli_thre)
            OB%typuse(i0, imes) = ctyp
            nbias_used(i0) = nbias_used(i0) + 1
          end if
        end if
50      continue
      end do
!
!! if one of the phases is zero, or the data is removed before
      do imes = 1, 4
        if (abs(OB%obs(i0, imes)) .lt. 1.d-9 .or. &
            abs(OB%obs(i0, imes)) .gt. 1.d12) then
          k = obs_used_index(imes)
          OB%obs(i0, imes) = obs(k)
          if (imes .le. 2) OB%lli(i0, imes) = iand(lli(k), lli_thre)
          if ('R' .eq. sysid(i) .and. imes .le. 2) cycle
          OB%typuse(i0, imes) = hobstyp(k)
          if (k .eq. 0) cycle
        end if
      end do
      if (any(abs(OB%obs(i0, 1:4)) .lt. 1.d-9) .or. &
          any(abs(OB%obs(i0, 1:4)) .gt. 1.d12)) then
        OB%obs(i0, 1:4) = 0.d0
        OB%lli(i0, 1:2) = 0
        nbias_used(i0) = -1
      end if
    end if
  end do
!
!! check PRN
  if (nprn0 .eq. 0) then
    OB%nprn = nprn
    do i = 1, nprn
      if (index(GNSS_PRIO, sysid(i)) .ne. 0) then
        OB%prn(i) = prn(i)
      else
        OB%prn(i) = ''
      end if
    end do
  else
    OB%nprn = nprn0
    OB%prn(1:nprn0) = prn0(1:nprn0)
  end if
  OB%dtrcv = dt
!
!! normal ending
  return
!
!! error
100 continue
  ierr = 1
  inquire (unit=lfn, name=name)
  write (*, '(a/a/a)') '***ERROR(rdrnxoi3): read file, '//trim(name), &
                       '   line :'//line(:), '   msg  :'//trim(msg)
  call exit(1)
!
!! come here on end of file
200 continue
  ierr = 2
  inquire (unit=lfn, name=name)
  do i = 1, MAXSAT
    OB%obs(i, 1:4) = 0.d0
    OB%lli(i, 1:2) = 0
  end do
  write (*, '(a/a)') '###WARNING(rdrnxoi3): end of file, ', trim(name)
  return
end subroutine
