!
!! rdrnxoi2.f90
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
!! Contributor: Maorong Ge, Jianghui Geng, Jihang Lin, Jing Zeng
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
subroutine rdrnxoi2(lfn, jd0, sod0, dwnd, nprn0, prn0, HD, OB, bias, nbias_used, ierr)
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
  integer*4     iflag, i, j, k, i0, isys, nline
  integer*4     ii, nobstype
  integer*1     lli(MAXTYP), lli_thre
  character*3   prn(MAXSAT)
  real*8        sec, ds, dt, obs(MAXTYP), vobs
  character*1   sysid(MAXSAT)
  character*80  line, cline, msg, name
  character*512 string
!! RINEX-2 Signal Priority
  logical*1     used_codeP(2)
  integer*4     prio_index
  integer*4     imes, ityp
  integer*4     obs_prio_index(6)       !! RINEX-3 signal code of each measurement type
  character*2   obs_styp_gnss(6)        !! measurement type (L/P) + frequency number(1-9)
  real*8        lambda(4)
! function used
  integer*4     modified_julday
  integer*4     pointer_string

!
!! initialize record
  ierr = 0
  line = ' '
  prn = ' '
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
    call rdrnxoh(lfn, HD, ioerr)
    goto 10
  end if
!
!! for some special event
  if (index(line, "COMMENT") .ne. 0) goto 10
!
!! number of satellite
  read (line(30:32), '(i3)', iostat=ioerr) nprn
  if (ioerr .ne. 0) then
    msg = 'read satellite number error.'
  end if
!
!! check the RINEX 2 event flag
  read (line(27:29), '(i3)', iostat=ioerr) iflag
  if (ioerr .ne. 0) then
    msg = 'read event flag error.'
    goto 100
  else if (iflag .gt. 1) then
    do i = 1, nprn
      read (lfn, '(a80)', iostat=ioerr, end=200) line
      if (line(61:80) .eq. 'ANTENNA: DELTA H/E/N') then
        msg = 'read internal antenna information error'
        read (line, '(3f14.4)', err=100) HD%h, HD%e, HD%n
      end if
    end do
    goto 10
  end if
  if (nprn .gt. MAXSAT) then
    msg = 'satellite number > maxsat'
  end if
  if (len_trim(msg) .ne. 0) goto 100
!
!! initialization
  OB%obs = 0.d0
  OB%lli = 0
!
!! format of the time tag line
  msg = 'read time & svn error'
  read (line, '(5i3,f11.7,2i3,12(a1,i2))', err=100) iy, im, id, ih, imi, sec, iflag, nprn
  if (nprn .gt. MAXSAT) then
    write (*, '(a,i3)') '***ERROR(rdrnxoi2): nprn > maxsat ', nprn
    call exit(1)
  end if
  read (line(68:80), '(f13.9)', iostat=ioerr) dt
  if (ioerr .ne. 0) dt = 0.d0
  read (line, '(32x,12(a3))', err=100) (prn(i), i=1, min(nprn, 12))
  do i = 1, min(nprn, 12)
    if (prn(i)(1:1) .eq. ' ' .and. prn(i)(3:3) .ne. ' ') prn(i)(1:1) = 'G'
    if (prn(i)(2:2) .eq. ' ' .and. prn(i)(3:3) .ne. ' ') prn(i)(2:2) = '0'
    sysid(i) = prn(i)(1:1)
  end do
!
!! check time
  if (im .le. 0 .or. im .gt. 12 .or. id .le. 0 .or. id .gt. 31 .or. ih .lt. 0 .or. ih .ge. 24 &
      .or. imi .lt. 0 .or. imi .gt. 60 .or. sec .lt. 0.d0 .or. sec .gt. 60.d0) then
    msg = 'epoch time incorrect'
    goto 100
  end if
  call yr2year(iy)
!
!! check on time tags. do not change the requested time if there is no data
  ds = 0.d0
  OB%jd = modified_julday(id, im, iy)
  OB%tsec = ih*3600.d0 + imi*60.d0 + sec
  if (jd0 .ne. 0) then
    ds = (jd0 - OB%jd)*86400.d0 + (sod0 - OB%tsec)
    if (ds .lt. -dwnd) then
      OB%jd = jd0
      OB%tsec = sod0
      backspace lfn
      OB%nprn = 0
      return
    else if (ds .gt. dwnd) then
      i = 0
      if (nprn .gt. 12) i = i + 1
      if (nprn .gt. 24) i = i + 1
      if (nprn .gt. 36) i = i + 1
      i = i + ((HD%nobstyp - 1)/5 + 1)*nprn
      do j = 1, i
        read (lfn, '(a)') line
      end do
      line = ' '
      goto 10
    end if
  end if
!
!! in case of mult lines for PRN number, continous lines must be read here.
  ii = 1
  do while (nprn .gt. 12*ii)
    read (lfn, '(32x,12(a3))', err=100) (prn(i), i=12*ii + 1, min(12*ii + 12, nprn))
    do i = 12*ii + 1, min(12*ii + 12, nprn)
      if ((prn(i)(1:1) .eq. ' ') .and. (prn(i)(3:3) .ne. ' ')) prn(i)(1:1) = 'G'
      if ((prn(i)(2:2) .eq. ' ') .and. (prn(i)(3:3) .ne. ' ')) prn(i)(2:2) = '0'
      sysid(i) = prn(i)(1:1)
    end do
    ii = ii + 1
  end do
!
!! read data
!! if more than 10 type 3 line should be merged to one
  do i = 1, nprn
    read (lfn, '(a80)', err=100, end=200) string
    cline = ' '
    ii = 1
    do while (HD%nobstyp .gt. 5*ii)
      read (lfn, '(a80)', err=100, end=200) cline
      string = string(1:80*ii)//cline
      cline = ' '
      ii = ii + 1
    end do
    isys = index(GNSS_PRIO, sysid(i))
    if (isys .eq. 0) cycle
!
!! get the index and the priority of each GNSS and signal type
    lambda = 1.d0
    if (index('GR', sysid(i)) .ne. 0) then
      if ('G' .eq. sysid(i)) then
        obs_prio_index(:) = index(obs_prio_G, 'W')
        obs_prio_index(5) = index(obs_prio_G, 'C')
        obs_prio_index(6) = index(obs_prio_G, 'L')
      else
        obs_prio_index(:) = index(obs_prio_G, 'W')
        obs_prio_index(5) = index(obs_prio_G, 'C')
        obs_prio_index(6) = index(obs_prio_G, 'L')
      end if
      do j = 1, 2
        k = idxfrq(isys, j)
        write (obs_styp_gnss(j + 0), '(a1,i1)') 'L', k
        write (obs_styp_gnss(j + 2), '(a1,i1)') 'P', k
        write (obs_styp_gnss(j + 4), '(a1,i1)') 'C', k
        lambda(j) = VLIGHT/FREQ_SYS(k, isys)
      end do
    elseif (index('ECJ', sysid(i)) .ne. 0) then
      if ('E' .eq. sysid(i)) then
        obs_prio_index(:) = index(obs_prio_E, 'X')
        obs_prio_index(5) = 0
        obs_prio_index(6) = 0
      elseif ('C' .eq. sysid(i)) then
        obs_prio_index(:) = index(obs_prio_C, 'I')
        obs_prio_index(5) = 0
        obs_prio_index(6) = 0
      elseif ('J' .eq. sysid(i)) then
        obs_prio_index(:) = index(obs_prio_J, 'L')
        obs_prio_index(5) = 0
        obs_prio_index(6) = 0
      end if
      do j = 1, 2
        k = idxfrq(isys, j)
        write (obs_styp_gnss(j + 0), '(a1,i1)') 'L', k
        write (obs_styp_gnss(j + 2), '(a1,i1)') 'C', k
        write (obs_styp_gnss(j + 4), '(a2)') 'XX'
        lambda(j) = VLIGHT/FREQ_SYS(k, isys)
      end do
    end if
!
!! check if the sallite is requested
    nobstype = min(HD%nobstyp, MAXTYP)
    read (string, '(30(f14.3,i1,1x))', err=100) (obs(j), lli(j), j=1, nobstype)
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
      used_codeP = .false.
      do j = 1, HD%nobstyp
        !
        !! check the index of prior signal & measurement type (L/C)
        do imes = 1, 6
          if (index(HD%obstyp(j), obs_styp_gnss(imes)) .ne. 0) then
            prio_index = obs_prio_index(imes)
            if (imes .le. 2) exit
            if (dabs(obs(j)) .gt. 1.d7) then
              if (imes .le. 4) then
                used_codeP(imes - 2) = .true.
              else
                if (used_codeP(imes - 4)) goto 50
              end if
              exit
            end if
          end if
          if (imes .eq. 6) goto 50
        end do
        if (imes .gt. 4) imes = imes - 2
        !
        !! GLONASS
        if (sysid(i) .eq. 'R') then
          if (HD%obstyp(j)(1:1) .ne. 'L') goto 20
          !
          !! no absolute bias for GLONASS phase observation
          OB%obs(i0, imes) = obs(j)
          if (imes .le. 2) OB%lli(i0, imes) = iand(lli(j), lli_thre)
          OB%typuse(i0, imes) = HD%obstyp(j)
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
          if (abs(bias(i0, ityp)%val(iepo) - 1.d9) .gt. 1.d-3) then
            vobs = bias(i0, ityp)%val(iepo)/lambda(imes)
            if (bias(i0, ityp)%length .gt. 1) then
              if (abs(bias(i0, ityp)%grd(iepo) - 1.d9) .gt. 1.d-3) then
                vobs = vobs + bias(i0, ityp)%grd(iepo)/lambda(imes) * ds
              end if
            end if
            OB%obs(i0, imes) = obs(j) - vobs
            if (imes .le. 2) OB%lli(i0, imes) = iand(lli(j), lli_thre)
            nbias_used(i0) = nbias_used(i0) + 1
            OB%typuse(i0, imes) = HD%obstyp(j)
            goto 50
          end if
        end if
        OB%obs(i0, imes) = obs(j)
        OB%typuse(i0, imes) = HD%obstyp(j)
        if (imes .le. 2) OB%lli(i0, imes) = iand(lli(j), lli_thre)
50      continue
      end do
      OB%itypuse(i0, imes) = -1
      if (abs(OB%tsec) .le. MAXWND .or. abs(OB%tsec - 864.d2) .le. MAXWND) then
        j = index(obs_prio_sys(isys), OB%typuse(i0, imes)(3:3))
        if (j .gt. 0) OB%itypuse(i0, imes) = j + (imes - 1) * 9
      end if
!
!! if one of the phases is zero, or the data is removed before
      if (any(OB%obs(i0, 1:4) .eq. 0.d0)) then
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
  write (*, '(a/a/a)') '***ERROR(rdrnxoi2): read file, '//trim(name), &
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
  write (*, '(a/a)') '###WARNING(rdrnxoi2): end of file, ', trim(name)
  return
end subroutine
