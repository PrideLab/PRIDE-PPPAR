!
!! rdrnxoi2.f90
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
!! Contributor: Maorong Ge, Jianghui Geng, Jihang Lin
!! 
!!
!!
!! purpose  : read one epoch data from a RINEX o-file
!!
!! parameter: lfn -- file unit
!!            jd0, sod0 --- julian day and second of day of the requested epoch
!!                          if they are zero, take the epoch the file pointer
!!                          points at.
!!            wind  --- window for time matching. If the obsersing time from the
!!                      rinex-file is close to the requested time within the
!!                      window, we take the data. Be careful with this parameter
!!                      when you are working sampling rate larger than 1Hz.
!!            nprn0,prn0 -- number of satellite and satellite PRNs are chosen
!!                          If nprn is zero, take all observation of the matched
!!                          epoch.
!!            HD -- rinex header structure
!!            OB -- observation structure
!!            dcb -- P1C1 from CODE
!!            bias -- biases for pppar
!!            ierr -- error code, end of file or read fil error
!!
!
subroutine rdrnxoi2(lfn, jd0, sod0, dwnd, nprn0, prn0, HD, OB, bias, ierr)
  implicit none
  include '../header/const.h'
  include '../header/rnxobs.h'

  integer*4 ierr, lfn, jd0, nprn0
  character*3 prn0(1:*)
  real*8 sod0, dwnd
  type(rnxhdr) HD
  type(rnxobr) OB
!
!! local
  logical*1 prior_gp1, prior_gp2, prior_rp1, prior_rp2
  integer*4 ioerr, iy, im, id, ih, imi, nprn
  character*3 prn(MAXSAT)
  integer*4 iflag, i, j, i0, nline, ii, nobstype
  real*8 sec, ds, dt, obs(MAXTYP), bias(MAXSAT, MAXTYP)
  character*1 sysid(MAXSAT)
  character*80 line, cline, msg, name
  character*1024 string
  ! R
  integer*4 frequency_glo_nu,prn_int
  real*8 :: FREQ1_R(-50:50),FREQ2_R(-50:50)
  integer*4 biasW_index_G,biasC_index_G,biasL_index_G,&
            biasC_index_R,biasP_index_R,&
            biasX_index_E,biasI_index_C,biasL_index_J
!
!! function used
  integer*4 modified_julday, pointer_string

  call frequency_glonass(FREQ1_R,FREQ2_R)
  ierr = 0
  line = ' '
  prn = ''
  biasW_index_G = index(obs_prio_G, 'W')
  biasC_index_G = index(obs_prio_G, 'C')
  biasL_index_G = index(obs_prio_G, 'L')
  biasC_index_R = index(obs_prio_R, 'C')
  biasP_index_R = index(obs_prio_R, 'P')
  biasX_index_E = index(obs_prio_E, 'X')
  biasI_index_C = index(obs_prio_C, 'I')
  biasL_index_J = index(obs_prio_J, 'L')

10 continue  ! next record
  read (lfn, '(a)', end=200) line
  msg = ' '
!
!! in case of multiple headers in file
  if (index(line, "RINEX VERSION / TYPE") .ne. 0) then
    backspace lfn
    call rdrnxoh(lfn, HD, ioerr)
    goto 10
  endif
!
!! number of satellite
  read (line(30:32), '(i3)', iostat=ioerr) nprn
  if (ioerr .ne. 0) then
    msg = 'read satellite number error.'
  endif
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
      endif
    enddo
    goto 10
  endif
  if (nprn .gt. MAXSAT) then
    msg = 'satellite number > maxsat'
  endif
  if (len_trim(msg) .ne. 0) goto 100
!
!! initialization
  do i = 1, MAXSAT
    do j = 1, 6
      OB%obs(i, j) = 0.d0
    enddo
  enddo
!
!! format of the time tag line
  msg = 'read time & svn error'
  read (line, '(5i3,f11.7,2i3,12(a1,i2))', err=100) iy, im, id, ih, imi, sec, iflag, nprn
  if (nprn .gt. MAXSAT) then
    write (*, '(a,i3)') '***ERROR(rdrnxoi2): nprn > maxsat ', nprn
    call exit(1)
  endif
  read (line(68:80), '(f13.9)', iostat=ioerr) dt
  if (ioerr .ne. 0) dt = 0.d0
  read (line, '(32x,12(a3))', err=100) (prn(i), i=1, min(nprn, 12))
  do i=1,min(nprn,12)
    if(prn(i)(1:1) .eq. ' ' .and. prn(i)(3:3) .ne. ' ')prn(i)(1:1)='G'
    if(prn(i)(2:2) .eq. ' ' .and. prn(i)(3:3) .ne. ' ')prn(i)(2:2)='0'
    sysid(i)=prn(i)(1:1)
  end do
!
!! check time
  if (im .le. 0 .or. im .gt. 12 .or. id .le. 0 .or. id .gt. 31 .or. ih .lt. 0 .or. ih .ge. 24 &
      .or. imi .lt. 0 .or. imi .gt. 60 .or. sec .lt. 0.d0 .or. sec .gt. 60.d0) then
    msg = 'epoch time incorrect'
    goto 100
  endif
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
      enddo
      line = ' '
      goto 10
    endif
  endif
!
!! in case of mult lines for PRN number, continous lines must be read here.
  ii = 1
  do while (nprn .gt. 12*ii)
    read (lfn, '(32x,12(a3))', err=100) (prn(i), i=12*ii + 1, min(12*ii + 12, nprn))
    do i=12*ii+1,min(12*ii+12,nprn)
      if(prn(i)(1:1) .eq. ' ' .and. prn(i)(3:3) .ne. ' ')prn(i)(1:1)='G'
      if(prn(i)(2:2) .eq. ' ' .and. prn(i)(3:3) .ne. ' ')prn(i)(2:2)='0'
      sysid(i)=prn(i)(1:1)
    end do
    ii = ii + 1
  end do
!
!! read data. if more than 10 type 3 line should be merged to one
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
!
!! check if the sallite is requested
    if(sysid(i).ne.'G' .and. sysid(i).ne.'R' .and. sysid(i).ne.'E' .and. sysid(i).ne.'C' .and. sysid(i).ne.'J') cycle
    read(prn(i),'(1x,i2)') prn_int
    nobstype = min(HD%nobstyp, MAXTYP)
    read (string, '(30(f14.3,2x))', err=100) (obs(j), j=1, nobstype)
    i0 = 0
    if (nprn0 .gt. 0) then
      do j = 1, nprn0
        if (prn0(j) .eq. prn(i)) i0 = j
      enddo
    else
      i0 = i
    endif
!
!! fill in the obs. structure
    if (i0 .ne. 0) then
      prior_gp1 = .false.
      prior_gp2 = .false.
      prior_rp1 = .false.
      prior_rp2 = .false.
      do j = 1, HD%nobstyp
        if (sysid(i).eq.'G') then
          if (HD%obstyp(j) .eq. 'L1') then
            OB%obs(i0, 1) = obs(j) - bias(i0, biasW_index_G)*FREQ1_G/VLIGHT
            OB%typuse(i0, 1) = HD%obstyp(j)
          else if (HD%obstyp(j) .eq. 'L2') then
            OB%obs(i0, 2) = obs(j) - bias(i0, biasW_index_G+9)*FREQ2_G/VLIGHT
            OB%typuse(i0, 2) = HD%obstyp(j)
          else if (HD%obstyp(j) .eq. 'P1' .and. dabs(obs(j)) .gt. 1.d7) then
            OB%obs(i0, 3) = obs(j) - bias(i0, biasW_index_G+9*2)
            prior_gp1 = .true.
            OB%typuse(i0, 3) = HD%obstyp(j)
          else if (HD%obstyp(j) .eq. 'C1' .and. dabs(obs(j)) .gt. 1.d7) then
            if (.not. prior_gp1) then
              OB%obs(i0, 3) = obs(j) - bias(i0, biasC_index_G+9*2)
              OB%typuse(i0, 3) = HD%obstyp(j)
            endif
          else if (HD%obstyp(j) .eq. 'P2' .and. dabs(obs(j)) .gt. 1.d7) then
            OB%obs(i0, 4) = obs(j) - bias(i0, biasW_index_G+9*3)
            prior_gp2 = .true.
            OB%typuse(i0, 4) = HD%obstyp(j)
          else if (HD%obstyp(j) .eq. 'C2' .and. dabs(obs(j)) .gt. 1.d7) then
            if (.not. prior_gp2) then
              OB%obs(i0, 4) = obs(j) - bias(i0, biasL_index_G+9*3)
              OB%typuse(i0, 4) = HD%obstyp(j)
            endif
          endif
        elseif (sysid(i).eq.'R') then
          frequency_glo_nu=OB%glschn(prn_int)
          if (HD%obstyp(j) .eq. 'L1') then
            OB%obs(i0, 1) = obs(j)
            OB%typuse(i0, 1) = HD%obstyp(j)
          else if (HD%obstyp(j) .eq. 'L2') then
            OB%obs(i0, 2) = obs(j)
            OB%typuse(i0, 2) = HD%obstyp(j)
          else if (HD%obstyp(j) .eq. 'P1' .and. dabs(obs(j)) .gt. 1.d7) then
            OB%obs(i0, 3) = obs(j) - bias(i0, biasP_index_R+9*2)
            prior_rp1 = .true.
            OB%typuse(i0, 3) = HD%obstyp(j)
          else if (HD%obstyp(j) .eq. 'C1' .and. dabs(obs(j)) .gt. 1.d7) then
            if (.not. prior_rp1) then
              OB%obs(i0, 3) = obs(j) - bias(i0, biasC_index_R+9*2)
              OB%typuse(i0, 3) = HD%obstyp(j)
            endif
          else if (HD%obstyp(j) .eq. 'P2' .and. dabs(obs(j)) .gt. 1.d7) then
            OB%obs(i0, 4) = obs(j) - bias(i0, biasP_index_R+9*3)
            prior_rp2 = .true.
            OB%typuse(i0, 4) = HD%obstyp(j)
          else if (HD%obstyp(j) .eq. 'C2' .and. dabs(obs(j)) .gt. 1.d7) then
            if (.not. prior_rp2) then
              OB%obs(i0, 4) = obs(j) - bias(i0, biasC_index_R+9*3)
              OB%typuse(i0, 4) = HD%obstyp(j)
            endif
          endif
        elseif (sysid(i).eq.'E') then
          if (HD%obstyp(j) .eq. 'L1') then
            OB%obs(i0, 1) = obs(j) - bias(i0, biasX_index_E)*FREQ1_E/VLIGHT
            OB%typuse(i0, 1) = HD%obstyp(j)
          else if (HD%obstyp(j).eq.'L5') then
            OB%obs(i0, 2) = obs(j) - bias(i0, biasX_index_E+9)*FREQ2_E/VLIGHT
            OB%typuse(i0, 2) = HD%obstyp(j)
          else if (HD%obstyp(j) .eq. 'C1' .and. dabs(obs(j)) .gt. 1.d7) then
            OB%obs(i0, 3) = obs(j) - bias(i0, biasX_index_E+9*2)
            OB%typuse(i0, 3) = HD%obstyp(j)
          else if (HD%obstyp(j) .eq. 'C5' .and. dabs(obs(j)) .gt. 1.d7) then
            OB%obs(i0, 4) = obs(j) - bias(i0, biasX_index_E+9*3)
            OB%typuse(i0, 4) = HD%obstyp(j)
          endif
        elseif (sysid(i).eq.'C') then
          if (HD%obstyp(j) .eq. 'L2') then
            OB%obs(i0, 1) = obs(j) - bias(i0, biasI_index_C)*FREQ1_C/VLIGHT
            OB%typuse(i0, 1) = HD%obstyp(j)
          else if (HD%obstyp(j).eq.'L6') then ! BDS-2 B1I&B3I
          !else if (HD%obstyp(j).eq.'L7') then ! BDS-2 B1I&B2I
            OB%obs(i0, 2) = obs(j) - bias(i0, biasI_index_C+9)*FREQ2_C/VLIGHT
            OB%typuse(i0, 2) = HD%obstyp(j)
          else if (HD%obstyp(j) .eq. 'C2' .and. dabs(obs(j)) .gt. 1.d7) then
            OB%obs(i0, 3) = obs(j) - bias(i0, biasI_index_C+9*2)
            OB%typuse(i0, 3) = HD%obstyp(j)
          else if (HD%obstyp(j) .eq. 'C6' .and. dabs(obs(j)) .gt. 1.d7) then ! BDS-2 B1I&B3I
          !else if (HD%obstyp(j) .eq. 'C7' .and. dabs(obs(j)) .gt. 1.d7) then ! BDS-2 B1I&B2I
            OB%obs(i0, 4) = obs(j) - bias(i0, biasI_index_C+9*3)
            OB%typuse(i0, 4) = HD%obstyp(j)
          endif
        elseif (sysid(i).eq.'J') then
          if (HD%obstyp(j) .eq. 'L1') then
            OB%obs(i0, 1) = obs(j) - bias(i0, biasL_index_J)*FREQ1_J/VLIGHT
            OB%typuse(i0, 1) = HD%obstyp(j)
          else if (HD%obstyp(j).eq.'L2') then
            OB%obs(i0, 2) = obs(j) - bias(i0, biasL_index_J+9)*FREQ2_J/VLIGHT
            OB%typuse(i0, 2) = HD%obstyp(j)
          else if (HD%obstyp(j) .eq. 'C1' .and. dabs(obs(j)) .gt. 1.d7) then
            OB%obs(i0, 3) = obs(j) - bias(i0, biasL_index_J+9*2)
            OB%typuse(i0, 3) = HD%obstyp(j)
          else if (HD%obstyp(j) .eq. 'C2' .and. dabs(obs(j)) .gt. 1.d7) then
            OB%obs(i0, 4) = obs(j) - bias(i0, biasL_index_J+9*3)
            OB%typuse(i0, 4) = HD%obstyp(j)
          endif
        endif
      enddo
!! if one of the phases is zero, or the data is removed before
      if (any(OB%obs(i0, 1:4) .eq. 0.d0)) then
        OB%obs(i0, 1:4) = 0.d0
      endif
    endif
  enddo

  if (nprn0 .eq. 0) then
    OB%nprn = nprn
    do i = 1, nprn
      if(sysid(i).eq.'G' .or. sysid(i).eq.'R' .or. sysid(i).eq.'E' .or. sysid(i).eq.'C' .or. sysid(i).eq.'J') then
        OB%prn(i) = prn(i)
      else
        OB%prn(i) = ''
      endif
    enddo
  else
    OB%nprn = nprn0
    do i = 1, nprn0
      OB%prn(i) = prn0(i)
    enddo
  endif
  OB%dtrcv = dt

!
!! normal ending
  return
!
!! error
100 continue
  ierr = 1
  inquire (unit=lfn, name=name)
  write (*, '(a/a/a)') '***ERROR(rdrnxoi2): read file, '//trim(name), '   line :'//line(1:80), ' &
    msg  :'//trim(msg)
  call exit(1)
!
!! come here on end of file
200 continue
  ierr = 2
  inquire (unit=lfn, name=name)
  do i = 1, MAXSAT
    OB%obs(i, 1:4) = 0.d0
  enddo
  write (*, '(a/a)') '###WARNING(rdrnxoi2): end of file, ', trim(name)
  return
end
