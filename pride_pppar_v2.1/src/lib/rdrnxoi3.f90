!
!! rdrnxoi3.f90
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
!! Contributor: Jianghui Geng, Songfeng Yang
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
!!            ierr -- error code, end of file or read fil error
!!
!
subroutine rdrnxoi3(lfn, jd0, sod0, dwnd, nprn0, prn0, HD, OB, bias, ierr)
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
  integer*4 ioerr, iy, im, id, ih, imi, nprn, l1, l2, p1, p2
  character*3 prn(MAXSAT)
  integer*4 iflag, i, j, i0, nline, ii, nobstype, iii
  real*8 sec, ds, dt, obs(MAXTYP), bias(MAXSAT, 36)
  character*1 sysid(MAXSAT)
  character*80 line, cline, msg, name
  ! R
  integer*4 frequency_glo_nu
  real*8 :: freq1_R(-50:50),freq2_R(-50:50)
!
!! RINEX-3 Signal Priority
  integer*4 obs_prio_index(4),phs_prio_index(4)
  integer*4 prio_index
  integer*4 biasW_index_G,biasX_index_E,biasI_index_C,biasL_index_J
  character*1024 string
  integer*4 nobstyp_tmp
!
!! function used
  integer*4 modified_julday, pointer_string
  character*3 prn_mat(MAXSAT)
  integer*4 prn_int
  
  call prn_matbld(prn_mat)
  call frequency_glonass(FREQ1_R,FREQ2_R)
  obs_prio_index = 0
  biasW_index_G = index(obs_prio_G, 'W')
  biasX_index_E = index(obs_prio_E, 'X')
  biasI_index_C = index(obs_prio_C, 'I')
  biasL_index_J = index(obs_prio_J, 'L')

  ierr = 0
  line = ' '
  prn = ''
10 continue  ! next record
  read (lfn, '(a)', end=200) line
  msg = ' '
!
!! start line
  if (line(1:1) .ne. '>') then   !!!! rinex 3.03
    goto 100
  endif
!
!! number of satellite
  read (line(33:35), '(i3)', iostat=ioerr) nprn
  if (ioerr .ne. 0) then
    msg = 'read satellite number error.'
  endif
  if (len_trim(msg) .ne. 0) goto 100
!
!! Check the RINEX 3 event flag
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
      endif
    enddo
    goto 10
  endif
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
  read (line, '(1x,i5,4i3,f11.7,i3,i3)', err=100) iy, im, id, ih, imi, sec, iflag, nprn
  read (line(42:56), '(f15.12)', iostat=ioerr) dt
  if (ioerr .ne. 0) dt = 0.d0
  if (nprn .gt. MAXSAT + 100) then
    write (*, '(a,i3)') '***ERROR(rdrnxoi): nprn > maxsat ', nprn
    call exit()
  endif
!
!! check time
  if (im .le. 0 .or. im .gt. 12 .or. id .le. 0 .or. id .gt. 31 .or. ih .lt. 0 .or. ih .ge. 24 &
      .or. imi .lt. 0 .or. imi .gt. 60 .or. sec .lt. 0.d0 .or. sec .gt. 60.d0) then
    msg = 'epoch time incorrect'
    goto 100
  endif
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
      i = nprn
      do j = 1, i
        read (lfn, '(a)') line
      enddo
      line = ' '
      goto 10
    endif
  endif

  ii = 0
!
!! read data. if more than 10 type 3 line should be merged to one
  do i = 1, nprn
    read (lfn, '(a)', err=100, end=200) string
    read (string,'(1x,i2)') prn_int
    if (string(1:1).eq.'R' .and. prn_int .gt. 30) cycle
    if (string(1:1).eq.'R' .and. (OB%glschn(prn_int) .gt. 50 .or. OB%glschn(prn_int) .lt. -50)) cycle
    if (string(1:1) .eq. 'G' .or. string(1:1).eq.'R' .or. string(1:1).eq.'E' .or. string(1:1).eq.'C' .or. string(1:1).eq.'J') then
      if(string(1:1) .eq. 'G') then
        nobstyp_tmp=HD%nobstyp3_G
      elseif(string(1:1) .eq. 'R') then
        nobstyp_tmp=HD%nobstyp3_R
      elseif(string(1:1) .eq. 'E') then
        nobstyp_tmp=HD%nobstyp3_E
      elseif(string(1:1) .eq. 'C') then
        nobstyp_tmp=HD%nobstyp3_C
      elseif(string(1:1) .eq. 'J') then
        nobstyp_tmp=HD%nobstyp3_J
      endif
      ii = ii + 1
      !
      !! check if the sallite is requested
      read (string, '(a3,50(f14.3,2x))', err=100) prn(ii), (obs(j), j=1, min(nobstyp_tmp, maxtyp))
      if(prn(ii)(1:1) .ne. ' ' .and. prn(ii)(2:2) .eq. ' ' .and. prn(ii)(3:3) .ne. ' ')prn(ii)(2:2)='0'
      sysid(ii)=prn(ii)(1:1)
      i0 = 0
      if (nprn0 .gt. 0) then
        do j = 1, nprn0
          if (prn0(j) .eq. prn(ii)) i0 = j
        enddo
      else
        i0 = ii
      endif
!
!! fill in the obs. structure
      if (i0 .ne. 0) then
        l1=0
        l2=0
        p1=0
        p2=0
        obs_prio_index = 0
        phs_prio_index = 0
        do j = 1, nobstyp_tmp
          prio_index = 0
          if (dabs(obs(j)) .lt. MAXWND) cycle
          if (nprn0 .gt. 0) then
            iii=pointer_string(MAXSAT,prn_mat,prn0(i0))
          else
            iii=pointer_string(MAXSAT,prn_mat,prn(i0))
          endif
          if (prn_mat(iii)(1:1).eq.'G') then
            prio_index = index(obs_prio_G, HD%obstyp3_G(j) (3:3))
            if(prio_index.eq.0) cycle
            if (HD%obstyp3_G(j) (1:2) .eq. 'L1') then
              if(bias(iii, prio_index) .eq. 0.d0) then
                if(phs_prio_index(1).lt.prio_index) then
                  l1=j
                  phs_prio_index(1)=prio_index
                endif
              else
                if(obs_prio_index(1).lt.prio_index) then
                  OB%obs(i0, 1) = obs(j) - bias(iii, prio_index)*freq1_G/vlight
                  obs_prio_index(1) = prio_index
                  OB%typuse(i0, 1) = HD%obstyp3_G(j)
                endif
              endif
            elseif (HD%obstyp3_G(j) (1:2) .eq. 'L2') then
              if(bias(iii, prio_index+9) .eq. 0.d0) then
                if(phs_prio_index(2).lt.prio_index) then
                  l2=j
                  phs_prio_index(2)=prio_index
                endif
              else
                if(obs_prio_index(2).lt.prio_index) then
                  OB%obs(i0, 2) = obs(j) - bias(iii, prio_index+9)*freq2_G/vlight
                  obs_prio_index(2) = prio_index
                  OB%typuse(i0, 2) = HD%obstyp3_G(j)
                endif
              endif
            elseif (HD%obstyp3_G(j) (1:2) .eq. 'C1') then
              if(bias(iii, prio_index+9*2).eq.0.d0) then
                if(phs_prio_index(3).lt.prio_index) then
                  p1=j
                  phs_prio_index(3)=prio_index
                endif
              else
                if(obs_prio_index(3).lt.prio_index) then
                  OB%obs(i0, 3) = obs(j) - bias(iii, prio_index+9*2)
                  obs_prio_index(3) = prio_index
                  OB%typuse(i0, 3) = HD%obstyp3_G(j)
                endif
              endif
            elseif (HD%obstyp3_G(j) (1:2) .eq. 'C2') then
              if(bias(iii, prio_index+9*3).eq.0.d0) then
                if(phs_prio_index(4).lt.prio_index) then
                  p2=j
                  phs_prio_index(4)=prio_index
                endif
              else
                if(obs_prio_index(4).lt.prio_index) then
                  OB%obs(i0, 4) = obs(j) - bias(iii, prio_index+9*3)
                  obs_prio_index(4) = prio_index
                  OB%typuse(i0, 4) = HD%obstyp3_G(j)
                endif
              endif
            endif
          elseif (prn_mat(iii)(1:1).eq.'R') then
            prio_index = index(obs_prio_R, HD%obstyp3_R(j) (3:3))
            if(prio_index.eq.0) cycle
            read(prn_mat(iii),'(1x,i2)') prn_int
            frequency_glo_nu=OB%glschn(prn_int)
            if (HD%obstyp3_R(j) (1:2) .eq. 'L1') then
              if (obs_prio_index(1) .lt. prio_index) then
                obs_prio_index(1) = prio_index
                OB%obs(i0, 1) = obs(j)
                OB%typuse(i0, 1) = HD%obstyp3_R(j)
              endif
            elseif (HD%obstyp3_R(j) (1:2) .eq. 'L2') then
              if (obs_prio_index(2) .lt. prio_index) then
                obs_prio_index(2) = prio_index
                OB%obs(i0, 2) = obs(j)
                OB%typuse(i0, 2) = HD%obstyp3_R(j)
              endif
            elseif (HD%obstyp3_R(j) (1:2) .eq. 'C1') then
              if(bias(iii, prio_index+9*2).eq.0.d0) then
                if(phs_prio_index(3).lt.prio_index) then
                  p1=j
                  phs_prio_index(3)=prio_index
                endif
              else
                if(obs_prio_index(3) .lt. prio_index) then
                  OB%obs(i0, 3) = obs(j) - bias(iii, prio_index+9*2)
                  obs_prio_index(3) = prio_index
                  OB%typuse(i0, 3) = HD%obstyp3_R(j)
                endif
              endif
            elseif (HD%obstyp3_R(j) (1:2) .eq. 'C2') then
              if(bias(iii, prio_index+9*3).eq.0.d0) then
                if(phs_prio_index(4).lt.prio_index) then
                  p2=j
                  phs_prio_index(4)=prio_index
                endif
              else
                if(obs_prio_index(4) .lt. prio_index) then
                  OB%obs(i0, 4) = obs(j) - bias(iii, prio_index+9*3)
                  obs_prio_index(4) = prio_index
                  OB%typuse(i0, 4) = HD%obstyp3_R(j)
                endif
              endif
            endif
          elseif (prn_mat(iii)(1:1).eq.'E') then
            prio_index = index(obs_prio_E, HD%obstyp3_E(j) (3:3))
            if(prio_index.eq.0) cycle
            if (HD%obstyp3_E(j) (1:2) .eq. 'L1') then
              if(bias(iii, prio_index) .eq. 0.d0) then
                if(phs_prio_index(1).lt.prio_index) then
                  l1=j
                  phs_prio_index(1)=prio_index
                endif
              else
                if(obs_prio_index(1).lt.prio_index) then
                  OB%obs(i0, 1) = obs(j) - bias(iii, prio_index)*freq1_E/vlight
                  obs_prio_index(1) = prio_index
                  OB%typuse(i0, 1) = HD%obstyp3_E(j)
                endif
              endif
            elseif (HD%obstyp3_E(j) (1:2) .eq. 'L5') then
              if(bias(iii, prio_index+9) .eq. 0.d0) then
                if(phs_prio_index(2).lt.prio_index) then
                  l2=j
                  phs_prio_index(2)=prio_index
                endif
              else
                if(obs_prio_index(2).lt.prio_index) then
                  OB%obs(i0, 2) = obs(j) - bias(iii, prio_index+9)*freq2_E/vlight
                  obs_prio_index(2) = prio_index
                  OB%typuse(i0, 2) = HD%obstyp3_E(j)
                endif
              endif
            elseif (HD%obstyp3_E(j) (1:2) .eq. 'C1') then
              if(bias(iii, prio_index+9*2).eq.0.d0) then
                if(phs_prio_index(3) .lt. prio_index) then
                  p1=j
                  phs_prio_index(3)=prio_index
                endif
              else
                if(obs_prio_index(3) .lt. prio_index) then
                  OB%obs(i0, 3) = obs(j) - bias(iii, prio_index+9*2)
                  obs_prio_index(3) = prio_index
                  OB%typuse(i0, 3) = HD%obstyp3_E(j)
                endif
              endif
            elseif (HD%obstyp3_E(j) (1:2) .eq. 'C5') then
              if(bias(iii, prio_index+9*3).eq.0.d0) then
                if(phs_prio_index(4) .lt. prio_index) then
                  p2=j
                  phs_prio_index(4)=prio_index
                endif
              else
                if(obs_prio_index(4) .lt. prio_index) then
                  OB%obs(i0, 4) = obs(j) - bias(iii, prio_index+9*3)
                  obs_prio_index(4) = prio_index
                  OB%typuse(i0, 4) = HD%obstyp3_E(j)
                endif
              endif
            endif
          elseif (prn_mat(iii)(1:1).eq.'C') then
            prio_index = index(obs_prio_C, HD%obstyp3_C(j) (3:3))
            if(prio_index.eq.0) cycle
            if (HD%obstyp3_C(j) (1:2) .eq. 'L2') then
              if(bias(iii, prio_index) .eq. 0.d0) then
                if(phs_prio_index(1).lt.prio_index) then
                  l1=j
                  phs_prio_index(1)=prio_index
                endif
              else
                if(obs_prio_index(1).lt.prio_index) then
                  OB%obs(i0, 1) = obs(j) - bias(iii, prio_index)*freq1_C/vlight
                  obs_prio_index(1) = prio_index
                  OB%typuse(i0, 1) = HD%obstyp3_C(j)
                endif
              endif
            elseif (HD%obstyp3_C(j) (1:2) .eq. 'L6') then
              if(bias(iii, prio_index+9) .eq. 0.d0) then
                if(phs_prio_index(2).lt.prio_index) then
                  l2=j
                  phs_prio_index(2)=prio_index
                endif
              else
                if(obs_prio_index(2).lt.prio_index) then
                  OB%obs(i0, 2) = obs(j) - bias(iii, prio_index+9)*freq2_C/vlight
                  obs_prio_index(2) = prio_index
                  OB%typuse(i0, 2) = HD%obstyp3_C(j)
                endif
              endif
            elseif (HD%obstyp3_C(j) (1:2) .eq. 'C2') then
              if(bias(iii, prio_index+9*2).eq.0.d0) then
                if(phs_prio_index(3) .lt. prio_index) then
                  p1=j
                  phs_prio_index(3)=prio_index
                endif
              else
                if(obs_prio_index(3) .lt. prio_index) then
                  OB%obs(i0, 3) = obs(j) - bias(iii, prio_index+9*2)
                  obs_prio_index(3) = prio_index
                  OB%typuse(i0, 3) = HD%obstyp3_C(j)
                endif
              endif
            elseif (HD%obstyp3_C(j) (1:2) .eq. 'C6') then
              if(bias(iii, prio_index+9*3).eq.0.d0) then
                if(phs_prio_index(4) .lt. prio_index) then
                  p2=j
                  phs_prio_index(4)=prio_index
                endif
              else
                if(obs_prio_index(4) .lt. prio_index) then
                  OB%obs(i0, 4) = obs(j) - bias(iii, prio_index+9*3)
                  obs_prio_index(4) = prio_index
                  OB%typuse(i0, 4) = HD%obstyp3_C(j)
                endif
              endif
            endif
          elseif (prn_mat(iii)(1:1).eq.'J') then
            prio_index = index(obs_prio_J, HD%obstyp3_J(j) (3:3))
            if(prio_index.eq.0) cycle
            if (HD%obstyp3_J(j) (1:2) .eq. 'L1') then
              if(bias(iii, prio_index) .eq. 0.d0) then
                if(phs_prio_index(1).lt.prio_index) then
                  l1=j
                  phs_prio_index(1)=prio_index
                endif
              else
                if(obs_prio_index(1).lt.prio_index) then
                  OB%obs(i0, 1) = obs(j) - bias(iii, prio_index)*freq1_J/vlight
                  obs_prio_index(1) = prio_index
                  OB%typuse(i0, 1) = HD%obstyp3_J(j)
                endif
              endif
            elseif (HD%obstyp3_J(j) (1:2) .eq. 'L2') then
              if(bias(iii, prio_index+9) .eq. 0.d0) then
                if(phs_prio_index(2).lt.prio_index) then
                  l2=j
                  phs_prio_index(2)=prio_index
                endif
              else
                if(obs_prio_index(2).lt.prio_index) then
                  OB%obs(i0, 2) = obs(j) - bias(iii, prio_index+9)*freq2_J/vlight
                  obs_prio_index(2) = prio_index
                  OB%typuse(i0, 2) = HD%obstyp3_J(j)
                endif
              endif
            elseif (HD%obstyp3_J(j) (1:2) .eq. 'C1') then
              if(bias(iii, prio_index+9*2).eq.0.d0) then
                if(phs_prio_index(3) .lt. prio_index) then
                  p1=j
                  phs_prio_index(3)=prio_index
                endif
              else
                if(obs_prio_index(3) .lt. prio_index) then
                  OB%obs(i0, 3) = obs(j) - bias(iii, prio_index+9*2)
                  obs_prio_index(3) = prio_index
                  OB%typuse(i0, 3) = HD%obstyp3_J(j)
                endif
              endif
            elseif (HD%obstyp3_J(j) (1:2) .eq. 'C2') then
              if(bias(iii, prio_index+9*3).eq.0.d0) then
                if(phs_prio_index(4) .lt. prio_index) then
                  p2=j
                  phs_prio_index(4)=prio_index
                endif
              else
                if(obs_prio_index(4) .lt. prio_index) then
                  OB%obs(i0, 4) = obs(j) - bias(iii, prio_index+9*3)
                  obs_prio_index(4) = prio_index
                  OB%typuse(i0, 4) = HD%obstyp3_J(j)
                endif
              endif
            endif
          endif
        enddo
!! if one of the phases is zero, or the data is removed before
        if(OB%obs(i0,1).eq.0.d0) then
          OB%obs(i0,1)=obs(l1)
          if(prn_mat(iii)(1:1).eq.'G') then
            OB%typuse(i0,1)=HD%obstyp3_G(l1)
          elseif(prn_mat(iii)(1:1).eq.'E') then
            OB%typuse(i0,1)=HD%obstyp3_E(l1)
          elseif(prn_mat(iii)(1:1).eq.'C') then
            OB%typuse(i0,1)=HD%obstyp3_C(l1)
          elseif(prn_mat(iii)(1:1).eq.'J') then
            OB%typuse(i0,1)=HD%obstyp3_J(l1)
          endif
        endif
        if(OB%obs(i0,2).eq.0.d0) then
          OB%obs(i0,2)=obs(l2)
          if(prn_mat(iii)(1:1).eq.'G') then
            OB%typuse(i0,2)=HD%obstyp3_G(l2)
          elseif(prn_mat(iii)(1:1).eq.'E') then
            OB%typuse(i0,2)=HD%obstyp3_E(l2)
          elseif(prn_mat(iii)(1:1).eq.'C') then
            OB%typuse(i0,2)=HD%obstyp3_C(l2)
          elseif(prn_mat(iii)(1:1).eq.'J') then
            OB%typuse(i0,2)=HD%obstyp3_J(l2)
          endif
        endif
        if(OB%obs(i0,3).eq.0.d0) then
          OB%obs(i0,3)=obs(p1)
          if(prn_mat(iii)(1:1).eq.'G') then
            OB%typuse(i0,3)=HD%obstyp3_G(p1)
          elseif(prn_mat(iii)(1:1).eq.'R') then
            OB%typuse(i0,3)=HD%obstyp3_R(p1)
          elseif(prn_mat(iii)(1:1).eq.'E') then
            OB%typuse(i0,3)=HD%obstyp3_E(p1)
          elseif(prn_mat(iii)(1:1).eq.'C') then
            OB%typuse(i0,3)=HD%obstyp3_C(p1)
          elseif(prn_mat(iii)(1:1).eq.'J') then
            OB%typuse(i0,3)=HD%obstyp3_J(p1)
          endif
        endif
        if(OB%obs(i0,4).eq.0.d0) then
          OB%obs(i0,4)=obs(p2)
          if(prn_mat(iii)(1:1).eq.'G') then
            OB%typuse(i0,4)=HD%obstyp3_G(p2)
          elseif(prn_mat(iii)(1:1).eq.'R') then
            OB%typuse(i0,4)=HD%obstyp3_R(p2)
          elseif(prn_mat(iii)(1:1).eq.'E') then
            OB%typuse(i0,4)=HD%obstyp3_E(p2)
          elseif(prn_mat(iii)(1:1).eq.'C') then
            OB%typuse(i0,4)=HD%obstyp3_C(p2)
          elseif(prn_mat(iii)(1:1).eq.'J') then
            OB%typuse(i0,4)=HD%obstyp3_J(p2)
          endif
        endif
        if (any(OB%obs(i0, 1:4) .eq. 0.d0)) then
          OB%obs(i0, 1:4) = 0.d0
        endif
      endif
    else
      cycle
    endif
  enddo

  if (nprn0 .eq. 0) then
    OB%nprn = ii
    do i = 1, nprn
      if (sysid(i) .eq. 'G' .or. sysid(i) .eq. 'R' .or. sysid(i) .eq. 'E' .or. sysid(i) .eq. 'C' .or. sysid(i) .eq. 'J') then
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
  write (*, '(a/a/a)') '***ERROR(rdrnxoi): read file, '//trim(name), '   line :'//line(1:80), ' &
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
  write (*, '(a/a)') '###WARNING(rdrnxoi): end of file, ', trim(name)
  return
end
