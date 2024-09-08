!
!! get_lsq_args.f90
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
!! purpose  : get arguments and read options for lsq
!! parameter:
!!            output : LCF  -- lsq configure options
!!                     SITE -- station infomation
!!                     SAT  -- satellite information
!!
!
subroutine get_lsq_args(LCF, SITE, OB, SAT, IM)
  implicit none
  include '../header/const.h'
  include '../header/orbit.h'
  include '../header/station.h'
  include '../header/rnxobs.h'
  include '../header/satellite.h'
  include 'lsqcfg.h'
  include '../header/ionex.h'

! common
  integer*4     idxfrq(MAXSYS, 2)
  common        idxfrq

! parameter
  type(lsqcfg)  LCF
  type(station) SITE
  type(rnxobr)  OB
  type(satellite) SAT(MAXSAT)
  type(ionex)   IM

! local
  logical*1     lexist
  integer*4     nargs, lfn, i, j, k, l, iunit, nfc, ierr
  integer*4     jd, iy, imon, id, ih, imin, idoy
  real*8        is, seslen, var
  character*1   tmpsys
  character*4   yyyy
  character*3   iprn
  character*3   ddd
  character*3   prn_str
  character*30  sesfil
  character*256 msg, key, obsfil
  character*256 path
  type(orbhdr)  OH

! function called
  integer*4     get_valid_unit
  integer*4     modified_julday
  integer*4     pointer_string
  integer*4     day_of_year
  real*8        timdif
  character*256 findkey, lower_string
  character*4   upper_string
!
!! read arguments
  nargs = iargc()
  if (nargs .lt. 2) then
    write (*, '(a)') 'Usage: lsq sesfil rinex_obs_file'
    call exit(4)
  end if
  call getarg(1, sesfil)
  lfn = get_valid_unit(10)
  open (lfn, file=sesfil, status='old', iostat=ierr)
  if (ierr .ne. 0) then
    write (*, '(2a)') '***ERROR(get_lsq_args): open file ', trim(sesfil)
    call exit(1)
  end if
  call getarg(2, obsfil)
!
!! start & stop time
  msg = 'Session time'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .eq. 'EMPTY') goto 100
  read (key, *, err=200) iy, imon, id, ih, imin, is, seslen
  LCF%jd0 = modified_julday(id, imon, iy)
  LCF%sod0 = ih*3600.d0 + imin*60.d0 + is
  LCF%seslen = seslen
  call timinc(LCF%jd0, LCF%sod0, seslen, LCF%jd1, LCF%sod1)
!
!! get GLONASS channel number
  call read_glschn(LCF%jd0, LCF%sod0, OB%glschn)
!
!! define file names
  call file_name(.false., 'orb', ' ', iy, imon, id, ih, LCF%flnorb)
  call file_name(.false., 'sck', ' ', iy, imon, id, ih, LCF%flnsck)
  call file_name(.false., 'rck', ' ', iy, imon, id, ih, LCF%flnrck)
  call file_name(.false., 'ztd', ' ', iy, imon, id, ih, LCF%flnztd)
  call file_name(.false., 'htg', ' ', iy, imon, id, ih, LCF%flnhtg)
  call file_name(.false., 'amb', ' ', iy, imon, id, ih, LCF%flnamb)
  call file_name(.false., 'res', ' ', iy, imon, id, ih, LCF%flnres)
  call file_name(.false., 'cst', ' ', iy, imon, id, ih, LCF%flncon)
  call file_name(.false., 'neq', ' ', iy, imon, id, ih, LCF%flnneq)
  call file_name(.false., 'vmf', ' ', iy, imon, id, ih, LCF%flnvmf)
  call file_name(.false., 'fcb', ' ', iy, imon, id, ih, LCF%flnfcb)
  call file_name(.false., 'tec', ' ', iy, imon, id, ih, LCF%flntec)
  call file_name(.false., 'att', ' ', iy, imon, id, ih, LCF%flnatt)
  call file_name(.false., 'lat', ' ', iy, imon, id, ih, LCF%flnlat)
!
!!
  msg = 'Frequency combination'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .ne. 'EMPTY') then
    do i = 1, MAXSYS
      j = index(key, GNSS_PRIO(i:i))
      if (j .le. 0) cycle
      read (key(j + 1:j + 1), '(i1)', err=200) idxfrq(i, 1)
      read (key(j + 2:j + 2), '(i1)', err=200) idxfrq(i, 2)
    end do
  end if
  i = index(GNSS_PRIO, 'G')
  if (any(idxfrq(i, :) .eq. 0)) idxfrq(i, :) = [1, 2]
  i = index(GNSS_PRIO, 'R')
  if (any(idxfrq(i, :) .eq. 0)) idxfrq(i, :) = [1, 2]
  i = index(GNSS_PRIO, 'E')
  if (any(idxfrq(i, :) .eq. 0)) idxfrq(i, :) = [1, 5]
  i = index(GNSS_PRIO, 'C')
  if (any(idxfrq(i, :) .eq. 0)) idxfrq(i, :) = [2, 6]
  i = index(GNSS_PRIO, 'J')
  if (any(idxfrq(i, :) .eq. 0)) idxfrq(i, :) = [1, 2]
!
!!
  msg = 'Satellite orbit'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .eq. 'EMPTY') goto 100
  read (key, *, err=200) LCF%flnorb_real
!
!!
  msg = 'Satellite clock'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .eq. 'EMPTY') goto 100
  read (key, *, err=200) LCF%flnsck_real
!
!!
  msg = 'ERP'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .eq. 'EMPTY') goto 100
  read (key, *, err=200) LCF%flnerp_real
!
!!
  msg = 'Quaternions'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .eq. 'EMPTY') goto 100
  read (key, *, err=200) LCF%flnatt_real
  inquire (file=LCF%flnatt, exist=LCF%attuse)
!
!!
  msg = 'Code/phase bias'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .eq. 'EMPTY') goto 100
  read (key, *, err=200) LCF%flnfcb_real
  inquire (file=LCF%flnfcb, exist=LCF%fcbuse)
!
!!
  msg = 'LEO quaternions'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .eq. 'EMPTY') goto 100
  read (key, *, err=200) LCF%flnlat_real
  inquire (file=LCF%flnlat, exist=LCF%latuse)
!
!! erp filename
  LCF%flnerp = 'igserp'
!
!! sampling rate
  msg = 'Interval'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .eq. 'EMPTY') goto 100
  read (key, *, err=200) LCF%dintv
!
!!
  msg = 'Strict editing'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .eq. 'EMPTY') goto 100
  LCF%edit = trim(key)
!
!! pre-eliminate bias
  msg = 'Ambiguity co-var'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .eq. 'EMPTY') goto 100
  LCF%lrmbias = .true.
  if (key(1:1) .eq. 'Y') LCF%lrmbias = .false.
!
!! Ambiguity constraints
  LCF%nconG  = 0
  LCF%nconE  = 0
  LCF%nconC2 = 0
  LCF%nconC3 = 0
  LCF%nconJ  = 0
!
!! Tide displacement
  msg = 'Tides'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .eq. 'EMPTY') goto 100
  LCF%tide = trim(key)
!
!! ZTD model
  msg = 'ZTD model'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .eq. 'EMPTY') goto 100
  LCF%ztdmod = trim(key)
!
!! Horizontal Troposphere Gradients
  msg = 'HTG model'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .eq. 'EMPTY') goto 100
  LCF%htgmod = trim(key)
!
!! Receiver clock model
  msg = 'RCK model'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .eq. 'EMPTY') goto 100
  LCF%rckmod = trim(key)
!
!!
  msg = 'Iono 2nd'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .eq. 'EMPTY') goto 100
  LCF%lioh = .false.
  if (key(1:1) .eq. 'Y') LCF%lioh = .true.
!
!! Common observation time
  msg = 'Ambiguity duration'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .eq. 'EMPTY') goto 100
  read (key, *, err=200) LCF%minsec_common
!
!! Cutoff elevation for AR
  msg = 'Cutoff elevation'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .eq. 'EMPTY') goto 100
  read (key, *, err=200) LCF%cutoff
!
!! PCO on Melbourne-Wubbena
  msg = 'PCO on wide-lane'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .eq. 'EMPTY') goto 100
  LCF%pcowl = .true.
  if (key(1:1) .eq. 'N') LCF%pcowl = .false.
!
!! Maximum deleted & Minimum saved & Validation tests
  msg = 'Critical search'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .eq. 'EMPTY') goto 100
  read (key, *, err=200) LCF%maxdel, LCF%minsav, LCF%chisq, LCF%ratio
!
!! Bias fixing
  msg = 'Widelane decision'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .eq. 'EMPTY') goto 100
  read (key, *, err=200) LCF%wl_maxdev, LCF%wl_maxsig, LCF%wl_alpha

  msg = 'Narrowlane decision'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .eq. 'EMPTY') goto 100
  read (key, *, err=200) LCF%nl_maxdev, LCF%nl_maxsig, LCF%nl_alpha
!
!! Remove observations without bias corrections 
  LCF%del_nobia = (index(LCF%flnfcb, "WUM0MGXRTS") .gt. 0)
!
!! Truncate the ambiguities at day-boundary 
  msg = 'Truncate at midnight'
  key = findkey(lfn, msg, ' ')
  if (key(1:1) .eq. 'Y') then
    LCF%amb_at_dbd = 'TRUNCATED'
  else if (key(1:1) .eq. 'C') then
    LCF%amb_at_dbd = 'CONNECTED'
  else
    LCF%amb_at_dbd = 'DEFAULT'
  end if
!
!! GNSS satellites
  rewind lfn
  msg = '+GNSS satellites'
  key = ' '
  do while (key(1:16) .ne. msg(1:16))
    read (lfn, '(a)', end=100) key
  end do
  LCF%nprn = 0
  call rdorbh(LCF%flnorb, iunit, OH)
!! check orbit span
  if (timdif(OH%jd0, OH%sod0, LCF%jd0, LCF%sod0) .gt. MAXWND) then
    LCF%jd0 = OH%jd0
    LCF%sod0 = OH%sod0
    write (*, '(a,i5,f8.1)') '###WARNING(get_lsq_args): data starting time delayed due to orbits ', &
      LCF%jd0, LCF%sod0
  end if
  if (timdif(OH%jd1, OH%sod1, LCF%jd1, LCF%sod1) .lt. -MAXWND) then
    LCF%jd1 = OH%jd1
    LCF%sod1 = OH%sod1
    write (*, '(a,i5,f8.1)') '###WARNING(get_lsq_args): data ending time truncated due to orbits ', &
      LCF%jd1, LCF%sod1
  end if
!! read satellites
  LCF%sys = ' '
  do while (key(1:16) .ne. '-GNSS satellites')
    read (lfn, '(a)', end=100) key
    if (key(1:1) .ne. ' ') cycle
    read (key, *, err=200) iprn, var
    if (pointer_string(OH%nprn, OH%prn, iprn) .eq. 0) cycle
    LCF%nprn = LCF%nprn + 1
    LCF%prn(LCF%nprn) = iprn
    SAT(LCF%nprn)%prn = iprn
    SAT(LCF%nprn)%var = var
    SAT(LCF%nprn)%iptatx = 0
    tmpsys = iprn(1:1)
    if (tmpsys .eq. 'C') then
      read (iprn(2:3), '(i2)') l
      if (l .gt. 17) tmpsys = '3'
    end if
    if (index(LCF%sys, tmpsys) .eq. 0) then
      LCF%sys = LCF%sys(1:len_trim(LCF%sys))//tmpsys
    end if
    if (iprn(1:1) .eq. 'G') then
      SAT(LCF%nprn)%typ = 'BLOCK'
    elseif (iprn(1:1) .eq. 'R') then
      SAT(LCF%nprn)%typ = 'GLONASS'
    elseif (iprn(1:1) .eq. 'E') then
      SAT(LCF%nprn)%typ = 'GALILEO'
    elseif (iprn(1:1) .eq. 'C') then
      SAT(LCF%nprn)%typ = 'BEIDOU'
    elseif (iprn(1:1) .eq. 'J') then
      SAT(LCF%nprn)%typ = 'QZSS'
    end if
    SAT(LCF%nprn)%xscf(1:3) = 0.d0
    SAT(LCF%nprn)%yscf(1:3) = 0.d0
    SAT(LCF%nprn)%zscf(1:3) = 0.d0
  end do
  close (iunit)
!
!! read ionosphere map
  IM%bradius = 6371.d0
  IM%height = 450.d0
  if (LCF%lioh) then
    call rdionex(LCF%flntec, LCF%nprn, LCF%prn, IM)
  end if
!
!! Stations
  rewind lfn
  msg = '+Station used'
  key = ' '
  do while (key(1:13) .ne. msg(1:13))
    read (lfn, '(a)', end=100) key
  end do
  do while (key(1:1) .ne. '-')
    read (lfn, '(a)', end=100) key
    if (key(1:1) .ne. ' ') cycle
    SITE%name = key(2:5)
    SITE%skd = key(7:8)
    SITE%map = key(10:12)
    read (key(13:), *, err=200) SITE%dclk0, SITE%qrck, &
      SITE%cutoff, SITE%dztd0, &
      SITE%qztd, SITE%dhtg0, SITE%qhtg, SITE%sigr, &
      SITE%sigp, SITE%pospd, (SITE%dx0(i), i=1, 3), SITE%rx0
    SITE%cutoff = SITE%cutoff*PI/180.d0
    SITE%undu = 0.d0
  !  SITE%rx0 = 2.d-1    ! process noise for piece-wise coordinates
    do i = 1, LCF%nprn
      SITE%first(i) = .true.
      SITE%prephi(i) = 0.d0
    end do
!! receiver antenna pointer to atx
    SITE%iptatx = 0
!! file names
    SITE%name = lower_string(SITE%name)
    OB%lfnrhd = 0
    call file_name(.false., 'log', 'SNAM='//SITE%name, iy, imon, id, ih, OB%rhdfil)
    call file_name(.false., 'pos', 'SNAM='//SITE%name, iy, imon, id, ih, LCF%flnpos)
    if (SITE%skd(1:1) .eq. 'P' .or. SITE%skd(1:1) .eq. 'K' .or. SITE%skd(1:1) .eq. 'L') then
      SITE%ikin = 0
      call file_name(.false., 'kin', 'SNAM='//SITE%name, iy, imon, id, ih, SITE%kinfil)
    end if

    SITE%iunit = 0
    SITE%imet = 0
    i = len_trim(obsfil)
    SITE%obsfil = obsfil(1:i)
    inquire (file=SITE%obsfil, exist=lexist)
    i = len_trim(SITE%obsfil)
    if (.not. lexist) then
      write (*, '(2a)') '###WARNING(get_lsq_args): observation not exist ', trim(SITE%obsfil)
      do j = 65, 90
        do k = 65, 90
          do l = 65, 90
            idoy = day_of_year(id, imon, iy)
            write (yyyy, '(i4.4)') iy
            write (ddd, '(i3.3)') idoy
            SITE%obsfil(i - 11:) = &
              upper_string(SITE%name)//"00"//char(j)//char(k)//char(l)//"_R_"//yyyy//ddd//"0000_01D_30S_MO.rnx"
            inquire (file=SITE%obsfil, exist=lexist)
            if (lexist) then
              goto 50
            end if
          end do
        end do
      end do
      write (*, '(2a)') '###WARNING(get_lsq_args): observation not exist ', &
        SITE%obsfil(:i - 12)//upper_string(SITE%name)//"00"//"XXX"//"_R_"//yyyy//ddd//"0000_01D_30S_MO.rnx"
      call exit(1)
    end if
!! receiver clock jump file
50  SITE%lfnjmp = 0
    path = '.'//SITE%name//'.jmp'
    inquire (file=path, exist=lexist)
    if (lexist) then
      SITE%lfnjmp = get_valid_unit(10)
      open (SITE%lfnjmp, file=path, status='old')
      write (*, '(2a)') '###INFO(get_lsq_args): read clock jump file: ', path
    end if
!! read position
    SITE%rlat = 0.d0
    SITE%rlon = 0.d0
    call read_position(LCF%flnpos, SITE%name, SITE%skd, SITE%x, SITE%dx0)
    if (all(SITE%x(1:3) .eq. 1.d0)) then
      write (*, '(a,a4)') '###WARNING(get_lsq_args): no position ', SITE%name
      call exit(1)
    else if (SITE%skd(1:1) .ne. 'L') then
      call xyzblh(SITE%x(1:3)*1.d3, 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, SITE%geod)
      SITE%geod(3) = SITE%geod(3)*1.d-3
      call rot_enu2xyz(SITE%geod(1), SITE%geod(2), SITE%rot_l2f)
      if (index(LCF%tide, 'OCEAN') .ne. 0) then
        call oceanload_coef(LCF%otluse, SITE%name, &
                            SITE%geod(1), SITE%geod(2), SITE%rlat, &
                            SITE%rlon, SITE%olc)
        call file_name(.false., 'otl', 'SNAM='//SITE%name, iy, imon, id, ih,SITE%otlfil)
        inquire(file=SITE%otlfil, exist=lexist)
        if (.not. LCF%otluse .and. lexist) then
          LCF%otluse = .true.
        else if (.not. LCF%otluse .and. .not. lexist) then
          LCF%otluse = .false.
          write (*, '(2a)') '###WARNING(get_lsq_args): no oceanload model(Zhang) for ', SITE%name
        end if
      end if
    end if
!
!! add vmf3 grid
    if (SITE%map(1:3) .eq. 'VM1') then
      call vmf1_grid(LCF%flnvmf, SITE)
    else if (SITE%map(1:3) .eq. 'VM3') then
      call vmf3_grid(LCF%flnvmf, SITE)
    end if
!
!! check availability of station
    seslen = LCF%dintv
    call read_obsrhd(jd, seslen, LCF%nprn, LCF%prn, OB)
    if (OB%ava_obs .lt. OB%rem_obs) then
      write (*, '(2a)') '###WARNING(get_lsq_args): bad observation quality ', trim(SITE%obsfil)
    end if
    if (timdif(LCF%jd0, LCF%sod0, jd, seslen) .lt. -MAXWND) then
      LCF%jd0 = jd
      LCF%sod0 = seslen
      write (*, '(a,i5,f8.1)') '###WARNING(get_lsq_args): data starting time delayed due to log ', &
        LCF%jd0, LCF%sod0
    end if
    call timinc(jd, seslen, OB%dtrcv, jd, seslen)
    if (timdif(LCF%jd1, LCF%sod1, jd, seslen) .gt. MAXWND) then
      LCF%jd1 = jd
      LCF%sod1 = seslen
      write (*, '(a,i5,f8.1)') '###WARNING(get_lsq_args): data ending time truncated due to log ', &
        LCF%jd1, LCF%sod1
    end if
    write (*, '(a4,1x,a4,1x,a,f7.1,a1)') 'STA:', SITE%name, trim(SITE%obsfil), &
      OB%ava_obs*1.d0/(OB%ava_obs + OB%rem_obs)*100.d0, '%'
    exit
  end do
  close (lfn)
  return

100 continue
  write (*, '(a,2(1x,a))') '***ERROR(get_lsq_args): find option', trim(msg), trim(key)
  call exit(1)
200 continue
  write (*, '(a,2(1x,a))') '***ERROR(get_lsq_args): read option', trim(msg), trim(key)
  call exit(1)
end subroutine
