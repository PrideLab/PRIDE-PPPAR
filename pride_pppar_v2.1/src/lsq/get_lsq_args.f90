!
!! get_lsq_args.f90
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
!! Contributor: Maorong Ge, Jianghui Geng, Songfeng Yang
!! 
!!
!!
!! purpose  : get arguments and read options for lsq
!! parameter:
!!   output : LCF  -- lsq configure options
!!            SITE -- station infomation
!!            SAT  -- satellite information
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

  type(lsqcfg) LCF
  type(station) SITE
  type(rnxobr) OB
  type(satellite) SAT(MAXSAT)
  type(ionex) IM
!
!! local
  logical*1 lexist
  integer*4 nargs, lfn, i, j, k, l, iunit, nfc, ierr
  character*3 iprn
  integer*4 i1, i2, i3
  integer*4 iy, imon, id, ih, imin, idoy
  character*4 yyyy
  character*3 ddd
  real*8 is, seslen, var
  character*30 sesfil
  character*256 msg, key, rnxpath
  character*256 path
  integer*4 flag, rename
  type(orbhdr) OH
  character*1 tmpsys_G, tmpsys_R, tmpsys_E, tmpsys_C, tmpsys_3, tmpsys_J
  integer*4 prn_int
  character*3 prn_str
!
!! function called
  integer*4 get_valid_unit, modified_julday, pointer_string, day_of_year
  real*8 timdif
  integer*4 rdglschn
  character*256 findkey, lower_string
  character*4 upper_string
!
!! read arguments
  nargs = iargc()
  if (nargs .lt. 1) then
    write (*, '(a)') 'Usage: lsq sesfil'
    call exit(4)
  endif
  call getarg(1, sesfil)
  lfn = get_valid_unit(10)
  open (lfn, file=sesfil, status='old', iostat=ierr)
  if (ierr .ne. 0) then
    write (*, '(2a)') '***ERROR(get_lsq_args): open file ', trim(sesfil)
    call exit(1)
  endif
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

  do i=1,30
    write(prn_str,'(a1,i2.2)') 'R',i
    OB%glschn(i)=rdglschn(prn_str,LCF%jd0,LCF%sod0)
  enddo
!
!! define file names
  call file_name(.false., 'orb', ' ', iy, imon, id, ih, LCF%flnorb)
  call file_name(.false., 'sck', ' ', iy, imon, id, ih, LCF%flnsck)
  call file_name(.false., 'rck', ' ', iy, imon, id, ih, LCF%flnrck)
  call file_name(.false., 'ztd', ' ', iy, imon, id, ih, LCF%flnztd)
  call file_name(.false., 'htg', ' ', iy, imon, id, ih, LCF%flnhtg)
  call file_name(.false., 'amb', ' ', iy, imon, id, ih, LCF%flnamb)
  call file_name(.false., 'res', ' ', iy, imon, id, ih, LCF%flnres)
  call file_name(.false., 'con', ' ', iy, imon, id, ih, LCF%flncon)
  call file_name(.false., 'neq', ' ', iy, imon, id, ih, LCF%flnneq)
  call file_name(.false., 'vmf', ' ', iy, imon, id, ih, LCF%flnvmf)
  call file_name(.false., 'fcb', ' ', iy, imon, id, ih, LCF%flnfcb)
  call file_name(.false., 'tec', ' ', iy, imon, id, ih, LCF%flntec)
  call file_name(.false., 'att', ' ', iy, imon, id, ih, LCF%flnatt)
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
  msg = 'Remove ambiguity'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .eq. 'EMPTY') goto 100
  LCF%lrmbias = .true.
  if (key(1:1) .ne. 'Y') LCF%lrmbias = .false.
!
!! Ambiguity fixing or not
  msg = 'Ambiguity fixing'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .eq. 'EMPTY') goto 100
  LCF%fcbuse = .true.
  if (key(1:1) .eq. 'N') LCF%fcbuse = .false.
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
!! Rinex directory
  msg = 'Rinex directory'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .eq. 'EMPTY') goto 100
  rnxpath = trim(key)
  i = len_trim(rnxpath)
  if (rnxpath(i:i) .ne. '/') then
    rnxpath(i + 1:i + 1) = '/'
  endif
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
!! GNSS satellites
  rewind lfn
  msg = '+GNSS satellites'
  key = ' '
  do while (key(1:16) .ne. msg(1:16))
    read (lfn, '(a)', end=100) key
  enddo
  LCF%nprn = 0
  call rdorbh(LCF%flnorb, iunit, OH)
!! check orbit span
  if (timdif(OH%jd0, OH%sod0, LCF%jd0, LCF%sod0) .gt. MAXWND) then
    LCF%jd0 = OH%jd0
    LCF%sod0 = OH%sod0
    write (*, '(a,i5,f8.1)') '###WARNING(get_lsq_args): data starting time truncated due to orbits ', &
      LCF%jd0, LCF%sod0
  endif
  if (timdif(OH%jd1, OH%sod1, LCF%jd1, LCF%sod1) .lt. -MAXWND) then
    LCF%jd1 = OH%jd1
    LCF%sod1 = OH%sod1
    write (*, '(a,i5,f8.1)') '###WARNING(get_lsq_args): data ending time truncated due to orbits ', &
      LCF%jd1, LCF%sod1
  endif
!! read satellites
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
    SAT(LCF%nprn)%sys = iprn(1:1)
    if(SAT(LCF%nprn)%sys .eq. 'G')then
      SAT(LCF%nprn)%typ = 'BLOCK'
    elseif(SAT(LCF%nprn)%sys .eq. 'R')then
      SAT(LCF%nprn)%typ='GLONASS'
    elseif(SAT(LCF%nprn)%sys .eq. 'E')then
      SAT(LCF%nprn)%typ='GALILEO'
    elseif(SAT(LCF%nprn)%sys .eq. 'C')then
      SAT(LCF%nprn)%typ='BEIDOU'
    elseif(SAT(LCF%nprn)%sys .eq. 'J')then
      SAT(LCF%nprn)%typ='QZSS'
    endif
    SAT(LCF%nprn)%xscf(1:3) = 0.d0
    SAT(LCF%nprn)%yscf(1:3) = 0.d0
    SAT(LCF%nprn)%zscf(1:3) = 0.d0
  enddo
  close (iunit)
!
!! read ionosphere map
  IM%bradius=6371.d0
  IM%height =450.d0
  if(LCF%lioh) then
    call rdionex(LCF%flntec, LCF%nprn, LCF%prn, IM)
  endif
!
!! Stations
  rewind lfn
  msg = '+Station used'
  key = ' '
  do while (key(1:13) .ne. msg(1:13))
    read (lfn, '(a)', end=100) key
  enddo
  do while (key(1:1) .ne. '-')
    read (lfn, '(a)', end=100) key
    if (key(1:1) .ne. ' ') cycle
    SITE%name = key(2:5)
    SITE%skd = key(7:8)
    SITE%map = key(10:12)
    read (key(13:), *, err=200) SITE%dclk0, SITE%cutoff, SITE%dztd0, &
      SITE%qztd, SITE%dhtg0, SITE%qhtg, SITE%sigr, &
      SITE%sigp, (SITE%dx0(i), i=1, 3)
    SITE%dclk0_G=SITE%dclk0
    SITE%dclk0_R=SITE%dclk0
    SITE%dclk0_E=SITE%dclk0
    SITE%dclk0_C=SITE%dclk0
    SITE%dclk0_3=SITE%dclk0
    SITE%dclk0_J=SITE%dclk0
    SITE%cutoff = SITE%cutoff*PI/180.d0
    SITE%undu = 0.d0
    SITE%rclock = 0.d0
    do i = 1, LCF%nprn
      SITE%first(i) = .true.
      SITE%prephi(i) = 0.d0
    enddo
!! receiver antenna pointer to atx
    SITE%iptatx = 0
!! file names
    SITE%name = lower_string(SITE%name)
    OB%lfnrhd = 0
    call file_name(.false., 'log', 'SNAM='//SITE%name, iy, imon, id, ih, OB%rhdfil)
    call file_name(.false., 'pos', 'SNAM='//SITE%name, iy, imon, id, ih, LCF%flnpos)
    if (index(SITE%skd, 'K') .ne. 0) then     ! for kinematic & pseudo-kinematic use
      SITE%ikin = 0
      call file_name(.false., 'kin', 'SNAM='//SITE%name, iy, imon, id, ih, SITE%kinfil)
    endif

    SITE%iunit = 0
    SITE%imet = 0
    i = len_trim(rnxpath)
    SITE%obsfil = rnxpath(1:i)
    call file_name(.false., 'rnxo', 'SNAM='//SITE%name, iy, imon, id, ih, SITE%obsfil(i + 1:))
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
            endif
          enddo
        enddo
      enddo
      write (*, '(2a)') '###WARNING(get_lsq_args): observation not exist ', &
                        SITE%obsfil(:i - 12)//upper_string(SITE%name)//"00"//"XXX"//"_R_"//yyyy//ddd//"0000_01D_30S_MO.rnx"
      call exit(1)
    endif
!! receiver clock jump file
50  SITE%lfnjmp = 0
    path = '.'//SITE%name//'.jmp'
    inquire(file=path, exist=lexist)
    if (lexist) then
      SITE%lfnjmp = get_valid_unit(10)
      open(SITE%lfnjmp, file=path, status='old')
      write(*,'(2a)') '###INFO(get_lsq_args): read clock jump file: ', path
    end if
!! read position
    SITE%rlat = 0.d0
    SITE%rlon = 0.d0
    call read_position(LCF%flnpos, SITE%name, SITE%skd, SITE%x, SITE%dx0)
    if (all(SITE%x(1:3) .eq. 1.d0)) then
      write (*, '(a,a4)') '###WARNING(get_lsq_args): no position ', SITE%name
      call exit(1)
    else
      call xyzblh(SITE%x(1:3)*1.d3, 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, SITE%geod)
      SITE%geod(3) = SITE%geod(3)*1.d-3
      call rot_enu2xyz(SITE%geod(1), SITE%geod(2), SITE%rot_l2f)
      if(index(LCF%tide,'OCEAN').ne.0) then
        call oceanload_coef(LCF%otluse, SITE%name, &
                            SITE%geod(1), SITE%geod(2), SITE%rlat, &
                            SITE%rlon, SITE%olc)
      endif
    endif
!
!! add vmf3 grid
    if (SITE%map(1:3) .eq. 'VM1') then
      call vmf1_grid(LCF%flnvmf, SITE)
    else if (SITE%map(1:3) .eq. 'VM3') then
      call vmf3_grid(LCF%flnvmf, SITE)
    endif
!
!! check availability of station
    call read_obsrhd(0, LCF%dintv, LCF%nprn, LCF%prn, OB)
    if (OB%ava_obs .lt. OB%rem_obs) then
      write (*, '(2a)') '###WARNING(get_lsq_args): bad observation quality ', trim(SITE%obsfil)
      write (*, '(2a,i12)') '###WARNING(get_lsq_args): ', 'obs available: ', OB%ava_obs
      !call exit(1)
    endif
    write (*, '(a4,1x,a4,1x,a,f7.1,a1)') 'STA:', SITE%name, trim(SITE%obsfil), &
      OB%ava_obs*1.d0/(OB%ava_obs + OB%rem_obs)*100.d0, '%'
  enddo
  close (lfn)
  !! add multisystem-GREC
  !
  tmpsys_G=''
  tmpsys_R=''
  tmpsys_E=''
  tmpsys_C=''
  tmpsys_3=''
  tmpsys_J=''
  do i=1,LCF%nprn
    if(LCF%prn(i)(1:1) .eq. 'G')then
      tmpsys_G='G'
    elseif(LCF%prn(i)(1:1) .eq. 'R')then
      tmpsys_R='R'
    elseif(LCF%prn(i)(1:1) .eq. 'E')then
      tmpsys_E='E'
    elseif(LCF%prn(i)(1:1) .eq. 'C')then
      read(LCF%prn(i)(2:3),'(i2)')prn_int
      if(prn_int .le. 17)then
        tmpsys_C='C'
      else
        tmpsys_3='3'
      endif
    elseif(LCF%prn(i)(1:1) .eq. 'J')then
      tmpsys_J='J'
    endif
  enddo
  LCF%sys=tmpsys_G//tmpsys_R//tmpsys_E//tmpsys_C//tmpsys_3//tmpsys_J
  LCF%sysnum=0
  if(tmpsys_G .eq. 'G')LCF%sysnum=LCF%sysnum+1
  if(tmpsys_R .eq. 'R')LCF%sysnum=LCF%sysnum+1
  if(tmpsys_E .eq. 'E')LCF%sysnum=LCF%sysnum+1
  if(tmpsys_C .eq. 'C')LCF%sysnum=LCF%sysnum+1
  if(tmpsys_3 .eq. '3')LCF%sysnum=LCF%sysnum+1
  if(tmpsys_J .eq. 'J')LCF%sysnum=LCF%sysnum+1
  !
  !! add multisystem-GREC
  return
100 continue
  write (*, '(3a)') '***ERROR(get_lsq_args): find option ', trim(msg), trim(key)
  call exit(1)
200 continue
  write (*, '(3a)') '***ERROR(get_lsq_args): read option ', trim(msg), trim(key)
  call exit(1)
end
