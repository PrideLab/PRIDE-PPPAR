!
!! get_arsig_args.f90
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
!! purpose  : get arguments
!! parameter:
!!    output: ACF -- fractional cycle biases
!
subroutine get_arsig_args(ACF)
  implicit none
  include '../header/const.h'
  include 'arscfg.h'

! common
  integer*4     idxfrq(MAXSYS, 2)
  common        idxfrq
! parameter
  type(arscfg)  ACF
! local
  integer*4     nargs, i, j, lfn, ierr
  integer*4     iy, imon, id, ih, im
  character*3   iprn
  real*8        is, seslen
  character*20  sesfil
  character*256 msg, key, ctrl
! function called
  integer*4     get_valid_unit
  integer*4     modified_julday
  character*256 findkey

!
!! read arguments
  nargs = iargc()
  if (nargs .eq. 0) then
    write (*, '(a)') "                                                            "
    write (*, '(a)') "arsig version 3.0, Wuhan University, Oct. 2024              "
    write (*, '(a)') "                                                            "
    write (*, '(a)') "Usage: arsig configfile                                     "
    write (*, '(a)') "                                                            "
    write (*, '(a)') "Description:                                                "
    write (*, '(a)') "  arsig is a module of PRIDE PPP-AR, is used for realizing  "
    write (*, '(a)') "  wide-lane and narrow-lane ambiguity fixed.                "
    write (*, '(a)') "                                                            "
    write (*, '(a)') "Required arguments:                                         "
    write (*, '(a)') "  configfile                                                "
    write (*, '(a)') "    PRIDE PPP-AR's config file.                             "
    write (*, '(a)') "                                                            "
    write (*, '(a)') "Note: Some dependent files need to be under folders.        "
    write (*, '(a)') "                                                            "
    write (*, '(a)') "Examples:                                                   "
    write (*, '(a)') "  arsig config_abmf                                         "
    write (*, '(a)') "                                                            "
    write (*, '(a)') "More details refer to PRIDE PPP-AR manual and repository    "
    write (*, '(a)') "  https://github.com/PrideLab/PRIDE-PPPAR/                  "
    write (*, '(a)') "                                                            "
    call exit(4)
  endif
  call getarg(1, sesfil)
  lfn = get_valid_unit(10)
  open (lfn, file=sesfil, status='old', iostat=ierr)
  if (ierr .ne. 0) then
    write (*, '(2a)') '***ERROR(get_arsig_args): open file ', trim(sesfil)
    call exit(1)
  endif
  ACF%flnwlf = ' '
!
!! frquency combination
  msg = 'Frequency combination'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .ne. 'EMPTY') then
    do i = 1, MAXSYS
      j = index(key, GNSS_PRIO(i:i))
      if (j .le. 0) cycle
      read (key(j+1:j+1), '(i1)', err=200) idxfrq(i, 1)
      read (key(j+2:j+2), '(i1)', err=200) idxfrq(i, 2)
    enddo
  endif
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
!! start & stop time
  msg = 'Session time'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .eq. 'EMPTY') goto 100
  read (key, *, err=200) iy, imon, id, ih, im, is, seslen
  ACF%jd0 = modified_julday(id, imon, iy)
  ACF%sod0 = ih*3600.d0 + im*60.d0 + is
  call timinc(ACF%jd0, ACF%sod0, seslen, ACF%jd1, ACF%sod1)
!
!! sampling rate
  msg = 'Interval'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .eq. 'EMPTY') goto 100
  read (key, *, err=200) ACF%dintv
!
!! define filename
  call file_name(.false., 'amb', ' ', iy, imon, id, ih, ACF%flnamb)
  call file_name(.false., 'neq', ' ', iy, imon, id, ih, ACF%flnneq)
  call file_name(.false., 'cst', ' ', iy, imon, id, ih, ACF%flncon)
!
!! GNSS satellites
  rewind lfn
  msg = '+GNSS satellites'
  key = ' '
  do while (key(1:16) .ne. msg(1:16))
    read (lfn, '(a)', end=100) key
  enddo
  ACF%nprn = 0
  do while (key(1:16) .ne. '-GNSS satellites')
    read (lfn, '(a)', end=100) key
    if (key(1:1) .ne. ' ') cycle
    read (key(2:4), *, err=200) iprn
    ACF%nprn = ACF%nprn + 1
    ACF%prn(ACF%nprn) = iprn
  enddo
!
!! ambiguity search or not
  ACF%lsearch = .false.
  msg = 'Ambiguity co-var'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .eq. 'EMPTY') goto 100
  if (index(key, 'NO') .eq. 0) ACF%lsearch = .true.
!
!! common observation time
  msg = 'Ambiguity duration'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .eq. 'EMPTY') goto 100
  read (key, *, err=200) ACF%minsec_common
!
!! cutoff elevation for AR
  msg = 'Cutoff elevation'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .eq. 'EMPTY') goto 100
  read (key, *, err=200) ACF%cutoff
!
!! maximum deleted & minimum saved & validation tests
  msg = 'Critical search'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .eq. 'EMPTY') goto 100
  read (key, *, err=200) ACF%maxdel, ACF%minsav, ACF%chisq, ACF%ratio
!
!! ambiguity fixing
  msg = 'Widelane decision'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .eq. 'EMPTY') goto 100
  read (key, *, err=200) ACF%wl_maxdev, ACF%wl_maxsig, ACF%wl_alpha
  msg = 'Narrowlane decision'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .eq. 'EMPTY') goto 100
  read (key, *, err=200) ACF%nl_maxdev, ACF%nl_maxsig, ACF%nl_alpha
!
!! verbose output or not 
  ACF%lverbose = .false.
  msg = 'Verbose output'
  key = findkey(lfn, msg, ' ')
  if (index(key, 'YES') .ne. 0) ACF%lverbose = .true.
!
!! Ambiguity validation 
  ACF%lambsvm = .false.
  msg = 'AI Ambiguity validation'
  key = findkey(lfn, msg, ' ')
  if (index(key, 'YES') .ne. 0) ACF%lambsvm = .true.

  close (lfn)
  return

100 write (*, '(3a)') '***ERROR(get_arsig_args): find option ', trim(msg), trim(key)
  call exit(1)

200 write (*, '(3a)') '***ERROR(get_arsig_args): read option ', trim(msg), trim(key)
  call exit(1)
end subroutine
