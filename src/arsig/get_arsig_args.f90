!
!! get_arsig_args.f90
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
!! purpose  : get arguments
!! parameter:
!!    output: ACF -- fractional cycle biases
!
subroutine get_arsig_args(ACF)
  implicit none
  include '../header/const.h'
  include 'arscfg.h'

  type(arscfg) ACF
!
!! local
  integer*4 nargs, i, lfn, iy, imon, id, ih, im, ierr
  character*3 iprn
  real*8 is, seslen
  character*20 sesfil
  character*256 msg, key, ctrl
!
!! function called
  integer*4 get_valid_unit, modified_julday
  character*256 findkey
!
!! read arguments
  nargs = iargc()
  if (nargs .eq. 0) then
    write (*, '(a)') 'Usage: arsig sesfil'
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
!! Ambiguity search or not
  ACF%lsearch = .false.
  msg = 'Ambiguity co-var'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .eq. 'EMPTY') goto 100
  if (index(key, 'NO') .eq. 0) ACF%lsearch = .true.
!
!! Common observation time
  msg = 'Ambiguity duration'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .eq. 'EMPTY') goto 100
  read (key, *, err=200) ACF%minsec_common
!
!! Cutoff elevation for AR
  msg = 'Cutoff elevation'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .eq. 'EMPTY') goto 100
  read (key, *, err=200) ACF%cutoff
!
!! Maximum deleted & Minimum saved & Validation tests
  msg = 'Critical search'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .eq. 'EMPTY') goto 100
  read (key, *, err=200) ACF%maxdel, ACF%minsav, ACF%chisq, ACF%ratio
!
!! Bias fixing
  msg = 'Widelane decision'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .eq. 'EMPTY') goto 100
  read (key, *, err=200) ACF%wl_maxdev, ACF%wl_maxsig, ACF%wl_alpha

  msg = 'Narrowlane decision'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .eq. 'EMPTY') goto 100
  read (key, *, err=200) ACF%nl_maxdev, ACF%nl_maxsig, ACF%nl_alpha

  close (lfn)
  return
100 write (*, '(3a)') '***ERROR(get_arsig_args): find option ', trim(msg), trim(key)
  call exit(1)
200 write (*, '(3a)') '***ERROR(get_arsig_args): read option ', trim(msg), trim(key)
  call exit(1)
end
