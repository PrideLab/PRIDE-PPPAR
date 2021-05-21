!
!! get_arsig_args.f90
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
!! purpose  : get arguments
!! parameter:
!!    output: FCB -- fractional cycle biases
!
subroutine get_arsig_args(FCB)
  implicit none
  include '../header/const.h'
  include 'abfcb.h'

  type(abfcb) FCB
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
  FCB%flnwlf = ' '
!
!! start & stop time
  msg = 'Session time'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .eq. 'EMPTY') goto 100
  read (key, *, err=200) iy, imon, id, ih, im, is, seslen
  FCB%jd0 = modified_julday(id, imon, iy)
  FCB%sod0 = ih*3600.d0 + im*60.d0 + is
  call timinc(FCB%jd0, FCB%sod0, seslen, FCB%jd1, FCB%sod1)
!
!! sampling rate
  msg = 'Interval'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .eq. 'EMPTY') goto 100
  read (key, *, err=200) FCB%dintv
!
!! define filename
  call file_name(.false., 'amb', ' ', iy, imon, id, ih, FCB%flnamb)
  call file_name(.false., 'neq', ' ', iy, imon, id, ih, FCB%flnneq)
  call file_name(.false., 'con', ' ', iy, imon, id, ih, FCB%flncon)
  call file_name(.false., 'fcb', ' ', iy, imon, id, ih, FCB%flnfcb)
!
!! GNSS satellites
  rewind lfn
  msg = '+GNSS satellites'
  key = ' '
  do while (key(1:16) .ne. msg(1:16))
    read (lfn, '(a)', end=100) key
  enddo
  FCB%nprn = 0
  do while (key(1:16) .ne. '-GNSS satellites')
    read (lfn, '(a)', end=100) key
    if (key(1:1) .ne. ' ') cycle
    read (key(2:4), *, err=200) iprn
    FCB%nprn = FCB%nprn + 1
    FCB%prn(FCB%nprn) = iprn
  enddo
!
!! Ambiguity search or not
  FCB%lsearch = .false.
  msg = 'Ambiguity fixing'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .eq. 'EMPTY') goto 100
  if (index(key, 'LAMBDA') .ne. 0) FCB%lsearch = .true.
!
!! Common observation time
  msg = 'Ambiguity duration'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .eq. 'EMPTY') goto 100
  read (key, *, err=200) FCB%minsec_common
!
!! Cutoff elevation for AR
  msg = 'Cutoff elevation'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .eq. 'EMPTY') goto 100
  read (key, *, err=200) FCB%cutoff
!
!! Maximum deleted & Minimum saved & Validation tests
  msg = 'Critical search'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .eq. 'EMPTY') goto 100
  read (key, *, err=200) FCB%maxdel, FCB%minsav, FCB%chisq, FCB%ratio
!
!! Bias fixing
  msg = 'Widelane decision'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .eq. 'EMPTY') goto 100
  read (key, *, err=200) FCB%wl_maxdev, FCB%wl_maxsig, FCB%wl_alpha

  msg = 'Narrowlane decision'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .eq. 'EMPTY') goto 100
  read (key, *, err=200) FCB%nl_maxdev, FCB%nl_maxsig, FCB%nl_alpha

  close (lfn)
  return
100 write (*, '(3a)') '***ERROR(get_arsig_args): find option ', trim(msg), trim(key)
  call exit(1)
200 write (*, '(3a)') '***ERROR(get_arsig_args): read option ', trim(msg), trim(key)
  call exit(1)
end
