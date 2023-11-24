!
!! get_sp3orb_args.f90
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
!! Contributor: Maorong Ge, Jianghui Geng, Jihang Lin
!! 
!!
!
subroutine get_sp3orb_args(sescfg, sp3fil, orbfil, erpfil, OH)
  implicit none
  include '../header/const.h'
  include '../header/orbit.h'

  character*(*) sescfg, sp3fil, orbfil, erpfil
  type(orbhdr) OH
!
!! local
  integer*4 nargs, lfn, i, iy, im, id, ih, nprn, ierr
  character*3 prn(MAXSAT)
  character*256 :: line = '', msg = '', key = ''
!
!! functions called
  integer*4 get_valid_unit, pointer_string
  character*256 :: findkey
!
!! initialization
  OH%nprn = 0
  OH%jd0 = 0
  OH%jd1 = 0
!
!! read command arguments
  nargs = iargc()
  if (nargs .eq. 0) then
    write (*, '(a)') 'Usage: sp3orb sp3fil -cfg sescfg'
    call exit(4)
  endif
  sp3fil = ' '
  call getarg(1, sp3fil)
  i = 2
  sescfg = ' '
  erpfil = 'igserp'
  do while (i .le. nargs)
    call getarg(i, line)
    i = i + 1
    if (line(1:4) .eq. '-cfg') then
      call getarg(i, sescfg)
      i = i + 1
    else if (line(1:4) .eq. '-erp') then
      erpfil = ' '
      call getarg(i, erpfil)
      i = i + 1
    endif
  enddo
!
!! read sp3 header
  nprn = 0    ! read all satellites
  call rdsp3h(sp3fil, OH%jd0, OH%sod0, OH%jd1, OH%sod1, OH%dintv, nprn, prn)
!
!! read session configure
  lfn = get_valid_unit(10)
  open (lfn, file=sescfg, status='old', iostat=ierr)
  if (ierr .ne. 0) then
    write (*, '(2a)') '***ERROR(get_sp3orb_args): open file ', trim(sescfg)
    call exit(1)
  endif
!
!! time tag
  msg = 'Session time'
  key = findkey(lfn, msg, ' ')
  if (key(1:5) .eq. 'EMPTY') goto 100
  read (key, *, err=200) iy, im, id, ih
  call yr2year(iy)
!
!! orbfil name
  orbfil = ' '
  call file_name(.false., 'orb', ' ', iy, im, id, ih, orbfil)
!
!! satellite infomation
  rewind lfn
  msg = '+GNSS satellites'
  key = ' '
  line = ' '
  do while (line(1:16) .ne. msg(1:16))
    read (lfn, '(a)', end=100) line
  enddo
  i = 0
  do while (key(1:1) .ne. '-')
    read (lfn, '(a)', end=100) key
    if (key(1:1) .ne. ' ') cycle
    if(key(2:2).ne.'G' .and. key(2:2).ne.'R' .and. key(2:2).ne.'E' .and. key(2:2).ne.'C' .and. key(2:2).ne.'J') cycle
    i = i + 1
    read (key, '(1x,a3)', err=200) OH%prn(i)
! if nonexistent in sp3, then remove request for this satellite
    if (pointer_string(nprn, prn, OH%prn(i)) .eq. 0) i = i - 1
  enddo
  OH%nprn = i
!
!! satellite type
  OH%sattyp = 'GNSS'

  close (lfn)
  return
100 continue
  write (*, '(3a)') '***ERROR(get_sp3orb_args): find option ', trim(msg), trim(key)
  call exit(1)
200 continue
  write (*, '(3a)') '***ERROR(get_sp3orb_args): read option ', trim(msg), trim(key)
  call exit(1)
end
