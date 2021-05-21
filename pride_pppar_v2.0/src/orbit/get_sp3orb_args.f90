!
!! get_sp3orb_args.f90
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
!! Contributor: Maorong Ge, Jianghui Geng
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
  integer*4 nargs, lfn, i, iy, im, id, ih, nprn, ierr, year, doy, mjd, lfnleap, lfnbrdc, taiutc, taiutc0(100000)
  character*4 yyyy
  character*3 ddd
  character*20 brdcfil
  character*3 prn(MAXSAT)
  real*8 mjderp0, mjderp1
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
    write (*, '(a)') 'Usage: sp3orb sp3fil -cfg sescfg [-erp erpfil]'
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
!! time span for the erp
  lfn = get_valid_unit(10)
  open (lfn, file=erpfil, status='old', iostat=ierr)
  if (ierr .ne. 0) then
    write (*, '(2a)') '***ERROR(get_sp3orb_args): open file ', trim(erpfil)
    call exit(1)
  endif
  line = ' '
  do while (index(line, 'MJD') .eq. 0 .or. index(line, 'UT1') .eq. 0 .or. index(line, 'LOD') .eq. 0)
    read (lfn, '(a)', end=50) line
  enddo
  read (lfn, '(a)') line
  mjderp0 = 0.d0
  mjderp1 = 0.d0
  do while (.true.)
    read (lfn, '(a)', end=50) line
    if (mjderp0 .eq. 0.d0) read (line, *) mjderp0
    read (line, *) mjderp1
  enddo
50 close (lfn)
!
!! read sp3 header
  nprn = 0    ! read all satellites
  call rdsp3h(sp3fil, OH%jd0, OH%sod0, OH%jd1, OH%sod1, OH%dintv, nprn, prn)
  if (OH%jd0 + OH%sod0/86400.d0 .lt. mjderp0) then
    OH%jd0 = int(mjderp0)
    OH%sod0 = (mjderp0 - OH%jd0)*86400.d0
  endif
  if (OH%jd1 + OH%sod1/86400.d0 .gt. mjderp1) then
    OH%jd1 = int(mjderp1)
    OH%sod1 = (mjderp1 - OH%jd1)*86400.d0
  endif
  do mjd=45150,59214
    if(mjd .gt. 45150 .and. mjd .le. 45515)then
      taiutc0(mjd)=21
    elseif(mjd .gt. 45515 .and. mjd .le. 46246)then
      taiutc0(mjd)=22
    elseif(mjd .gt. 46246 .and. mjd .le. 47160)then
      taiutc0(mjd)=23
    elseif(mjd .gt. 47160 .and. mjd .le. 47891)then
      taiutc0(mjd)=24
    elseif(mjd .gt. 47891 .and. mjd .le. 48256)then
      taiutc0(mjd)=25
    elseif(mjd .gt. 48256 .and. mjd .le. 48803)then
      taiutc0(mjd)=26
    elseif(mjd .gt. 48803 .and. mjd .le. 49168)then
      taiutc0(mjd)=27
    elseif(mjd .gt. 49168 .and. mjd .le. 49533)then
      taiutc0(mjd)=28
    elseif(mjd .gt. 49533 .and. mjd .le. 50082)then
      taiutc0(mjd)=29
    elseif(mjd .gt. 50082 .and. mjd .le. 50629)then
      taiutc0(mjd)=30
    elseif(mjd .gt. 50629 .and. mjd .le. 51178)then
      taiutc0(mjd)=31
    elseif(mjd .gt. 51178 .and. mjd .le. 53735)then
      taiutc0(mjd)=32
    elseif(mjd .gt. 53735 .and. mjd .le. 54831)then
      taiutc0(mjd)=33
    elseif(mjd .gt. 54831 .and. mjd .le. 56108)then
      taiutc0(mjd)=34
    elseif(mjd .gt. 56108 .and. mjd .le. 57203)then
      taiutc0(mjd)=35
    elseif(mjd .gt. 57203 .and. mjd .le. 57753)then
      taiutc0(mjd)=36
    elseif(mjd .gt. 57753 .and. mjd .le. 59214)then
      taiutc0(mjd)=37
    endif
  enddo
  lfnleap = get_valid_unit(10)
  open (lfnleap, file='leap.sec')
  write (lfnleap,'(a9)') "+leap sec"
  do mjd=OH%jd0,OH%jd1
    if (mjd .le. 59214) then !2020366
        write (lfnleap,'(2i10)') mjd, taiutc0(mjd)
    else
      call mjd2doy(mjd,year,doy)
      write (yyyy, '(i4.4)') year
      write (ddd, '(i3.3)') doy
      lfnbrdc = get_valid_unit(10)
      brdcfil = "brdc"//ddd//"0."//yyyy(3:4)//"n"
      open (lfnbrdc, file=brdcfil, status='old', iostat=ierr)
      if (ierr .ne. 0) then
        write (*, '(2a)') '***ERROR(get_sp3orb_args): open file ', trim(brdcfil)
        call exit(1)
      endif
      line = ''
      taiutc = 0
      do while (index(line, 'LEAP SECONDS') .eq. 0)
        read (lfnbrdc, '(a)') line
      enddo
      read (line, *) taiutc
      if (taiutc .ne. 0) then
        taiutc = taiutc + 19
        write (lfnleap,'(2i10)') mjd, taiutc
      endif
      close (lfnbrdc)
    endif
  enddo
  write (lfnleap,'(a9)') "-leap sec"
  close (lfnleap)
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
