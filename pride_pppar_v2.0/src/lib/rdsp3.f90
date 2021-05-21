!
!! rdsp3h.f90
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
!! Contributor: Maorong Ge
!! 
!!
!!
!! purpose   : read IGS sp3 orbit file
!!
!! parameters: fln -- sp3 file name
!!             jd0,sod0,jd1,sod1,dintv -- start and stop time and interval
!!             nprn,prn -- number of satellites and satellite list
!!             jd,sod -- request epoch time
!!             xsp3   -- satellite position and velocity
!!
!
subroutine rdsp3h(fln, jd0, sod0, jd1, sod1, dintv, nprn, prn)
  implicit none
  include '../header/const.h'

  integer*4 jd, jd1, nprn
  character*3 prn(1:*)
  real*8 sod, sod1, dintv, xsp3(6, MAXSAT)
  character*(*) fln
!
!! local
  logical*1 lfirst
  integer*4 iy, im, id, ih, imin, nwk, mjd, i, jj, k, l, lfn, ierr, iflag
  character*3 j
  integer*4 nepo, jd0, nprn0, jdx
  character*3 prn0(170)
  real*8 fmjd, sec, sod0, x, y, z, t, dt, sodx
  character*80 :: line = ''
  character*1 :: sys(170), cid = ''
!
!! function called
  integer*4 modified_julday, pointer_string, get_valid_unit
  real*8 timdif

  data lfirst/.true./
  save lfirst, nprn0, lfn
!
!! first line with start time, # of epochs...
  lfn = 0
  if (lfirst) then
    lfn = get_valid_unit(10)
    open (lfn, file=fln, status='old', iostat=ierr)
    if (ierr .ne. 0) then
      write (*, '(2a)') '***ERROR(rdsp3h): open file ', trim(fln)
      call exit(1)
    endif
    jd0 = 0
    sod0 = 0.d0
!
!! scan sp3 file to load all gnss satellites
    do while (.true.)
!
!! second line with interval
      read (lfn, '(a)', end=33) line
      if (line(1:2) .ne. '#c' .and. line(1:2) .ne. '#d') cycle ! look for headers of sp3
      read (lfn, '(3x,i4,17x,f14.8,1x,i5,1x,f15.13)') nwk, dintv, mjd, fmjd
!
!! third to 12th lines with satellite prns
      read (lfn, '(2x,i4,3x,17(a3),9(/,9x,17(a3)))') nprn0, ( prn0(i), i=1, 170)
      do i=1,170
        sys(i)=prn0(i)(1:1)
      end do
      if (nprn0 .gt. 170) then
        write (*, '(a,i3)') '***ERROR(rdsp3h): too many satellites in SP3 file ', nprn0
        call exit(1)
      endif
!
!! only GNSS satellites are allowed
      do i = 1, nprn0
        if((sys(i).eq.'G'.or.sys(i).eq.'R'.or.sys(i).eq.'E'.or.sys(i).eq.'C'.or.sys(i).eq.'J')&
           .and.pointer_string(nprn,prn,prn0(i)).eq.0)then 
          nprn = nprn + 1
          prn(nprn) = prn0(i)
        endif
      enddo
    enddo
!
!! copy nprn to nprn0 for rdsp3i
33  nprn0 = nprn
!
!! time tags
    rewind lfn
    do while (.true.)
      read (lfn, '(a)', end=5) line
      if (line(1:1) .eq. '*') then
        read (line(2:), *) iy, im, id, ih, imin, sec
        if (jd0 .eq. 0) then
          jd0 = modified_julday(id, im, iy)
          sod0 = ih*3600.d0 + imin*60.d0 + sec
        endif
      endif
    enddo
5   rewind lfn
    jd1 = modified_julday(id, im, iy)
    sod1 = ih*3600.d0 + imin*60.d0 + sec
    lfirst = .false.
    return
  endif
!
!! read one record in SP3 file
!
  Entry rdsp3i(jd, sod, nprn, prn, xsp3, iflag)
  if (lfn .eq. 0) then
    write (*, '(a)') '***ERROR(rdsp3i): sp3fil not open '
    call exit(1)
  endif
!
!! check if all requested satellites in sp3-file
  iflag = 1
  do jj = 1, nprn0
    do i = 1, 3
      xsp3(i, jj) = 1.d15
    enddo
  enddo
!
!! read the required record
  k = 0
  line = ' '
  do while (line(1:1) .ne. '*')
    read (lfn, '(a)', end=100) line
  enddo
10 continue
  read (line(2:), *) iy, im, id, ih, imin, sec
  call yr2year(iy)
  jdx = modified_julday(id, im, iy)
  sodx = ih*3600.d0 + imin*60.d0 + sec
  dt = timdif(jdx, sodx, jd, sod)
  if (dabs(dt) .lt. MAXWND) then
    iflag = 0
    line = ' '
    do while (line(1:1) .ne. '*')
      read (lfn, '(a)', end=100) line
      if (line(1:1) .ne. 'P') cycle
      read (line(2:), '(a1)') cid
      read (line(2:), '(a3)') j
      if(cid.ne.'G'.and.cid.ne.'R'.and.cid.ne.'E'.and.cid.ne.'C'.and.cid.ne.'J') cycle
      read (line(5:), '(3f14.6)') x, y, z
!
!! nprn <= 0, get all the satellite in the file, otherwise, get data of the requested satellites
      if (nprn .gt. 0) then
        k = pointer_string(nprn, prn, j)
      else
        k = k + 1
        prn(k) = j
      endif
!
!! k == 0, this satellite is not in the requested list
      if (k .ne. 0) then
        xsp3(1, k) = x
        xsp3(2, k) = y
        xsp3(3, k) = z
      endif
    enddo
  else if (dt .lt. -MAXWND) then
    line = ' '
    do while (line(1:1) .ne. '*')
      read (lfn, '(a)', end=100) line
    enddo
    goto 10
  endif
  backspace lfn

100 continue
  if (nprn .le. 0) nprn = k
!
!! check if all prns have data
  do i = 1, nprn
    if (dabs(xsp3(1, i)) + dabs(xsp3(2, i)) + dabs(xsp3(3, i)) .eq. 0.d0) then
      xsp3(1, i) = 1.d15
      xsp3(2, i) = 1.d15
      xsp3(3, i) = 1.d15
    endif
  enddo
  return
!
!! close file
  Entry rdsp3c()
  close (lfn)
  lfn = 0
  return
end
