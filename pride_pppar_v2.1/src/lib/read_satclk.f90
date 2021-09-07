!
!! read_satclk.f90
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
!! purpose   : read and interpolate satellite clock for rinex clock files
!!
!! parameter :
!!    input  : clkfil -- clock file name
!!             iprn   -- satellite prn
!!             jd,sod -- time requested
!!    output : jdc,sodc -- referece time
!!             x0     -- clock error
!!             x1     -- clock drift
!!             iflag  -- flag for clock existing
!!
!
subroutine read_satclk(clkfil, iprn, jd, sod, jdc, sodc, x0, x1, iflag)
  implicit none
  include '../header/const.h'

  integer*4 jd, jdc, iflag
  character*3 iprn
  real*8 sod, sodc, x0, x1
  character*(*) clkfil
!
!! local
  logical*1 lfirst
  integer*4 i, j, k, lfn, nprn(2), ierr
  character*3 ii
  character*3 prn(maxsat, 2)
  integer*4 jdf(2), jdx, iy, imon, id, ih, im
  real*8 sec, sodf(2), sodx, dt1, dt2, alpha
  real*8 dintv, a0(2, maxsat), a1(2, maxsat), a2(2, maxsat), c(6)
  character*256 line
!
!! function used
  integer*4 get_valid_unit, modified_julday, pointer_string
  real*8 timdif

  data lfirst/.true./
  save lfirst, lfn, jdf, sodf, nprn, prn, a0, a1, a2, dintv
  
!
!! first enter and open clock file
  iflag = 0
  if (lfirst) then
    lfirst = .false.
    lfn = get_valid_unit(10)
    open (lfn, file=clkfil, iostat=ierr)
    if (ierr .ne. 0) then
      write (*, '(2a)') '***ERROR(read_satclk) : open file ', trim(clkfil)
      call exit(1)
    endif
    dintv = 30.d0
    line = ' '
    do while (index(line, 'END OF HEADER') .eq. 0)
      read (lfn, '(a)', end=200) line
      if (index(line, 'INTERVAL') .ne. 0) read (line, *) dintv
    enddo
!
!! read two epoch clocks
    k = 1
    jdf = 0
    nprn = 0
    do while (k .le. 2)
      line = ' '
      do while (line(1:3) .ne. 'AS ')
        read (lfn, '(a)', end=200, err=100) line
      enddo
      do while (line(1:3) .eq. 'AS ')
        c = 0.d0
        read (line(4:), *, iostat=ierr) ii, iy, imon, id, ih, im, sec, j, c(1), c(2)
!
!! check time tag
        call yr2year(iy)
        jdx = modified_julday(id, imon, iy)
        sodx = ih*3600.d0 + im*60.d0 + sec
        if (jdf(k) .eq. 0) then
          jdf(k) = jdx
          sodf(k) = sodx
        else if (timdif(jdx, sodx, jdf(k), sodf(k)) .gt. MAXWND) then
          backspace lfn
          goto 202
        endif
!
!! read continuing lines
        if (j .eq. 4) then
          read (lfn, '(a)', end=200, err=100) line
          read (line, *, iostat=ierr) c(3), c(4)
        else if (j .eq. 6) then
          read (lfn, '(a)', end=200, err=100) line
          read (line, *, iostat=ierr) c(3), c(4), c(5), c(6)
        endif
!
!! set prn and clocks
        nprn(k) = nprn(k) + 1
        prn(nprn(k), k) = ii
        a0(k, nprn(k)) = c(1)
        a1(k, nprn(k)) = c(3)
        a2(k, nprn(k)) = c(5)
        read (lfn, '(a)', end=200, err=100) line
      enddo
202   k = k + 1
    enddo
    dintv = timdif(jdf(2), sodf(2), jdf(1), sodf(1))
    write (*, '(a,f7.1)') '%%%MESSAGE(read_satclk): satellite clock interval ', dintv
  endif
!
!! check time tag
10 dt1 = timdif(jd, sod, jdf(1), sodf(1))
  dt2 = timdif(jd, sod, jdf(2), sodf(2))
  if (dt1 .lt. -MAXWND) then
    write (*, '(a)') '***ERROR(read_satclk) : read_satclk t < trefclk '
    iflag = 1
    return
  else if (dt2 .gt. MAXWND) then
!
!! transfer clocks
    nprn(1) = nprn(2)
    do i = 1, nprn(2)
      prn(i, 1) = prn(i, 2)
      a0(1, i) = a0(2, i)
      a1(1, i) = a1(2, i)
      a2(1, i) = a2(2, i)
    enddo
    jdf(1) = jdf(2)
    sodf(1) = sodf(2)
    jdf(2) = 0
    nprn(2) = 0
!
!! read next epoch clocks
    line = ' '
    do while (line(1:3) .ne. 'AS ')
      read (lfn, '(a)', end=200, err=100) line
    enddo
    do while (line(1:3) .eq. 'AS ')
      c = 0.d0
      read (line(4:), *, iostat=ierr) ii, iy, imon, id, ih, im, sec, j, c(1), c(2)
!
!! check time tag
      call yr2year(iy)
      jdx = modified_julday(id, imon, iy)
      sodx = ih*3600.d0 + im*60.d0 + sec
      if (jdf(2) .eq. 0) then
        jdf(2) = jdx
        sodf(2) = sodx
      else if (timdif(jdx, sodx, jdf(2), sodf(2)) .gt. MAXWND) then
        backspace lfn
        exit
      endif
!
!! read continuing lines
      if (j .eq. 4) then
        read (lfn, '(a)', end=200, err=100) line
        read (line, *, iostat=ierr) c(3), c(4)
      else if (j .eq. 6) then
        read (lfn, '(a)', end=200, err=100) line
        read (line, *, iostat=ierr) c(3), c(4), c(5), c(6)
      endif
!
!! set prn and clocks
      nprn(2) = nprn(2) + 1
      prn(nprn(2), 2) = ii
      a0(2, nprn(2)) = c(1)
      a1(2, nprn(2)) = c(3)
      a2(2, nprn(2)) = c(5)
      read (lfn, '(a)', iostat=ierr) line
      if (ierr .ne. 0) exit
    enddo
    goto 10
  else
!
!! check prn exist
    j = pointer_string(nprn(1), prn(1, 1), iprn)
    k = pointer_string(nprn(2), prn(1, 2), iprn)
!
!! set reference time
    if (dt1 .lt. 0.5d0*dintv) then
      jdc = jdf(1)
      sodc = sodf(1)
      if (j .ne. 0 .and. (a1(1, j) .ne. 0.d0 .or. k .eq. 0)) then
        x0 = a0(1, j)
        x1 = a1(1, j)
        return
      endif
    else
      jdc = jdf(2)
      sodc = sodf(2)
      if (k .ne. 0 .and. (a1(2, k) .ne. 0.d0 .or. j .eq. 0)) then
        x0 = a0(2, k)
        x1 = a1(2, k)
        return
      endif
    endif
!
!! interpolation
    if (j .ne. 0 .and. k .ne. 0) then
      alpha = (a0(2, k) - a0(1, j))/(timdif(jdf(2), sodf(2), jdf(1), sodf(1)))
      x0 = a0(1, j) + alpha*dt1
      x1 = 0.d0
    else
      iflag = 1
    endif
  endif

  return
100 write (*, '(2a)') '###WARNING(read_satclk) : read file ', trim(clkfil)
  iflag = 1
  return
200 write (*, '(2a)') '***ERROR(read_satclk) : end of file ', trim(clkfil)
  call exit(1)
end
