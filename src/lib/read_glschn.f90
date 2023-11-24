!
!! read_glschn.f90
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
!! Contributor: Jihang Lin
!! 
!!
!!
!! purpose  : read GLONASS channel number from sat_parameters or glonass_chn
!! parameter:
!!    input : mjd, sod -- MJD and seconds of day
!!    output: glschn -- GLONASS channel numbers
!
subroutine read_glschn(mjd, sod, glschn)
  implicit none
  include '../header/const.h'

  integer*4 mjd
  real*8 sod
  integer*4 glschn(MAXSAT_R)
!
!! local
  integer*4 lfnchn, ierr
  integer*4 isat, chn
  integer*4 year0, doy0, sod0
  integer*4 year1, doy1, sod1
  real*8 mjd0, mjd1
  character*256 fchn
  character*256 line
  character*3 svn
!
!! function used
  integer*4 get_valid_unit
  integer*4 modified_julday
!
!! read sat_parameters
  fchn = 'sat_parameters'
  lfnchn = get_valid_unit(10)
  open(lfnchn, file=fchn, status='old', iostat=ierr)
  if (ierr .ne. 0) then
    write (*, '(2a)') '###WARNING(read_glschn): unable to open file ', trim(fchn)
    call exit(1)
  endif
  read (lfnchn, '(a)') line
  do while (.true.)
    if (index(line, '+prn_indexed') .ne. 0) then
      glschn = 0
      do while (.true.)
 50     read (lfnchn, '(a)', iostat=ierr) line
        if (ierr .ne. 0) then
          backspace lfnchn
          exit
        endif
        if (line(1:1) .eq. '#') cycle
        if (line(2:2) .ne. 'R') cycle
        read (line, '(1x,a3,8x,2(i4,i3,1x,i5,1x),27x,i3)', err=200) &
          svn, year0, doy0, sod0, year1, doy1, sod1, chn
        read (svn, '(1x,i2)') isat
        if (isat .le. 0 .or. isat .gt. MAXSAT_R) then
          write (*, '(2a)') '###WARNING(read_glschn): invalid svn ', svn
          cycle
        endif
        mjd0 = modified_julday(doy0, 0, year0) + sod0/864.d2
        if ((year1+doy1) .eq. 0) then
          mjd1 = 99999
        else
          mjd1 = modified_julday(doy1, 0, year1) + sod1/864.d2
        endif
        if (mjd+sod/864.d2 .lt. mjd0 .or. mjd+sod/864.d2 .gt. mjd1) cycle
        glschn(isat) = chn
      enddo
    endif
    read (lfnchn, '(a)', end=100) line
  enddo
!
!! check end tag
100 continue
  if (index(line, '-prn_indexed') .eq. 0) then
    write (*, '(2a)') '***ERROR(read_glschn): end of file ', trim(fchn)
    call exit(1)
  else
    goto 150
  endif
!
!! check channel numbers
150 continue
  do isat = 1, MAXSAT_R
    if (abs(glschn(isat)) .gt. 50) then
      write (*, '(a,i3)') '###WARNINF(read_glschn): invalid GLONASS channel number ', glschn(isat)
      glschn(isat) = 99999
    endif
  enddo
  close(lfnchn)
  return

200 continue
  write (*, '(2a)') '***ERROR(read_glschn): read file ', trim(fchn)
  call exit(1)

300 continue
  write (*, '(a)') '***ERROR(read_glschn): no sat_parameters'
end subroutine
