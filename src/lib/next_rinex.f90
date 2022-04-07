!
!! next_rinex.f90
!!
!!    Copyright (C) 2022 by Wuhan University
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
!! Contributor: Jihang Lin
!! 
!!
!!
!! purpose  : read one epoch data from a RINEX o-file
!!
!! parameter: lfn_this -- input file unit
!!            lfn_next -- successor file unit of input lfn,
!!                        default as 0 if no valid successor file
!!
!
subroutine next_rinex(lfn_this, lfn_next)
  implicit none

  integer*4 lfn_this
  integer*4 lfn_next
!
!! local
  character*256 rinex_dir
  character*256 name_this
  character*256 name_next
  character*1   rinex_ver
  integer*4 nlen
  integer*4 doy, year, mjd
  integer*4 ind1, ind2
  integer*4 ierr
  logical*4 is_alive
!
!! function used
  integer*4  get_valid_unit
  integer*4  modified_julday

  lfn_next = 0
  inquire(unit=lfn_this, name=name_this, exist=is_alive) 
  if (name_this .eq. '') then
    write (*, '(2a)') '###WARNING(next_rinex): open ', trim(name_this) 
    return
  endif
!
!! split directory and base
  ind1=index(name_this, '/',  BACK = .TRUE.)
  ind2=index(name_this, '\\', BACK = .TRUE.)
  if (ind1 .ne. 0) then
    rinex_dir = trim(name_this(1:ind1))
    name_this = trim(name_this(ind1+1:))
  elseif (ind2 .ne. 0) then
    rinex_dir = trim(name_this(1:ind2+1))
    name_this = trim(name_this(ind2+2:))
  endif
!
!! try writing successor file name
  nlen = len_trim(name_this)
  select case (nlen)
    case (12)
      rinex_ver = '2'
      if (name_this(8:9) .ne. "0.") then
        goto 110
      elseif (index('ongplm', name_this(12:12)) .ne. 0) then
        read (name_this(1:nlen), '(4x,i3,2x,i2,1x)', err=110) doy, year
        if (year < 80) then
          year = year + 2000
        else
          year = year + 1900
        endif
        mjd = modified_julday(doy, 0, year)
        call mjd2doy(mjd+1, year, doy)
        write (name_next, '(a4,i0.3,a2,i0.2,a1)') &
          name_this(1:4), doy, name_this(8:9), mod(year,100), name_this(12:12)
      else
        goto 120
      endif
    case (34,38)
      rinex_ver = '3'
      if (name_this(nlen-3:nlen) .ne. ".rnx") then
        goto 110
      elseif (name_this(20:28) .ne. "0000_01D_") then
        goto 110
      elseif (index('ONM', name_this(nlen-4:nlen-4)) .ne. 0) then
        read (name_this(1:nlen), '(12x,i4,i3,19x)', err=110) year, doy
        mjd = modified_julday(doy, 0, year)
        call mjd2doy(mjd+1, year, doy)
        write (name_next, '(a12,i0.4,i0.3,a)') &
          name_this(1:12), year, doy, name_this(20:nlen)
      else
        goto 120
      endif
    case default
      write (*, '(2a)') '###WARNING(next_rinex): unrecognized naming of RINEX files: ', trim(name_this)
      return
  endselect

  name_next = trim(rinex_dir)//trim(name_next)
  
  lfn_next = get_valid_unit(10)
  open(lfn_next, file=name_next, status='old', iostat=ierr)
  if (ierr .ne. 0) then
    write(*,'(2a)') '###WARNING(next_rinex): open ', trim(name_next)
    lfn_next = 0
    return
  else
    write(*,'(2a)') '%%%MESSAGE(next_rinex): open ', trim(name_next)
  endif

  return

110 continue
  write (*,'(2a)') '###WARNING(next_rinex): illegal naming of RINEX-'//rinex_ver//' files: ', trim(name_this)
  return

120 continue
  write (*,'(2a)') '###WARNING(next_rinex): unrecognized RINEX-'//rinex_ver//' file type: ', name_this(nlen-3:nlen)
  return
end

