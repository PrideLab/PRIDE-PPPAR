!
!! next_rinex.f90
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
!! Contributor: Jihang Lin, Jing Zeng
!! 
!!
!!
!! purpose  : read one epoch data from a RINEX o-file
!!
!! parameter: lfn_this -- input file unit
!!            lfn_next -- successor file unit of input lfn,
!!                        default as 0 if no valid successor file
!!            mjd_end  -- loop to match end date of RINEX file
!!
!
subroutine next_rinex(lfn_this, lfn_next, mjd_end)
  implicit none

  integer*4 lfn_this
  integer*4 lfn_next
  integer*4 mjd_end 
!
!! local
  character*256 rinex_dir
  character*256 name_this
  character*256 name_next
  character*1   rinex_ver
  integer*4 nlen
  integer*4 doy, year, mjd1, mjd2
  integer*4 ind1, ind2
  integer*4 ierr
  logical*1 is_alive
!
!! function used
  integer*4  get_valid_unit
  integer*4  modified_julday

  lfn_next = 0
  name_next = ''
  inquire(unit=lfn_this, name=name_this, exist=is_alive) 
100 continue
  if (name_this .eq. '') return
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
  end if
!
!! try writing successor file name
  nlen = len_trim(name_this)
  select case (nlen)
    case (12)
      rinex_ver = '2'
      if (name_this(9:9) .ne. ".") goto 110
      if (index('ongplm', name_this(12:12)) .eq. 0) goto 120
      read (name_this(1:nlen), '(4x,i3,2x,i2,1x)', err=110) doy, year
      if (year < 80) then
        year = year + 2000
      else
        year = year + 1900
      end if
      mjd1 = modified_julday(doy, 0, year)
      if (name_this(8:8) .eq. "0") then
        mjd2 = mjd1 + 1
        call mjd2doy(mjd2, year, doy)
        write (name_next, '(a4,i0.3,a2,i0.2,a1)') &
          name_this(1:4), doy, name_this(8:9), mod(year,100), name_this(12:12)
      elseif (name_this(8:8) .ge. 'a' .and. name_this(8:8) .le. 'x') then
        ind1 = iachar(name_this(8:8)) - 97
        ind2 = ind1 + 1
        mjd2 = mjd1 + ind2/24
        ind2 = ind2 - (mjd2-mjd1)*24
        call mjd2doy(mjd2, year, doy)
        write (name_next, '(a4,i0.3,a1,a1,i0.2,a1)') &
          name_this(1:4), doy, achar(ind2+97), '.', mod(year,100), name_this(12:12)
      end if
    case (34,38)
      rinex_ver = '3'
      if (name_this(nlen-3:nlen) .ne. ".rnx" .and. name_this(nlen-3:nlen) .ne. ".RNX") goto 110
      if (index('ONM', name_this(nlen-4:nlen-4)) .eq. 0) goto 120
      read (name_this(1:nlen), '(12x,i4,i3,i2)', err=110) year, doy, ind1
      mjd1 = modified_julday(doy, 0, year)
      if (index(name_this(1:nlen), '_01D_') .ne. 0) then
        mjd2 = mjd1 + 1
        call mjd2doy(mjd2, year, doy)
        write (name_next, '(a12,i0.4,i0.3,a)') &
          name_this(1:12), year, doy, name_this(20:nlen)
      elseif (index(name_this(1:nlen), '_01H_') .ne. 0) then
        ind2 = ind1 + 1
        mjd2 = mjd1 + ind2/24
        ind2 = ind2 - (mjd2-mjd1)*24
        call mjd2doy(mjd2, year, doy)
        write (name_next, '(a12,i0.4,i0.3,i0.2,a)') &
          name_this(1:12), year, doy, ind2, name_this(22:nlen)
      end if
    case default
      write (*, '(2a)') '###WARNING(next_rinex): unrecognized naming of RINEX files: ', trim(name_this)
      return
  end select
  if (mjd2 .gt. mjd_end) return

  name_next = trim(rinex_dir)//trim(name_next)
  
  lfn_next = get_valid_unit(10)
  open (lfn_next, file=name_next, status='old', iostat=ierr)
  if (ierr .ne. 0) then
    write (*,'(2a)') '###WARNING(next_rinex): open ', trim(name_next)
    lfn_next = 0
    name_this = name_next
    goto 100
  else
    write (*,'(2a)') '%%%MESSAGE(next_rinex): open ', trim(name_next)
  end if

  return

110 continue
  write (*,'(2a)') '###WARNING(next_rinex): illegal naming of RINEX-'//rinex_ver//' files: ', trim(name_this)
  return

120 continue
  write (*,'(2a)') '###WARNING(next_rinex): unrecognized RINEX-'//rinex_ver//' file type: ', name_this(nlen-3:nlen)
  return
end subroutine
