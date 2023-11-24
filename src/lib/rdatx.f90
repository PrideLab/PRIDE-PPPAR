!
!! rdatx.f90
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
!! purpose   : read antenna atx file
!! parameter :
!!    input  : fjd_beg,fjd_end -- time span
!!    output : ATX%-- antenna correction
!
subroutine rdatx(fjd_beg, fjd_end, ATX, snxcode)
  implicit none
  include '../header/const.h'
  include '../header/antatx.h'

  real*8        fjd_beg, fjd_end
  type(antatx)  ATX
  character*60  snxcode
!
!! local
  logical*1     lfirst, lfound
  integer*4     i0, i, j, k, lfn, ierr
  integer*4     iy, imon, id, ih, im
  integer*4     ncol, nrow, icol, irow
  real*8        sod, fjd0, fjd1
  character*1   atxtyp
  character*512 line
  character*5   lable
!
!! function called
  integer*4     get_valid_unit, modified_julday

  data lfirst/.true./
  save lfirst, lfn, atxtyp

  if (lfirst) then
    lfirst = .false.
    lfn = get_valid_unit(10)
    open (lfn, file='abs_igs.atx', status='old', iostat=ierr)
    if (ierr .ne. 0) then
      write (*, '(a,i4)') '***ERROR(rdatx): open antenna file abs_igs.atx', ierr
      call exit(1)
    end if
!
!! read header
    line = ' '
    do while (index(line, 'END OF HEADER') .ne. 61)
      read (lfn, '(a)', end=100) line
      if (index(line, 'PCV TYPE / REFANT') .eq. 61) then
        read (line, '(a1)') atxtyp
      end if
    end do
  end if
  if (atxtyp .ne. 'A') then
    write (*, '(a)') '***ERROR(rdatx): absolute antenna mode in abs_igs.atx'
    call exit(1)
  end if
!
!! find the right antenna
10 continue
  rewind lfn
  lfound = .false.
  do while (.true.)
    line = ' '
    do while (index(line, 'TYPE / SERIAL NO') .ne. 61)
      read (lfn, '(a)', end=100) line
    end do
    if (ATX%antnam(1:5) .eq. 'BLOCK'   .and. ATX%antnam(1:5) .ne. line(1:5) .or.  &
        ATX%antnam(1:7) .eq. 'GLONASS' .and. ATX%antnam(1:7) .ne. line(1:7) .or.  &
        ATX%antnam(1:7) .eq. 'GALILEO' .and. ATX%antnam(1:7) .ne. line(1:7) .or.  &
        ATX%antnam(1:6) .eq. 'BEIDOU'  .and. ATX%antnam(1:6) .ne. line(1:6) .or.  &
        ATX%antnam(1:4) .eq. 'QZSS'    .and. ATX%antnam(1:4) .ne. line(1:4) .or.  &
        ATX%antnam(1:5) .ne. 'BLOCK'   .and. ATX%antnam(1:7) .ne. 'GLONASS' .and. &
        ATX%antnam(1:7) .ne. 'GALILEO' .and. ATX%antnam(1:6) .ne. 'BEIDOU'  .and. ATX%antnam(1:4) .ne. 'QZSS' .and. &
        ATX%antnam .ne. line(1:20) .or. ATX%antnum .ne. line(21:40)) cycle
    ATX%antnam = line(1:20)
    ATX%dazi = 0.d0
    ATX%dzen = 0.d0
    ATX%nfreq = 0
    ATX%neu = 0.d0
    ATX%pcv(1:50, 0:200, 1:10, 1:5) = 0.d0
    fjd0 = 0.d0
    fjd1 = 1.d10
    do while (index(line, 'START OF FREQUENCY') .ne. 61)
      read (lfn, '(a)', end=100) line
      if (index(line, 'DAZI') .eq. 61) then
        read (line, *, err=200) ATX%dazi
      else if (index(line, 'ZEN1 / ZEN2 / DZEN') .eq. 61) then
        read (line, *, err=200) ATX%zen1, ATX%zen2, ATX%dzen
      else if (index(line, 'VALID FROM') .eq. 61) then
        read (line, *, err=200) iy, imon, id, ih, im, sod
        fjd0 = modified_julday(id, imon, iy) + ih/24.d0 + im/1440.d0 + sod/86400.d0
      else if (index(line, 'VALID UNTIL') .eq. 61) then
        read (line, *, err=200) iy, imon, id, ih, im, sod
        fjd1 = modified_julday(id, imon, iy) + ih/24.d0 + im/1440.d0 + sod/86400.d0
      else if (index(line, '# OF FREQUENCIES') .eq. 61) then
        read (line, *, err=200) ATX%nfreq
        if (ATX%nfreq .gt. 2) then
        end if
      else if (index(line, 'SINEX CODE') .eq. 61) then
        if (len_trim(snxcode) .gt. 0 .and. &
            snxcode(1:len_trim(snxcode)) .ne. line(1:len_trim(line(1:60)))) then
          write (*, '(3a)') '###WARNING(rdatx): inconsistent antex SINEX code: ', snxcode(1:len_trim(snxcode)), &
            line(1:len_trim(line(1:60)))
        else
          snxcode = line(1:len_trim(line(1:60)))
        end if
      end if
    end do
    if ((fjd0 - fjd_beg)*86400.d0 .gt. MAXWND .or. (fjd_end - fjd1)*86400.d0 .gt. MAXWND) cycle
    lfound = .true.
    exit
  end do
100 continue
  backspace lfn
  if (.not. lfound) then
    if (ATX%antnam(1:5) .ne. 'BLOCK'   .and. ATX%antnam(1:7) .ne. 'GLONASS' .and. &
        ATX%antnam(1:7) .ne. 'GALILEO' .and. ATX%antnam(1:6) .ne. 'BEIDOU'  .and. ATX%antnam(1:4) .ne. 'QZSS') then
      write (*, '(3a)') '###WARNING(rdatx): Antenna not found ', ATX%antnam, ATX%antnum
      if (ATX%antnam(17:20) .eq. 'NONE') then
        ATX%dazi = 5.d0
        ATX%zen1 = 0.d0
        ATX%zen2 = 90.d0
        ATX%dzen = 5.d0
        ATX%neu = 0.d0
        ATX%pcv(1:50, 0:200, 1:10, 1:5) = 0.d0
        write (*, '(a)') '$$$MESSAGE(rdatx): zero antenna applied instead'
      else
        ATX%antnam(17:20) = 'NONE'
        write (*, '(a)') '$$$MESSAGE(rdatx): dome NONE is tried. Search again ...'
        goto 10
      end if
      return
    else
      write (*, '(3a)') '###WARNING(rdatx): atx not found ', ATX%antnam, ATX%antnum
      return
    end if
  end if
!
!! read offset and phase center variation
  ncol = 0
  nrow = 0
  if (ATX%dazi .ne. 0.d0) nrow = int(360.d0/ATX%dazi) + 1
  if (ATX%dzen .ne. 0.d0) ncol = int((ATX%zen2 - ATX%zen1)/ATX%dzen) + 1
  ATX%frq = ''
  ATX%sys_multi = ''
  do i = 1, ATX%nfreq
!
!! find frequency
    line = ' '
    do while (line(61:78) .ne. 'START OF FREQUENCY')
      read (lfn, '(a)') line
      if (index(line, 'END OF ANTENNA') .ne. 0) then
        write (*, '(a)') '***ERROR(rdatx): frquency not found '
        call exit(1)
      end if
    end do
    i0 = index(GNSS_PRIO, line(4:4))
    if (i0 .eq. 0) cycle
    ATX%sys_multi(i0:i0) = line(4:4)
    read (line(5:6), *) k
    ATX%frq(k, i0) = line(4:6)
!
!! reading
    irow = 0
    do while (.true.)
      read (lfn, '(a)') line
      if (line(61:76) .eq. 'END OF FREQUENCY') exit
      if (line(61:77) .eq. 'NORTH / EAST / UP') then
        read (line, *, err=200) (ATX%neu(j, k, i0), j=1, 3)
!! convert unit from mm to m
        do j = 1, 3
          ATX%neu(j, k, i0) = ATX%neu(j, k, i0)*1.d-3
        end do
      else if (line(4:8) .eq. 'NOAZI') then
        read (line(9:), *, err=200) (ATX%pcv(icol, 0, k, i0), icol=1, ncol)
      else
        irow = irow + 1
        read (line(9:), *, err=200) (ATX%pcv(icol, irow, k, i0), icol=1, ncol)
      end if
    end do
    if (irow .ne. nrow) then
      write (*, '(a)') '***ERROR(rdatx): pcv lost '
      call exit(1)
    end if
!! convert unit from mm to m
    do irow = 0, nrow
      do icol = 1, ncol
        ATX%pcv(icol, irow, k, i0) = ATX%pcv(icol, irow, k, i0)*1.d-3
      end do
    end do
  end do

  return

200 continue
  write (*, '(2a)') '***ERROR(rdatx): read file ', trim(line)
  call exit(1)
end subroutine
