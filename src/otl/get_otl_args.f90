!
!! get_otl_args.f90
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
!! Contributor: Yingda Deng
!
subroutine get_otl_args(OTL)
  implicit none
  include 'otl.h'
  type(otls),intent(out) :: OTL ! ocean thermal loading settings
!
!! local
  integer*4 nargs, i, n, ymdhm(5), stat
  real*8 sec, mjd0, mjd
  character*256 :: line
  character*100 :: buffer(6)
!
!! functions called
  integer*4 get_valid_unit
!
!! read command arguments
  nargs = iargc()
  if (nargs .lt. 14) then
    write (*, '(a)') "                                                            "
    write (*, '(a)') "otl version 3.0, based on Scientific Computation Platform   "
    write (*, '(a)') "  for Geophysical Geodesy. Wuhan University, Oct. 2024      "
    write (*, '(a)') "                                                            "
    write (*, '(a)') "Usage: otl -b lat -l lon -h hgt -s ts -e te -i ti -o output "
    write (*, '(a)') "                                                            "
    write (*, '(a)') "Description:                                                "
    write (*, '(a)') "  otl is a module of PRIDE PPP-AR, used for calculating the "
    write (*, '(a)') "  tide correction component of the station.                 "
    write (*, '(a)') "                                                            "
    write (*, '(a)') "Required arguments:                                         "
    write (*, '(a)') "  -b lat                                                    "
    write (*, '(a)') "    station latitude. (Unit: degree)                        "
    write (*, '(a)') "  -l lon                                                    "
    write (*, '(a)') "    station longitude. (Unit: degree)                       "
    write (*, '(a)') "  -h hgt                                                    "
    write (*, '(a)') "    station ellipsoid height. (Unit: meter)                 "
    write (*, '(a)') "  -s ts                                                     "
    write (*, '(a)') "    start time. ts = ""year/month/day hour:minute:second""  "
    write (*, '(a)') "  -e te                                                     "
    write (*, '(a)') "    end time. te = ""year/month/day hour:minute:second""    "
    write (*, '(a)') "  -i ti                                                     "
    write (*, '(a)') "    sample interval. (Unit: second)                         "
    write (*, '(a)') "  -o output                                                 "
    write (*, '(a)') "    output file name.                                       "
    write (*, '(a)') "                                                            "
    write (*, '(a)') "Note: Some dependent files need to be under folders.        "
    write (*, '(a)') "                                                            "
    write (*, '(a)') "Examples:                                                   "
    write (*, '(a,$)') "  otl -b 50.7980638 -l 4.3585607 -h 159.6967 -s 2024/10/1 "
    write (*, '(a)') "00:00:00.00 -e 2024/10/1 23:59:59.00 -i 30 -o otl_brux      "
    write (*, '(a)') "                                                            "
    write (*, '(a)') "More details refer to PRIDE PPP-AR manual and repository    "
    write (*, '(a)') "  https://github.com/PrideLab/PRIDE-PPPAR/                  "
    write (*, '(a)') "                                                            "
    call exit(1)
  endif
  i = 1
  do while (i .le. nargs)
    call getarg(i, line)
    i = i + 1
    if (line(1:2) .eq. '-b') then
      call getarg(i, line)
      read (line, *) OTL%plat
    elseif (line(1:2) .eq. '-l') then
      call getarg(i, line)
      read (line, *) OTL%plon
    elseif (line(1:2) .eq. '-h') then
      call getarg(i, line)
      read (line, *) OTL%phgt
    elseif (line(1:2) .eq. '-s') then
      call getarg(i, line)
      call split_string(.false., line, ' ', ' ', '/', n, buffer)
      read(buffer(1:3),*) ymdhm(1), ymdhm(2),ymdhm(3)
      i=i+1
      call getarg(i, line)
      call split_string(.false., line, ' ', ' ', ':', n, buffer)
      read(buffer(1:3),*) ymdhm(4), ymdhm(5),sec
      call CAL2JD(ymdhm(1), ymdhm(2), ymdhm(3), mjd0, mjd, stat)
      OTL%mjd0=INT(mjd)
      OTL%sod0=ymdhm(4)*3600+ymdhm(5)*60+sec
    elseif (line(1:2) .eq. '-e') then
      call getarg(i, line)
      call split_string(.false., line, ' ', ' ', '/', n, buffer)
      read(buffer(1:3),*) ymdhm(1), ymdhm(2),ymdhm(3)
      i=i+1
      call getarg(i, line)
      call split_string(.false., line, ' ', ' ', ':', n, buffer)
      read(buffer(1:3),*) ymdhm(4), ymdhm(5),sec
      call CAL2JD(ymdhm(1), ymdhm(2), ymdhm(3), mjd0, mjd, stat)
      OTL%mjd1=INT(mjd)
      OTL%sod1=ymdhm(4)*3600+ymdhm(5)*60+sec
    elseif (line(1:2) .eq. '-i') then
      call getarg(i, line)
      read (line, *) OTL%interval    
    elseif (line(1:2) .eq. '-o') then
      call getarg(i, line)
      OTL%otlname=line
    endif
    i = i + 1
  enddo

end

