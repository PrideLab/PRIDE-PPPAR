!
!! read_mhm.f90
!!
!!    Copyright (C) 2023 by Wuhan University
!!
!!    This program belongs to PRIDE PPP-AR which is an open source software:
!!    you can redistribute it and/or modify it under the terms of the GNU
!!    General Public License (version 3) as published by the Free Software
!Foundation.
!!
!!    This program is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!!    GNU General Public License (version 3) for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with this program. If not, see <https://www.gnu.org/licenses/>.
!!
!! Contributor: Honghai Zhang
!!
!!
!!
!! purpose   : read the mhm model file
!! parameters:
!!             mhm          -- MHM model
!!             sit_name     -- station name
!
subroutine read_mhm(mhm,sit_name)
  implicit none
  include '../header/const.h'
  include '../header/station.h'
  
  integer :: i, j, k, ios
  integer :: gps_s = 0, gps_e = 0, gal_s = 0, gal_e = 0, bds_s = 0, bds_e = 0, gln_s = 0, gln_e = 0, qzs_s = 0, qzs_e = 0
  integer :: col1, col2, line_num
  real :: col3, col7, col4, col5, col6
  character(len=200) :: flnmhm, line, modelname
  logical :: end_of_file 
  real*8 :: mhm(360, 90, 2, 5)
  integer*4 :: get_valid_unit, lfn
  character*4 sit_name

  mhm=0.0d0
  end_of_file = .false.
!
  lfn = get_valid_unit(600)
  modelname='./mhm_'
  flnmhm = trim(modelname) // trim(sit_name)
  write (*, *)  trim(flnmhm) 
!
  open(lfn, file=trim(flnmhm), status='old', action='read', iostat=ios)
  if (ios /= 0) then
    write (*, '(2a)') '###WARNING(read_mhm): open file ', trim(flnmhm)
    end_of_file = .true.
	return
  endif

  line_num=0
  do while (.not. end_of_file)
    read(lfn, '(A)', iostat=ios) line
    
	if (ios /= 0) then
        end_of_file = .true.
        exit
    end if
	
    line_num = line_num + 1
    if (index(line, "START OF") > 0 .and. index(line, "GPS") > 0 ) then
      gps_s = line_num
    end if
    if (index(line, "END OF") > 0 .and. index(line, "GPS") > 0 ) then
      gps_e = line_num
    end if    
    if (index(line, "START OF") > 0 .and. index(line, "GAL") > 0 ) then
      gal_s = line_num
    end if 
    if (index(line, "END OF") > 0 .and. index(line, "GAL") > 0 ) then
      gal_e = line_num
    end if 
    if (index(line, "START OF") > 0 .and. index(line, "BDS") > 0 ) then
      bds_s = line_num
    end if 
    if (index(line, "END OF") > 0 .and. index(line, "BDS") > 0 ) then
      bds_e = line_num
    end if 
    if (index(line, "START OF") > 0 .and. index(line, "GLO") > 0 ) then
      gln_s = line_num
    end if 
    if (index(line, "END OF") > 0 .and. index(line, "GLO") > 0 ) then
      gln_e = line_num
    end if
  
  end do 
  
  rewind(lfn)
  
  if(gps_s > 0) then
    do i = 1, gps_s + 1
      read(lfn, '(A)', iostat=ios) line
      
    end do

    do i = gps_s + 2, gps_e - 1 
      read(lfn, '(A)', iostat=ios) line
     
      read(line, *) col1, col2, col3, col4, col5, col6, col7
      mhm(col1+1, col2+1, 1, 1) = col3
      mhm(col1+1, col2+1, 2, 1) = col7
    end do
  end if 

  rewind(lfn)

  if(gal_s > 0) then
    do i = 1, gal_s + 1
      read(lfn, '(A)', iostat=ios) line
    end do

    do i = gal_s + 2, gal_e - 1
      read(lfn, '(A)', iostat=ios) line
      read(line, *) col1, col2, col3, col4, col5, col6, col7
      mhm(col1+1, col2+1, 1, 3) = col3
      mhm(col1+1, col2+1, 2, 3) = col7
    end do
  end if

  rewind(lfn)

  if(bds_s > 0) then
    do i = 1, bds_s + 1
      read(lfn, '(A)', iostat=ios) line
    end do

    do i = bds_s + 2, bds_e - 1
      read(lfn, '(A)', iostat=ios) line
      read(line, *) col1, col2, col3, col4, col5, col6, col7
      mhm(col1+1, col2+1, 1, 4) = col3
      mhm(col1+1, col2+1, 2, 4) = col7
    end do
  end if

  rewind(lfn)

  if(gln_s > 0) then
    do i = 1, gln_s + 1
      read(lfn, '(A)', iostat=ios) line
    end do

    do i = gln_s + 2, gln_e - 1
      read(lfn, '(A)', iostat=ios) line
      read(line, *) col1, col2, col3, col4, col5, col6, col7
      mhm(col1+1, col2+1, 1, 2) = col3
      mhm(col1+1, col2+1, 2, 2) = col7
    end do
  end if

  rewind(lfn)

  if(qzs_s > 0) then
    do i = 1, qzs_s + 1
      read(lfn, '(A)', iostat=ios) line
    end do

    do i = qzs_s + 2, qzs_e - 1
      read(lfn, '(A)', iostat=ios) line
      read(line, *) col1, col2, col3, col4, col5, col6, col7
      mhm(col1+1, col2+1, 1, 5) = col3
      mhm(col1+1, col2+1, 2, 5) = col7
    end do
  end if

  close(lfn)
  return
end subroutine
