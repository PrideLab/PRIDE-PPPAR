!
!! oceanload_coef.f90
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
!! purpose  : get oceanload coefficients
!! parameter:
!!    input : sname   -- station name
!!            lat,lon -- geodetic position of station
!!            rlat,rlon -- reference position corresponding to olc
!!    output: olc     -- oceanload coefficients
!
subroutine oceanload_coef(name, lat, lon, rlat, rlon, olc)
  implicit none
  include '../header/const.h'

  character*4 name
  real*8 olc(11, 6), lat, lon, rlat, rlon
!
!! local
  logical*1 lfirst
  integer*4 i, j, lfn, ierr
  character*256 line
  real*8 lat_deg, lon_deg, deg2rad, tol_d, tol_lat, tol_lon, dlat, dlon
!
!! function called
  integer*4 get_valid_unit

  data lfirst/.true./
  save lfirst, lfn, deg2rad, tol_lat
!
!! open oceanload file
  if (lfirst) then
    lfirst = .false.
    lfn = get_valid_unit(10)
    open (lfn, file='oceanload', status='old', iostat=ierr)
    if (ierr .ne. 0) then
      write (*, '(a)') '***ERROR(oceanload_coef): open file oceanload '
      call exit(1)
    endif
    deg2rad = datan(1.d0)/45.d0
    tol_d = 1.d4
    tol_lat = tol_d/ERAD/1.d3
  endif
  tol_lon = tol_lat/dcos(lat)
!
!! whether new coef should be looked for
  if (rlat .ne. 0.d0 .and. rlon .ne. 0.d0) then
    dlon = dabs(rlon - lon)
    dlat = dabs(rlat - lat)
    if (dlat .le. tol_lat .and. dlon .le. tol_lon) return
  endif
!
!! look for a new one
  do j = 1, 6
    do i = 1, 11
      olc(i, j) = 0.d0
    enddo
  enddo

  rewind lfn
  do while (.true.)
    read (lfn, '(a)', end=100) line
    i = index(trim(line), 'lon/lat:')
    if (i .ne. 0) then
      read (line(i + 8:), *) lon_deg, lat_deg
      if (lon_deg .lt. 0.d0) lon_deg = lon_deg + 360.d0
      dlon = dabs(lon_deg*deg2rad - lon)
      dlat = dabs(lat_deg*deg2rad - lat)
      if (dlat .le. tol_lat .and. dlon .le. tol_lon) then
        rlon = lon_deg*deg2rad       ! saved for next time
        rlat = lat_deg*deg2rad
        do j = 1, 6
          read (lfn, *, iostat=ierr) (olc(i, j), i=1, 11)
          if (ierr .ne. 0) then
            write (*, '(a)') '***ERROR(oceanload_coef): read file oceanload '
            call exit(1)
          endif
! convert to radian
          if (j .ge. 4) then
            do i = 1, 11
              olc(i, j) = olc(i, j)*deg2rad
            enddo
          endif
        enddo
        return
      endif
    endif
  enddo

100 write (*, '(2a)') '###WARNING(oceanload_coef): no oceanload coefficients for ', name
  return
end
