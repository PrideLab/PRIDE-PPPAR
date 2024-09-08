!
!!oceanload.f90
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
program oceanload
  implicit none
  include 'otl.h'
  include '../header/const.h'
  integer*4 :: fmjd, mjds_utc, mjde_utc  ! find mjd, lfn_otl
  real*8 :: fsod, sods_utc, sode_utc     ! find sod
  integer*4 :: lfnotl, lfnotl_30min
  real*8 :: enu(3)      ! result of de dn du
  type(otls) :: OTL     ! ocean thermal loading settings
  integer*4, parameter :: DIM = 10
  real*8, parameter :: INT = 1800.d0
  ! function called
  real*8, external :: timdif, taiutc
  integer*4, external :: get_valid_unit
  !
  !! Get arguments from command line
  call get_otl_args(OTL)
  !
  !! Convert mjd to utc time
  call timinc(OTL%mjd0,OTL%sod0,19.d0-taiutc(OTL%mjd0),mjds_utc,sods_utc)
  call timinc(OTL%mjd1,OTL%sod1,19.d0-taiutc(OTL%mjd1),mjde_utc,sode_utc)
  !
  !! Generate reference results
  lfnotl_30min = get_valid_unit(10)
  open(lfnotl_30min,file="tmp_otl_30min",status="replace")
  call Otdloalharmsynth(OTL%plon,OTL%plat,OTL%phgt,mjds_utc,sods_utc,mjde_utc,sode_utc,INT,DIM,lfnotl_30min)
  !
  !! Do lagrange interpolation
  lfnotl = get_valid_unit(10)
  open(lfnotl,file=OTL%otlname,access='stream',form='unformatted',status="replace")
  fmjd=mjds_utc
  fsod=sods_utc
  do while (timdif(OTL%mjd1, OTL%sod1, fmjd, fsod) .gt. -MAXWND)
    rewind(lfnotl_30min)
    call lag_interp_otl(OTL, fmjd, fsod, DIM, lfnotl_30min, INT, enu)
    write(lfnotl) fmjd,fsod,enu(1:3)
    call timinc(fmjd, fsod, OTL%interval, fmjd, fsod)
  enddo
  close(lfnotl_30min,status='delete')
  close(lfnotl)

end program
