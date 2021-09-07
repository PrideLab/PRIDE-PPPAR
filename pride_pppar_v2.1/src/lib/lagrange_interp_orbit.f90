!
!! lagrange_interp_orbit.f90
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
!! Contributor: Maorong Ge, Shuyin Mao
!! 
!!
!!
!! purpose  : lagrange interplolation of satellite position and partial derivatives
!!            based on table data
!!
!! parameter: iprn -- satellite PRN
!!            lcheck_only -- read table file head and check if the table is ok
!!                          for desired interpolation. Time, and return ref. time of ICS
!!            lpos  -- .true. for position
!!            lvel  -- .true. for velocity
!!            lpartial -- .true. for partial derivatives
!!            jd,sod  -- time of the requested point
!!            x,v,part  -- interpolated position, velocity and partial derivatives
!!
!
subroutine lagrange_interp_orbit(orbfil, lpos, lvel, jd, sod, iprn, x, v)
  implicit none
  include '../header/const.h'
  include '../header/orbit.h'

  logical*1 lpos, lvel
  integer*4 jd
  character*3 iprn
  real*8 sod, x(1:*), v(1:*)
  character*(*) orbfil
  type(orbhdr) OH
!
!! local
  integer*4, parameter :: MAXDGR = 10, MAXPNT = 2*(MAXDGR + 1)

  type orb_int_tab
    integer*4 lunit, ndgr
    integer*4 irec_inmemory, nrec_inmemory
    integer*4 jds, jde, nrec
    real*8 sods, sode
    real*8 dintv
  end type

  type orb_int_data
    real*8 table(MAXPNT, 3)
  end type
  type(orb_int_tab) LI
  type(orb_int_data) DT(MAXSAT)
!
!! local
  logical*1 first, update_alpha, update_memory
  integer*4 it, nprn
  character*3 prn(MAXSAT)
  integer*4 ierr, ipos, i, j, k, ivar, jpos, jepo, sepo, irec
  real*8 s, u, s2, u2, sum1, sum2, epo
  real*8 coeff(MAXPNT)
!
!! function
  integer*4 get_valid_unit, pointer_string
  real*8 timdif

  data first/.true./
  save first, LI, DT, nprn, prn
!
!! first time call
  if (first) then
    first = .false.
!
!! read header of the orbit file
    call rdorbh(orbfil, LI%lunit, OH)
    LI%jds = OH%jd0
    LI%sods = OH%sod0
    LI%dintv = OH%dintv
    LI%jde = OH%jd1
    LI%sode = OH%sod1
    LI%nrec = nint(timdif(OH%jd1, OH%sod1, OH%jd0, OH%sod0)/OH%dintv) + 1
    nprn = OH%nprn
    do i = 1, nprn
      prn(i) = OH%prn(i)
    enddo
    LI%ndgr = 8
    if (LI%ndgr .gt. MAXDGR) then
      write (*, '(a,2i3)') '***ERROR(lagrange_interp_orbit): LI%ndgr larger than MAXDGR', LI%ndgr, MAXDGR
      call exit(1)
    endif

    LI%irec_inmemory = 1
    LI%nrec_inmemory = LI%ndgr+1
    do irec = 1, LI%nrec_inmemory
      read (LI%lunit, end=300) k, ((DT(i)%table(irec, j), j=1, 3), i=1, nprn)
    enddo
  endif
!
!! check time
  if (timdif(jd, sod, LI%jds, LI%sods) .lt. -1.d-6 .or. timdif(jd, sod, LI%jde, LI%sode) .gt. 1.d-6) then
    !write (*, '(a,i5,f9.1)') '###WARNING(everett_interp_orbit): arc not cover epoch ', jd, sod
    !x(1:3) = 1.d15
    !v(1:3) = 1.d15
    !call exit(1)
  endif
!
!! which satellite
  it = pointer_string(nprn, prn, iprn)
  if (it .eq. 0) then
    write (*, '(a3)') '***ERROR(everett_interp_orbit): in table no data for satellite ', iprn
    call exit(1)
  endif
  
  !
  !! interpolation epoch
  epo=timdif(jd, sod, LI%jds, LI%sods)/LI%dintv+1.d0
  !
  !! start epoch of the interpolation arc
  sepo=nint(epo-LI%ndgr/2+1d-5/LI%dintv)
  
  if(sepo.lt.1) sepo=1
  if(sepo+LI%ndgr.gt.LI%nrec) sepo=LI%nrec-LI%ndgr
  
  !
  !! check if we have enough data in memory.
  !
  !! forwards
  if(LI%irec_inmemory + LI%nrec_inmemory - 1 .lt. sepo+LI%ndgr .and. &
     LI%irec_inmemory + LI%nrec_inmemory - 1 .lt. LI%nrec) then
    do while (LI%irec_inmemory .lt. sepo)
      do k = 1, nprn
        do j = 1, 3
          do i = 1, LI%nrec_inmemory-1
            DT(k)%table(i, j) = DT(k)%table(i + 1, j)
          enddo
        enddo
      enddo
      read (LI%lunit, end=300) k, ((DT(i)%table(LI%nrec_inmemory, j), j=1, 3), i=1, nprn)
      LI%irec_inmemory = LI%irec_inmemory + 1
    enddo
  !
  !! backwards
  else if(LI%irec_inmemory .gt. sepo) then
    do while (LI%irec_inmemory .gt. sepo)
      backspace LI%lunit
      LI%irec_inmemory = LI%irec_inmemory - 1
    enddo
    do i = 1, LI%nrec_inmemory
      backspace LI%lunit
    enddo
    do irec = 1, LI%nrec_inmemory
      read (LI%lunit, end=300) k, ((DT(i)%table(irec, j), j=1, 3), i=1, nprn)
    enddo
  endif  
  
  !
  !! if missing data
  if (any(DT(it)%table(1:LI%nrec_inmemory, 1) .eq. 1.d15)) then
    x(1:3) = 1.d15
    v(1:3) = 1.d15
    return
  endif
  !
  !!
  jepo=sepo-LI%irec_inmemory
  !
  !!
  call lagrange_coeff(LI%ndgr,epo-sepo,.true.,coeff)
  !
  !! interpolation
  do ivar = 1, 3
    x(ivar) = 0.d0
    v(ivar) = 0.d0
    do k = 1, LI%ndgr+1
      x(ivar)=x(ivar)+coeff(k)*DT(it)%table(jepo+k,ivar)
    enddo
    if (ivar .le. 3 .and. lvel) then
      do k = 1, LI%ndgr+1
      v(ivar) = v(ivar)+coeff(LI%ndgr+1+k)*DT(it)%table(jepo+k,ivar)
      enddo
      v(ivar)=v(ivar)/LI%dintv
    endif
  enddo
  return

300 continue
  write (*, '(a)') '***ERROR(everett_interpolate_orbit): end of file '
  call exit(1)
!
!! reset
  Entry lagrange_reset()
  close (LI%lunit)
  first = .true.
  return
end
