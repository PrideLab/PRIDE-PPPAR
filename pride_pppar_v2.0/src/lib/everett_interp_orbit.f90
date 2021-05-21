!
!! everett_interp_orbit.f90
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
!! purpose  : everett interplolation of satellite position and partial derivatives
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
subroutine everett_interp_orbit(orbfil, lpos, lvel, jd, sod, iprn, x, v)
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
  integer*4, parameter :: MAXDGR = 7, MAXPNT = 2*(MAXDGR + 1)

  type orb_int_tab
    integer*4 lunit, ndgr
    integer*4 irec_inmemory, nrec_inmemory, irec_middle_alpha
    integer*4 jds, jde, nrec
    real*8 sods, sode
    real*8 dintv
    real*8 ec(0:MAXDGR, 0:MAXDGR)
  end type

  type orb_int_data
    real*8 table(MAXPNT, 3)
    real*8 alpha(3, 0:MAXDGR)
    real*8 beta(3, 0:MAXDGR)
  end type
  type(orb_int_tab) EI
  type(orb_int_data) DT(MAXSAT)
!
!! local
  logical*1 first, update_alpha, update_memory
  integer*4 it, nprn
  character*3 prn(MAXSAT)
  integer*4 ierr, ipos, i, j, k, ivar, jpos
  real*8 s, u, s2, u2, sum1, sum2
!
!! function
  integer*4 get_valid_unit, pointer_string
  real*8 timdif

  data first/.true./
  save first, EI, DT, nprn, prn
!
!! first time call
  if (first) then
    first = .false.
!
!! read header of the orbit file
    call rdorbh(orbfil, EI%lunit, OH)
    EI%jds = OH%jd0
    EI%sods = OH%sod0
    EI%dintv = OH%dintv
    EI%jde = OH%jd1
    EI%sode = OH%sod1
    EI%nrec = nint(timdif(OH%jd1, OH%sod1, OH%jd0, OH%sod0)/OH%dintv) + 1
    nprn = OH%nprn
    do i = 1, nprn
      prn(i) = OH%prn(i)
    enddo
    EI%ndgr = 6
    if (EI%ndgr .gt. MAXDGR) then
      write (*, '(a,2i3)') '***ERROR(everett_interp_orbit): EI%ndgr larger than MAXDGR', EI%ndgr, MAXDGR
      call exit(1)
    endif

    EI%irec_inmemory = 0
    EI%nrec_inmemory = 0
    EI%irec_middle_alpha = 0
    call everett_coeff(MAXDGR, EI%ndgr, EI%ec)
  endif
!
!! check time
  if (timdif(jd, sod, EI%jds, EI%sods) .lt. -1.d-6 .or. timdif(jd, sod, EI%jde, EI%sode) .gt. 1.d-6) then
    write (*, '(a,i5,f9.1)') '***ERROR(everett_interp_orbit): arc not cover epoch ', jd, sod
    x(1:3) = 1.d15
    v(1:3) = 1.d15
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
!! interpolation center
  ipos = int((jd - EI%jds)*(86400.d0/EI%dintv) + (sod - EI%sods)/EI%dintv) + 1
  if (ipos .le. EI%ndgr) then   ! Geng. 9/29/2006 'lt'-->'le'
    ipos = EI%ndgr + 1           ! Geng. 7/02/2006 interpolation at beginning is allowed
  else if (ipos + EI%ndgr + 1 .gt. EI%nrec) then
    ipos = EI%irec_middle_alpha
  endif
!
!! check if we have to update alpha and beta
  update_alpha = ipos .gt. EI%irec_middle_alpha
!
!! check if we have enough data in memory.
!! for interpolation, 2(n+1) points are needed. In case of a multi-points table,
!! at least, two lines should be in memory
  update_memory = update_alpha
  do while (update_memory)
    update_memory = EI%nrec_inmemory .lt. 2*EI%ndgr + 2 .or. &
                    ipos + EI%ndgr + 1 .gt. EI%irec_inmemory + EI%nrec_inmemory
    if (.not. update_memory) exit
    if (EI%nrec_inmemory .ge. 2*EI%ndgr + 2) then
      EI%nrec_inmemory = EI%nrec_inmemory - 1
      EI%irec_inmemory = EI%irec_inmemory + 1
      do k = 1, nprn
        do j = 1, 3
          do i = 1, EI%nrec_inmemory
            DT(k)%table(i, j) = DT(k)%table(i + 1, j)
          enddo
        enddo
      enddo
    endif
    read (EI%lunit, end=300) k, ((DT(i)%table(EI%nrec_inmemory + 1, j), j=1, 3), i=1, nprn)
    EI%nrec_inmemory = EI%nrec_inmemory + 1
  enddo
!
!! update alpha and beta
  if (update_alpha) then
    EI%irec_middle_alpha = ipos
    jpos = ipos - EI%irec_inmemory
    do i = 1, nprn
!! if missing, alpha and beta then unavailable
      if (any(DT(i)%table(1:EI%nrec_inmemory, 1) .eq. 1.d15)) then
        DT(i)%alpha(1, 0:EI%ndgr) = 0.d0
        cycle
      endif
      do ivar = 1, 3
        do k = 0, EI%ndgr
          DT(i)%alpha(ivar, k) = 0.d0
          DT(i)%beta(ivar, k) = 0.d0
          do j = 0, EI%ndgr
            DT(i)%alpha(ivar, k) = DT(i)%alpha(ivar, k) + EI%ec(k, j)* &
                                   (DT(i)%table(jpos - j, ivar) + DT(i)%table(jpos + j, ivar))
            DT(i)%beta(ivar, k) = DT(i)%beta(ivar, k) + EI%ec(k, j)* &
                                  (DT(i)%table(jpos - j + 1, ivar) + DT(i)%table(jpos + j + 1, ivar))
          enddo
        enddo
      enddo
    enddo
  endif
!
!! if missing ..., send a note to upper-level routine
  if (all(DT(it)%alpha(1, 0:EI%ndgr) .eq. 0.d0)) then
    x(1:3) = 1.d15
    v(1:3) = 1.d15
    return
  endif
!
!! interpolation
  s = (jd - EI%jds)*(86400.d0/EI%dintv) + (sod - EI%sods)/EI%dintv + 1.d0 - EI%irec_middle_alpha
  u = s - 1.d0
  s2 = s*s
  u2 = u*u

  do ivar = 1, 3
    sum1 = DT(it)%beta(ivar, EI%ndgr)
    sum2 = DT(it)%alpha(ivar, EI%ndgr)
    do k = EI%ndgr, 1, -1
      sum1 = sum1*s2 + DT(it)%beta(ivar, k - 1)
      sum2 = sum2*u2 + DT(it)%alpha(ivar, k - 1)
    enddo
    sum1 = sum1*s - sum2*u
    x(ivar) = sum1
    if (ivar .le. 3 .and. lvel) then
      sum1 = DT(it)%beta(ivar, EI%ndgr)*(2.d0*EI%ndgr + 1.d0)
      sum2 = DT(it)%alpha(ivar, EI%ndgr)*(2.d0*EI%ndgr + 1.d0)
      do k = EI%ndgr, 1, -1
        sum1 = sum1*s2 + DT(it)%beta(ivar, k - 1)*(2.d0*(k - 1.d0) + 1.d0)
        sum2 = sum2*u2 + DT(it)%alpha(ivar, k - 1)*(2.d0*(k - 1.d0) + 1.d0)
      enddo
      v(ivar) = (sum1 - sum2)/EI%dintv
    endif
  enddo
  return

300 continue
  write (*, '(a)') '***ERROR(everett_interpolate_orbit): end of file '
  call exit(1)
!
!! reset
  Entry everett_reset()
  close (EI%lunit)
  first = .true.
  return
end
