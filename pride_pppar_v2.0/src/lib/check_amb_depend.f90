!
!! check_amb_depend.f90
!!
!!    Copyright (C) 2021 by Wuhan University
!!
!!    This program is free software: you can redistribute it and/or modify
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
!! purpose  : check if a new ambiguity is dependant of the already selected
!!            set expressed by a set of orth. unit vector Ei.
!!
!!     With a set of orth. unit vector of the oneway amb space, Ei i=1,namb,
!!     any candidate can be expressed by a = sum (ci * Ei). Assume we already
!!     defined the first ndef Ei, the new vector a is dependent on Ei i=1,..ndef,
!!     if for any i>ndef, ci=0.
!!     Take dot on both side of the expression, we have
!!        a dot a = sum(ci * Ei) dot sum(ci * Ei)
!!                = sum (ci**2)
!!     For testing, if <a,a>= sum (ci**2), i=1,ndef, a is dependent.
!!     Otherwise, the to be added E(ndef+1) is unit(a - sum(ci * Ei))
!!
!! parameter:
!!            namb     -- # of total one way amb
!!            ndef     -- # of candidates already selected
!!            pt2ow    -- position of the 2 oneway amb
!!            depend   -- .true. dependant
!!
!
subroutine check_amb_depend(namb, ndef, pt2ow, depend)
  implicit none
  include '../header/const.h'

  logical*1 depend
  integer*4 namb, ndef, pt2ow(1:*)
!
!! local
  logical*1 lfirst
  integer*4 i, j, ierr
  integer*4, pointer :: idx(:)
  real*8, pointer :: e(:), c(:)  ! orth. unit vectors represent the selected candidates
  real*8 operator(2), c_dot, EPS

  data lfirst, EPS/.true., 1.d-12/, operator/1.d0, -1.d0/
  save lfirst, e, c, idx, EPS
!
!! initialization
  if (lfirst) then
    lfirst = .false.
    nullify (e); nullify (c); nullify (idx)
  endif
!
!! memory allocation. namb denotes total number of one-way ambiguities
  if (ndef .eq. 0) then  !! reallocation
    if (associated(e)) deallocate (e)
    if (associated(c)) deallocate (c)
    if (associated(idx)) deallocate (idx)
    if (namb .gt. 0) then
      allocate (e(namb*namb), stat=ierr)
      if (ierr .ne. 0) then
        write (*, '(a,i10)') '***ERROR(check_amb_depend): memory allocation for e ', namb
        call exit(1)
      endif
      allocate (c(namb), stat=ierr)
      if (ierr .ne. 0) then
        write (*, '(a)') '***ERROR(check_amb_depend): memory allocation for c '
        call exit(1)
      endif
      allocate (idx(namb), stat=ierr)
      if (ierr .ne. 0) then
        write (*, '(a)') '***ERROR(check_amb_depend): memory allocation for idx '
        call exit(1)
      endif
      do i = 1, namb
        idx(i) = (i - 1)*namb
      enddo
    else !! only deallocate memory
      return
    endif
  endif
!
!! compute Ci. Ci is saved for adding new e(ndef+1) if it is independent
  c_dot = 0.d0
  do i = 1, ndef
    c(i) = 0.d0
    do j = 1, 2
      c(i) = c(i) + e(idx(i) + pt2ow(j))*operator(j)
    enddo
    c_dot = c_dot + c(i)*c(i)
  enddo
!
!! 2 == dot a
  depend = .false.
  if (dabs(c_dot - 2) .lt. EPS) then
    depend = .true.
    return
  endif
!
!! a new independent vector
  if (ndef + 1 .gt. namb) then
    write (*, '(a)') '***ERROR(check_amb_depend): redundant independent vector '
    call exit(1)
  endif
  e(idx(ndef + 1) + 1:idx(ndef + 1) + namb) = 0.d0
  do j = 1, 2
    e(idx(ndef + 1) + pt2ow(j)) = operator(j)
  enddo
!
!! a - sum (ci * Ei)
  if (c_dot .gt. EPS) then
    do i = 1, ndef
      if (dabs(c(i)) .lt. EPS) cycle
      do j = 1, namb
        e(idx(ndef + 1) + j) = e(idx(ndef + 1) + j) - c(i)*e(idx(i) + j)
      enddo
    enddo
!
!! normalize new indenpendent vector
    c_dot = 0.d0
    do j = 1, namb
      c_dot = c_dot + e(idx(ndef + 1) + j)**2
    enddo
    c_dot = dsqrt(c_dot)
    do j = 1, namb
      e(idx(ndef + 1) + j) = e(idx(ndef + 1) + j)/c_dot
    enddo
  else
!! first vector
    c_dot = dsqrt(2*1.d0)
    do j = 1, 2
      e(idx(ndef + 1) + pt2ow(j)) = e(idx(ndef + 1) + pt2ow(j))/c_dot
    enddo
  endif

  ndef = ndef + 1
  return
end
