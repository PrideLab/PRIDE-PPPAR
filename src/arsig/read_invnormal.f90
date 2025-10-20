!
!! read_invnormal.f90
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
!! purpose  : read inversed normal matrix (Q matrix)
!! parameter:
!!    input : ACF -- fractional part of initial phases
!!    output: QN  -- inversed normal matrix
!!            AS  -- ambiguity station struct
!
subroutine read_invnormal(ACF, PM, QN, AS)
  implicit none
  include '../header/const.h'
  include 'invnor.h'
  include 'arscfg.h'
  include 'ambsit.h'

! common
  integer*4     idxfrq(MAXSYS, 2)
  common        idxfrq
! parameter
  type(arscfg)  ACF
  type(ambsit)  AS
  type(pest)    PM(1:*)
  type(invm)    QN
! local
  integer*4     lfn, i0, i, j, isat, nprn, ierr
  character*3   prn(MAXSAT)
  real*8        xrwl, xswl, elev
  real*8        f1(MAXSYS)
  character*20  anttyp
! function called
  integer*4     get_valid_unit
  integer*4     pointer_string

!
!! frequency
  do i0 = 1, MAXSYS
    f1(i0) = FREQ_SYS(idxfrq(i0, 1), i0)
    if (f1(i0) .eq. 0.d0) goto 100
  end do
!
!! open file
  lfn = get_valid_unit(10)
  open (lfn, file=ACF%flnneq, form='unformatted', status='old', iostat=ierr)
  if (ierr .ne. 0) then
    write (*, '(2a)') '***ERROR(read_invnormal): open file ', trim(ACF%flnneq)
    call exit(1)
  end if
  AS%now = 0
  ACF%nobs = 0
  read (lfn) AS%name
  read (lfn) ACF%fcbnprn, (ACF%fcbprn(i), i=1, ACF%fcbnprn)
  read (lfn) nprn, (prn(i), i=1, nprn)
  read (lfn) QN%ntot, QN%vtpv, QN%frdm
  if (QN%ntot .gt. 3 + MAXOW_ST) then
    write (*, '(a)') '***ERROR(read_invnormal): too many parameters '
    call exit(1)
  end if
  allocate (QN%invx(1:QN%ntot, 1:QN%ntot), stat=ierr)
  if (ierr .ne. 0) then
    write (*, '(a)') '***ERROR(read_invnormal): memory allocation invx '
    call exit(1)
  end if
  allocate (QN%idq(QN%ntot), stat=ierr)
  do i = 1, QN%ntot
    QN%idq(i) = (i - 1)*QN%ntot
  end do
!
!! read parameters
  QN%nxyz = 0
  do i = 1, QN%ntot
    read (lfn) PM(i)%pname, isat, xrwl, PM(i)%xest, PM(i)%ptime(1:2), xswl, elev
    PM(i)%psat = isat
    if (PM(i)%pname(1:4) .ne. 'AMBC') then
      QN%nxyz = QN%nxyz + 1
      cycle
    end if
    PM(i)%psat = pointer_string(ACF%nprn, ACF%prn, prn(isat))
    if (PM(i)%psat .eq. 0) then
      write (*, '(a,a3)') '***ERROR(read_invnormal): satellite not exist ', prn(isat)
      call exit(1)
    end if
    if (elev .le. ACF%cutoff .or. xswl*3.d0 .gt. 0.2d0 .or. (PM(i)%ptime(2) - PM(i)%ptime(1))*86400.d0 .lt. ACF%minsec_common) cycle
    AS%now = AS%now + 1
    if (AS%now .gt. MAXOW_ST) then
      write (*, '(a)') '***ERROR(read_invnormal): too many one-way ambiguities'
      call exit(1)
    end if
    AS%isat(AS%now) = PM(i)%psat
    i0 = index(GNSS_PRIO, ACF%prn(PM(i)%psat)(1:1))
    PM(i)%xest = PM(i)%xest*f1(i0)/VLIGHT
    AS%xamb(AS%now) = PM(i)%xest
    AS%xrwl(AS%now) = xrwl
    AS%xswl(AS%now) = xswl
    AS%iepc(1, AS%now) = nint(((PM(i)%ptime(1) - ACF%jd0)*86400.d0 - ACF%sod0)/ACF%dintv) + 1
    AS%iepc(2, AS%now) = nint(((PM(i)%ptime(2) - ACF%jd0)*86400.d0 - ACF%sod0)/ACF%dintv) + 1
    ACF%nobs = AS%iepc(2, AS%now) - AS%iepc(1, AS%now) + 1 + ACF%nobs
    AS%ipt(AS%now) = i       ! pointer to Q matrix
  end do
!
!! read inversed normal matrix
  QN%invx(1:QN%ntot, 1:QN%ntot) = 0.d0
  read (lfn) ((QN%invx(i, j), i=j, QN%ntot), j=1, QN%ntot)   ! lower part
  close (lfn)
  return

100 continue
  write (*, '(a,5(1x,a,2i1))') '***ERROR(lsq_add_obs): invalid frequency number:', &
    'G', idxfrq(index(GNSS_PRIO, 'G'), 1:2), &
    'R', idxfrq(index(GNSS_PRIO, 'R'), 1:2), &
    'E', idxfrq(index(GNSS_PRIO, 'E'), 1:2), &
    'C', idxfrq(index(GNSS_PRIO, 'C'), 1:2), &
    'J', idxfrq(index(GNSS_PRIO, 'J'), 1:2)
  call exit(1)
end subroutine
