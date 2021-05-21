!
!! read_ambiguity.f90
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
!! Contributor: Maorong Ge, Jianghui Geng, Songfeng Yang
!! 
!!
!!
!! purpose  : read un-differenced ambiguity estimates
!! parameter:
!!    input : FCB -- fractional cycle biases
!!    output: AS  -- ambiguity estimates at a station
!
subroutine read_ambiguity(FCB, AS)
  implicit none
  include '../header/const.h'
  include 'abfcb.h'
  include 'ambsit.h'

  type(abfcb) FCB
  type(ambsit) AS
!
!! local
  integer*4 i, j, lfn, ierr
  character*3 iprn
  real*8 t(2), elev, xamb, xrwl, xrms, xswl
  character*4 name
  character*256 line
!
!! function called
  integer*4 get_valid_unit, pointer_string
!
!! open ambiguity file
  lfn = get_valid_unit(10)
  open (lfn, file=FCB%flnamb, status='old', iostat=ierr)
  if (ierr .ne. 0) then
    write (*, '(2a)') '***ERROR(read_ambiguity): open file ', trim(FCB%flnamb)
    call exit(1)
  endif
!
!! read un-differenced ambiguity estimates
  AS%now = 0
  do while (.true.)
    read (lfn, '(a)', iostat=ierr) line
    if (ierr .ne. 0) exit
    read (line, '(a4,a4,2f22.6,2f18.10,2f9.4,f6.1)') name, iprn, xamb, xrwl, t(1:2), xrms, xswl, elev
    if (elev .le. FCB%cutoff .or. xswl*3.d0 .gt. 0.2d0 .or. (t(2) - t(1))*86400.d0 .lt. FCB%minsec_common) cycle
    if (pointer_string(FCB%nprn, FCB%prn, iprn) .eq. 0) then
      write (*, '(a,a3)') '***ERROR(read_ambiguity): satellite not exist ', iprn
      call exit(1)
    endif
!
!! check whether a new station
    AS%name = name
!
!! add a new eligible ambiguity
    AS%now = AS%now + 1
    if (AS%now .gt. MAXOW_ST) then
      write (*, '(2a)') '***ERROR(read_ambiguity): too many un-differenced ambiguity estimates ', AS%name
      call exit(1)
    endif
    AS%isat(AS%now) = pointer_string(FCB%nprn, FCB%prn, iprn)
    AS%xamb(AS%now) = xamb
    AS%xrwl(AS%now) = xrwl
    AS%xrms(AS%now) = xrms
    AS%xswl(AS%now) = xswl
    AS%iepc(1, AS%now) = nint(((t(1) - FCB%jd0)*86400.d0 - FCB%sod0)/FCB%dintv) + 1
    AS%iepc(2, AS%now) = nint(((t(2) - FCB%jd0)*86400.d0 - FCB%sod0)/FCB%dintv) + 1
  enddo
  close (lfn)

  return
end
