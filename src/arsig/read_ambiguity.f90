!
!! read_ambiguity.f90
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
!! Contributor: Maorong Ge, Jianghui Geng, Songfeng Yang
!! 
!!
!!
!! purpose  : read un-differenced ambiguity estimates
!! parameter:
!!    input : ACF -- fractional cycle biases
!!    output: AS  -- ambiguity estimates at a station
!
subroutine read_ambiguity(ACF, AS)
  implicit none
  include '../header/const.h'
  include 'arscfg.h'
  include 'ambsit.h'

  type(arscfg) ACF
  type(ambsit) AS
!
!! local
  integer*4 i, j, lfn, ierr, isat
  character*3 iprn, prn(MAXSAT+20)
  real*8 t(2), elev, xamb, xrwl, xrms, xswl
  character*256 line
!
!! function called
  integer*4 get_valid_unit, pointer_string
!
!! open ambiguity file
  lfn = get_valid_unit(10)
  open (lfn, file=ACF%flnamb, status='old', iostat=ierr)
  if (ierr .ne. 0) then
    write (*, '(2a)') '***ERROR(read_ambiguity): open file ', trim(ACF%flnamb)
    call exit(1)
  endif
!
!! read un-differenced ambiguity estimates
!
!! read header
  isat = 1
  do i = 1, MAXSAT
    ACF%fcbprn(i) = ''
  enddo
  prn = ''
  ACF%fcbnprn = 0
  do while (.true.)
    read (lfn, '(a)', end=100) line
    if (index(line, 'STATION') .ne. 0) then
      read (line, '(a4)') AS%name
    elseif (index(line, '# OF AMB RESOLVABLE SAT') .ne. 0) then
      read (line, '(i5)') ACF%fcbnprn
    elseif (index(line, 'AMB RESOLVABLE SATELLITES') .ne. 0) then
      read (line, '(15(a3,1x))') prn(isat:isat + 14)
      do i = isat, isat + 15
        if (prn(i) .eq. '') then
          exit
        else
          ACF%fcbprn(i) = prn(i)
        endif
      enddo
      isat = i
    elseif (index(line, 'END OF HEADER') .ne. 0) then
      exit
    endif
  enddo
  AS%now = 0
  ACF%nobs = 0
  do while (.true.)
    read (lfn, '(a)', iostat=ierr) line
    if (ierr .ne. 0) exit
    if (line(1:1) .eq. '*') cycle
    read (line, '(a4,2f22.6,2f18.10,2f9.4,f6.1)') iprn, xamb, xrwl, t(1:2), xrms, xswl, elev
    if (elev .le. ACF%cutoff .or. xswl*3.d0 .gt. 0.2d0 .or. (t(2) - t(1))*86400.d0 .lt. ACF%minsec_common) cycle
    if (pointer_string(ACF%nprn, ACF%prn, iprn) .eq. 0) then
      write (*, '(a,a3)') '***ERROR(read_ambiguity): satellite not exist ', iprn
      call exit(1)
    endif
!
!! add a new eligible ambiguity
    AS%now = AS%now + 1
    if (AS%now .gt. MAXOW_ST) then
      write (*, '(a)') '***ERROR(read_ambiguity): too many un-differenced ambiguity estimates'
      call exit(1)
    endif
    AS%isat(AS%now) = pointer_string(ACF%nprn, ACF%prn, iprn)
    AS%xamb(AS%now) = xamb
    AS%xrwl(AS%now) = xrwl
    AS%xrms(AS%now) = xrms
    AS%xswl(AS%now) = xswl
    AS%iepc(1, AS%now) = nint(((t(1) - ACF%jd0)*86400.d0 - ACF%sod0)/ACF%dintv) + 1
    AS%iepc(2, AS%now) = nint(((t(2) - ACF%jd0)*86400.d0 - ACF%sod0)/ACF%dintv) + 1
    ACF%nobs = AS%iepc(2, AS%now) - AS%iepc(1, AS%now) + 1 + ACF%nobs
  enddo
  close (lfn)

  return
100 write (*, '(2a)') '***ERROR(read_ambiguity): end of file ', trim(ACF%flnamb)
  call exit(1)
end
