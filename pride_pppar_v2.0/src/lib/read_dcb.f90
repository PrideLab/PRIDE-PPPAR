!
!! read_dcb.f90
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
!! Contributor: Jianghui Geng, Songfeng Yang
!! 
!!
!!
!! purpose  : read dcb bias
!! parameter:
!! output   : dcb -- dcb bias
!
subroutine read_dcb(dcb)
  implicit none
  include '../header/const.h'

  integer*4 lfnp1, lfnp2
  real*8 dcb(MAXSAT, 2)
!
!! local variables
  integer*4 ierr
  character*3 prn0,prn_mat(MAXSAT)
  integer*4 i
  character*256 line
!
!! function used
  integer*4 get_valid_unit,pointer_string
  call prn_matbld(prn_mat)
!
!! dcb(P1C1)
  lfnp1 = get_valid_unit(10)
  open (lfnp1, file='P1C1.dcb', status='old', iostat=ierr)
  if (ierr .ne. 0) then
    write (*, '(2a)') '###WARNING(lsq): open file ', 'P1C1.dcb'
    goto 50
  endif
  i=0
  do while (.true.)
    read (lfnp1, '(a)', end=50) line
    if((line(1:1).eq.'G' .or. line(1:1).eq.'R') .and. line(2:3).ne.'  ') then
      read (line, '(a3)') prn0
      i=pointer_string(MAXSAT,prn_mat,prn0)
      read (line, '(a3,20x,f12.3)') prn0, dcb(i, 1)
    endif
  enddo
!
!! dcb(P2C2)
50 lfnp2 = get_valid_unit(10)
  open (lfnp2, file='P2C2.dcb', status='old', iostat=ierr)
  if (ierr .ne. 0) then
    write (*, '(2a)') '###WARNING(lsq): open file ', 'P2C2.dcb'
    goto 100
  endif
  i=0 
  do while (.true.)
    read (lfnp2, '(a)', end=100) line
    if((line(1:1).eq.'G' .or. line(1:1).eq.'R') .and. line(2:3).ne.'  ') then
      read (line, '(a3)') prn0
      i=pointer_string(MAXSAT,prn_mat,prn0)
      read (line, '(a3,20x,f12.3)') prn0, dcb(i, 2)
    endif
  enddo
100 close (lfnp1)
  close (lfnp2)
  return
end
