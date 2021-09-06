!
!! enucov.f90
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
!! Contributor: Yuanxin Pan
!!
!!
program enucov
  implicit none
!! local
  logical*1 lfirst
  integer*4 lfnpos, lfncov, ierr, nargs, nobs
  real*8 xyz(3), blh(3), xyzref(3), dx(3), enu(3), sig0, mjd
  real*8 xyzcov(3,3), enucov1(3,3), rot(3,3), F(3,3)
  character*4   site, site0
  character*256 posfile, buf
!! function called
  integer*4 get_valid_unit

  nargs = iargc()
  if (nargs .lt. 1) then
    write (*, '(a)') 'Usage: enucov pos_file (per site)'
    call exit(1)
  endif
  call getarg(1, posfile)

  lfnpos = get_valid_unit(10)
  open(lfnpos, file=posfile, status='old', iostat=ierr)
  if (ierr .ne. 0) then
    write(*,'(2a)') 'error: no such file: ', trim(posfile)
    call exit(1)
  endif

  !! output enu & enucov
  lfncov = get_valid_unit(10)
  open(lfncov, file='enucov.out', status='replace', iostat=ierr)

  lfirst = .true.
  xyzref = 0.d0
  !! read xyz & xyzcov
  do while (.true.)
    read(lfnpos, '(a)', iostat=ierr) buf
    if (ierr .ne. 0) exit
    if (len_trim(buf) .eq. 0) cycle
    read(buf, *) site, mjd, xyz, xyzcov(1,1),xyzcov(2,2),xyzcov(3,3), &
         xyzcov(1,2),xyzcov(1,3),xyzcov(2,3), sig0, nobs

    if (lfirst) then
      lfirst = .false.
      site0 = site
      xyzref = xyz
      !! compute blh & enucov
      call xyzblh(xyz, 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, blh)
      call rot_enu2xyz(blh(1),blh(2),rot)
      F = transpose(rot)
    endif

    if (site .ne. site0) cycle

    dx = xyz - xyzref
    call matmpy(F, dx, enu, 3,3,1)

    xyzcov(2,1) = xyzcov(1,2)
    xyzcov(3,1) = xyzcov(1,3)
    xyzcov(3,2) = xyzcov(2,3)
    call matmpy(F, xyzcov, enucov1, 3, 3, 3)
    call matmpy(enucov1, rot, enucov1, 3, 3, 3)
    enucov1 = sig0*sig0*enucov1;
    
    write(lfncov, '(a4,f12.4,3f15.4,6e25.14,i15)') site, mjd, enu, enucov1(1,1), enucov1(2,2), &
          enucov1(3,3), enucov1(1,2), enucov1(1,3), enucov1(2,3), nobs
  enddo

  close(lfnpos)
  close(lfncov)
end program
