!
!! sp3orb.f90
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
!! Contributor: Maorong Ge, Jianghui Geng
!! 
!!
!
!
!! Orbit generator
!
program sp3orb
  implicit none
  include '../header/const.h'
  include '../header/orbit.h'

  character*128 :: sescfg = '', sp3fil = '', orbfil = '', erpfil = ''
  integer*4 lfnorb, jd
  real*8 x(6, MAXSAT), dt, sod, dif
!
!! for rotation matrix, only mat_trans is used
  real*8 :: mate2j(3, 3), rmte2j(3, 3), gast, xpole, ypole

  integer*4 i, j, k, ierr, iflag
  type(orbhdr) OH
!
!! function called
  integer*4 get_valid_unit
  real*8 timdif
!
!!initial
  do i = 1, MAXSAT
    do j = 1, 6
      x(j, i) = 0.d0
    enddo
  enddo
  do i = 1, 3
    do j = 1, 3
      mate2j(j, i) = 0.d0
      rmte2j(j, i) = 0.d0
    enddo
  enddo
!
!! get arguements
  call get_sp3orb_args(sescfg, sp3fil, orbfil, erpfil, OH)
  write (*, '(a,190a3)') ' GNSS Satellite: ', (OH%prn(i), i=1, OH%nprn)
!
!! open temp orbit file
  lfnorb = get_valid_unit(10)
  open (lfnorb, file=orbfil, form='unformatted')
!
!! write header of orbit file
  write (lfnorb) OH
!
!! read position from sp3 file
  k = 0
  jd = OH%jd0
  sod = OH%sod0
  do while (timdif(jd, sod, OH%jd1, OH%sod1) .lt. MAXWND)
    call rdsp3i(jd, sod, OH%nprn, OH%prn, x, iflag)
    k = k + 1
    if (iflag .eq. 1) then
      write (*, '(a,i6,f8.2)') '###WARNING(sp3orb): epoch lost ', jd, sod
    else
      call ef2int(erpfil, jd, sod, mate2j, rmte2j, gast, xpole, ypole)
      do i = 1, OH%nprn
        if (x(1, i) .eq. 1.d15) then
          write (*, '(a,i6,f9.2,a4)') '###WARNING(sp3orb): satellite lost ', jd, sod, OH%prn(i)
          cycle
        endif
        call matmpy(mate2j, x(1, i), x(1, i), 3, 3, 1)
      enddo
    endif
    write (lfnorb) k, ((x(j, i), j=1, 3), i=1, OH%nprn)
    call timinc(jd, sod, OH%dintv, jd, sod)
  enddo
  call rdsp3c()
  close (lfnorb)
end
