!
!! sp3orb.f90
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
!! Contributor: Maorong Ge, Jianghui Geng, Jihang Lin
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
  logical*1 lrepeat
  type(orbhdr) OH
!
!! function called
  integer*4 get_valid_unit
  real*8 timdif
!
!!initial
  x = 0.d0
  mate2j = 0.d0
  rmte2j = 0.d0
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
  k = 1
  jd = OH%jd0
  sod = OH%sod0
  do while (timdif(jd, sod, OH%jd1, OH%sod1) .lt. MAXWND)
    lrepeat = .false.
101 continue
    call rdsp3i(jd, sod, OH%nprn, OH%prn, x, iflag)
    if (iflag .eq. 1) then
      if (lrepeat) goto 102
      write (*, '(a,i6,f8.2)') '###WARNING(sp3orb): epoch lost ', jd, sod
    else
      call ef2int(erpfil, jd, sod, mate2j, rmte2j, gast, xpole, ypole)
      do i = 1, OH%nprn
        if (x(1, i) .eq. 1.d15) then
          write (*, '(a,i6,f9.2,a4)') '###WARNING(sp3orb): satellite lost ', jd, sod, OH%prn(i)
          cycle
        end if
        call matmpy(mate2j, x(1, i), x(1, i), 3, 3, 1)
      end do
    end if
    write (lfnorb) k, lrepeat, ((x(j, i), j=1, 3), i=1, OH%nprn)
    if (.not. (jd .eq. OH%jd0 .and. sod .eq. OH%sod0) .and. &
        .not. (jd .eq. OH%jd1 .and. sod .eq. OH%sod1) .and. &
        .not. lrepeat .and. sod .eq. 0.d0) then
      lrepeat = .true.
      goto 101
    end if
102 continue
    call timinc(jd, sod, OH%dintv, jd, sod)
    k = k + 1
  end do
  call rdsp3c()
  close (lfnorb)
end program
