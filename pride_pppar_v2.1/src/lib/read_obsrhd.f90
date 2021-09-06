!
!! read_obsrhd.f90
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
!! Contributor: Maorong Ge, Jianghui Geng, Songfeng Yang
!! 
!!
!!
!! purpose  : read rinex health diagnose information
!! parameter:
!!    input : jd,sod   -- requested time
!!            nprn,prn -- satellite
!!    output: OB       -- rinex observation struct
!
subroutine read_obsrhd(jd, sod, nprn, prn, OB)
  implicit none
  include '../header/const.h'
  include '../header/rnxobs.h'

  integer*4 jd, nprn
  character*3 prn(1:*)
  real*8 sod
  type(rnxobr) OB
!
!! local
  integer*4 i, iy, imon, id, ih, im, jdx, jde, ierr
  character*3 iprn
  real*8 dt, sec, sodx, sode, dintv
  character*256 line
!
!! function called
  integer*4 get_valid_unit, pointer_string, modified_julday
  real*8 timdif
!
!! header
  if (OB%rhdfil(1:1) .eq. ' ') return
  if (OB%rhdfil(1:1) .ne. ' ' .and. OB%lfnrhd .eq. 0) then
    OB%lfnrhd = get_valid_unit(10)
    open (unit=OB%lfnrhd, file=OB%rhdfil, status='OLD', iostat=ierr)
    if (ierr .ne. 0) then
      write (*, '(2a)') '***ERROR(read_obsrhd): open file ', trim(OB%rhdfil)
      call exit(1)
    endif
    line = ' '
    OB%amb_epo = 0
    OB%amb_tot = 0
    OB%ava_obs = 0
    OB%rem_obs = 0
    do while (index(line, 'END OF HEADER') .eq. 0)
      read (OB%lfnrhd, '(a)') line
      if (index(line, 'INTERVAL') .ne. 0) then
        read (line, *, err=100) dintv
      else if (index(line, 'AMB MAX/TOT/NEW') .ne. 0) then
        read (line, *, err=100) OB%amb_epo, OB%amb_tot
      else if (index(line, 'EPO AVA/REM/NEW') .ne. 0) then
        read (line, *, err=100) OB%ava_obs, OB%rem_obs
      else if (index(line, 'EFF EPO/SUM/NEW') .ne. 0) then
        read (line, *, err=100) OB%act_obs, OB%sum_obs
      endif
    enddo
    if (dintv .gt. sod) then
      write (*, '(a)') '***ERROR(read_obsrhd): rhd interval not small enough '
      call exit(1)
    endif
    if (OB%ava_obs .ne. 0 .and. OB%amb_epo .eq. 0) then
      write (*, '(a)') '***ERROR(read_obsrhd): MAX AMB not found'
      call exit(1)
    endif
    return
  endif
!
!! read epoch by epoch
  do while (.true.)
    read (OB%lfnrhd, '(a)', end=50) line
    if (line(1:3) .ne. 'TIM') cycle
    read (line, '(3x,i5,4i3,f11.7)') iy, imon, id, ih, im, sec
    jdx = modified_julday(id, imon, iy)
    sodx = ih*3600.d0 + im*60.d0 + sec
    dt = timdif(jdx, sodx, jd, sod)
    if (dt .gt. MAXWND) then
      backspace OB%lfnrhd
      exit
    else
      line = ' '
      do while (.true.)
        read (OB%lfnrhd, '(a)', end=50) line
        if (line(1:3) .eq. 'TIM') exit
        read (line, '(a3)') iprn
        i = pointer_string(nprn, prn, iprn)
        if (i .eq. 0) cycle
        if (line(61:63) .eq. 'DEL') then
          if (dt .lt. -MAXWND) cycle
          OB%obs(i, 1) = 0.d0
          OB%obs(i, 2) = 0.d0
          OB%obs(i, 3) = 0.d0
          OB%obs(i, 4) = 0.d0
        else if (line(61:63) .eq. 'AMB') then
          read (line, '(31x,i5,4i3,f11.7)') iy, imon, id, ih, im, sec
          jde = modified_julday(id, imon, iy)
          sode = ih*3600.d0 + im*60.d0 + sec
          if (timdif(jde, sode, jd, sod) .gt. -MAXWND) then
            OB%flag(i) = 1
            OB%lifamb(i, 1) = jdx + sodx/86400.d0
            OB%lifamb(i, 2) = jde + sode/86400.d0
          endif
        endif
      enddo
      backspace OB%lfnrhd
    endif
  enddo
50 backspace OB%lfnrhd
  return
100 write (*, '(a)') '***ERROR(read_obsrhd): read file ', trim(OB%rhdfil)
  call exit(1)
end
