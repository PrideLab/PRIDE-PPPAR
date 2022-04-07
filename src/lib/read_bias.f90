!
!! read_bias.f90
!!
!!    Copyright (C) 2022 by Wuhan University
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
!! Contributor: Jianghui Geng, Songfeng Yang, Xingyu Chen, Jihang Lin
!! 
!!
!!
!! purpose  : read absolute code/phase biases for ambiguity resolution in PPP
!! parameter:
!!    input : OSB -- absolute code/phase biases
!
subroutine read_bias(flnosb, nprn, prn, bias, osbjd0, osbjd1)
  implicit none
  include '../header/const.h'
!
!! local
  integer*4 nprn, iprn, ntyp, ityp, isft
  integer*4 lfn, jd0, jd1
  integer*4 iy, imon, id, ih, im, ierr
  integer*4 iyear0, idoy0, isod0, imon0, id0
  integer*4 iyear1, idoy1, isod1, imon1, id1
  real*8 sod0, sod1
  real*8 osbjd0, osbjd1
  real*8 val, bias(MAXSAT, 36)
  character*2 styp(0:3)
  character*3 cprn, prn(1:*)
  character*3 ctyp
  character*20 flnosb
  character*256 line
!
!! function used
  integer*4 get_valid_unit, modified_julday, pointer_string
!
!! open file
  lfn = get_valid_unit(500)
  open (lfn, file=flnosb, status='old', iostat=ierr)
  if (ierr .ne. 0) then
    write (*, '(2a)') '###WARNING(read_bias): open file ', trim(flnosb)
    bias=0.d0
    return
  endif
!
!! read header
  read (lfn, '(a)') line
  do while (.true.)
    if (index(line, '+BIAS/SOLUTION') .ne. 0) then
      !
      !! read fractional parts per satellite
      do while (.true.)
 50     read (lfn, '(a)', iostat=ierr) line
        if (ierr .ne. 0) then
          backspace lfn
          exit
        endif
        if (index(line, 'OSB') .eq. 0) cycle
        if (line(16:19) .ne. '') cycle
        !
        !! read bias
        read (line, '(2(11x,a3),7x,2(i4,1x,i3,1x,i5,1x))', err=200) &
            cprn, ctyp, iyear0, idoy0, isod0, iyear1, idoy1, isod1
        if (cprn .eq. '') cycle
        iprn = pointer_string(nprn, prn, cprn)
        if (iprn .eq.  0) cycle
        call yeardoy2monthday(iyear0, idoy0, imon0, id0)
        jd0 = modified_julday(id0, imon0, iyear0)
        call yeardoy2monthday(iyear1, idoy1, imon1, id1)
        jd1 = modified_julday(id1, imon1, iyear1)
        if (jd1 + isod1/86400.d0 .le. osbjd0 .or. jd0 + isod0/86400.d0 .ge. osbjd1) cycle
        read (line, '(70x,f21.15)') val
        !
        !! get index of bias type
        if (cprn(1:1) .eq. 'G') then
          ityp = index(obs_prio_G,ctyp(3:3))
          styp = (/'L1', 'L2', 'C1', 'C2'/)
        elseif (cprn(1:1) .eq. 'R') then
          !! read only code bias
          ityp = index(obs_prio_R,ctyp(3:3))
          styp = (/'XX', 'XX', 'C1', 'C2'/)
        elseif (cprn(1:1) .eq. 'E') then
          ityp = index(obs_prio_E,ctyp(3:3))
          styp = (/'L1', 'L5', 'C1', 'C5'/)
        elseif (cprn(1:1) .eq. 'C') then
          ityp = index(obs_prio_C,ctyp(3:3))
          styp = (/'L2', 'L6', 'C2', 'C6'/)
        elseif (cprn(1:1) .eq. 'J') then
          ityp = index(obs_prio_J,ctyp(3:3))
          styp = (/'L1', 'L2', 'C1', 'C2'/)
        endif
        if (ityp .eq. 0) cycle
        do isft = 0, 3
          if (index(ctyp, styp(isft)) .ne. 0) exit
          if (isft .eq. 3) goto 50
        enddo
        !
        !! conversion from nanosecond to meter
        bias(iprn, ityp+9*isft) = val*VLIGHT*1d-9
      enddo
    endif
    read (lfn, '(a)', end=100) line
  enddo

100 continue
  if (index(line, '%=ENDBIA') .eq. 0) then 
    write (*, '(2a)') '***ERROR(read_bias): end of file ', trim(flnosb)
    call exit(1)
  else
    close (lfn)
  endif
!
!! complete phase biases
  do iprn = 1, nprn
    cprn = prn(iprn)
    if (cprn .eq. '') then
      cycle
    elseif (cprn(1:1) .eq. 'G') then
      ntyp = len_trim(obs_prio_G)
    elseif (cprn(1:1) .eq. 'R') then
      ntyp = len_trim(obs_prio_R)
    elseif (cprn(1:1) .eq. 'E') then
      ntyp = len_trim(obs_prio_E)
    elseif (cprn(1:1) .eq. 'C') then
      ntyp = len_trim(obs_prio_C)
    elseif (cprn(1:1) .eq. 'J') then
      ntyp = len_trim(obs_prio_J)
    endif
    do isft = 0, 1
      !
      !! assumes all phase biases are equal then find substitute
      val = 0.d0
      do ityp = 1, ntyp
        val = bias(iprn, ityp+9*isft)
        if (val .ne. 0.d0) exit
      enddo
      !
      !! complete but pass over when code biases not exists
      do ityp = 1, ntyp
        if (bias(iprn, ityp+9*isft) .eq. 0.d0) then
          if (bias(iprn, ityp+9*(isft+2)) .eq. 0.d0) cycle
          bias(iprn, ityp+9*isft) = val
        endif
      enddo
    enddo
  enddo

  return

200 continue
  write (*, '(2a)') '***ERROR(read_bias): read file ', trim(flnosb)
  call exit(1)

end subroutine
