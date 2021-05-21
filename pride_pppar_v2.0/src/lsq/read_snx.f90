!
!! read_fcb.f90
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
!! Contributor: Jianghui Geng, Songfeng Yang, Xingyu Chen
!! 
!!
!!
!! purpose  : read fractional-cycle biases for ambiguity resolution in PPP
!! parameter:
!!    input : FCB -- fractional part of uncalibrated phase delays
!
subroutine read_snx(flnfcb, bias)
  implicit none
  include '../header/const.h'

!
!! local
  integer*4 i, j, lfn, jd0, jd1, iy, imon, id, ih, im, ierr
  character*3 iprn, jprn
  character*3 prn_mat(MAXSAT)
  integer*4 ii
  real*8 sod0, sod1, is, nfcb, wfcb, l1(32), l2(32), p1(32), p2(32)
  real*8 bias(MAXSAT, 4), lamdw, lamdn, g, g1
  character*256 line
  character*20 flnfcb
!
!! function used
  integer*4 get_valid_unit, modified_julday, pointer_string
  real*8 timdif
  real*8 f1, f2, a1, a2, bc, dw, dn, aw1, aw2, an1, an2

  call prn_matbld(prn_mat)
!
!! open file
  lfn = get_valid_unit(500)
  open (lfn, file=flnfcb, status='old', iostat=ierr)
  if (ierr .ne. 0) then
    write (*, '(2a)') '***WARNING(read_fcb): open file ', trim(flnfcb)
    bias=0.d0
    return
  endif
!
!! read header
  read (lfn, '(a)') line
  do while (index(line, '%=ENDBIA') .eq. 0)
    if (index(line, '+BIAS/SOLUTION') .ne. 0) then
      !
      !! read fractional parts per satellite
      do while (.true.)
        read (lfn, '(a)', iostat=ierr) line
        if (ierr .ne. 0) exit
        if (index(line, '-BIAS/SOLUTION') .ne. 0) goto 300
        if (index(line, 'OSB') .eq. 0) cycle
        !! read bias
        read (line, '(11x,a3)', err=200) iprn
        if (iprn .eq. '') cycle
        ii=pointer_string(MAXSAT,prn_mat,iprn)
        if (ii .eq. 0) cycle
        if (line(16:19) .ne. '') cycle
        if (iprn(1:1) .eq. 'G') then
          if (index(line, 'L1C') .ne. 0) then
            read (line, '(70x,f21.5)') bias(ii, 1)
            bias(ii, 1) = bias(ii, 1)*VLIGHT*1d-9
          elseif (index(line, 'L2W') .ne. 0) then
            read (line, '(70x,f21.5)') bias(ii, 2)
            bias(ii, 2) = bias(ii, 2)*VLIGHT*1d-9
          elseif (index(line, 'C1W') .ne. 0) then
            read (line, '(70x,f21.4)') bias(ii, 3)
            bias(ii, 3) = bias(ii, 3)*VLIGHT*1d-9
          elseif (index(line, 'C2W') .ne. 0) then
            read (line, '(70x,f21.4)') bias(ii, 4)
            bias(ii, 4) = bias(ii, 4)*VLIGHT*1d-9
          endif
        elseif (iprn(1:1) .eq. 'E') then
          if (index(line, 'L1C') .ne. 0) then
            read (line, '(70x,f21.5)') bias(ii, 1)
            bias(ii, 1) = bias(ii, 1)*VLIGHT*1d-9
          elseif (index(line, 'L5Q') .ne. 0) then
            read (line, '(70x,f21.5)') bias(ii, 2)
            bias(ii, 2) = bias(ii, 2)*VLIGHT*1d-9
          elseif (index(line, 'C1C') .ne. 0) then
            read (line, '(70x,f21.4)') bias(ii, 3)
            bias(ii, 3) = bias(ii, 3)*VLIGHT*1d-9
          elseif (index(line, 'C5Q') .ne. 0) then
            read (line, '(70x,f21.4)') bias(ii, 4)
            bias(ii, 4) = bias(ii, 4)*VLIGHT*1d-9
          endif
        elseif (iprn(1:1) .eq. 'C') then
          if (index(line, 'L2I') .ne. 0) then
            read (line, '(70x,f21.5)') bias(ii, 1)
            bias(ii, 1) = bias(ii, 1)*VLIGHT*1d-9
          elseif (index(line, 'L6I') .ne. 0) then
            read (line, '(70x,f21.5)') bias(ii, 2)
            bias(ii, 2) = bias(ii, 2)*VLIGHT*1d-9
          elseif (index(line, 'C2I') .ne. 0) then
            read (line, '(70x,f21.4)') bias(ii, 3)
            bias(ii, 3) = bias(ii, 3)*VLIGHT*1d-9
          elseif (index(line, 'C6I') .ne. 0) then
            read (line, '(70x,f21.4)') bias(ii, 4)
            bias(ii, 4) = bias(ii, 4)*VLIGHT*1d-9
          endif
        elseif (iprn(1:1) .eq. 'J') then
          if (index(line, 'L1C') .ne. 0) then
            read (line, '(70x,f21.5)') bias(ii, 1)
            bias(ii, 1) = bias(ii, 1)*VLIGHT*1d-9
          elseif (index(line, 'L2L') .ne. 0) then
            read (line, '(70x,f21.5)') bias(ii, 2)
            bias(ii, 2) = bias(ii, 2)*VLIGHT*1d-9
          elseif (index(line, 'C1C') .ne. 0) then
            read (line, '(70x,f21.4)') bias(ii, 3)
            bias(ii, 3) = bias(ii, 3)*VLIGHT*1d-9
          elseif (index(line, 'C2L') .ne. 0) then
            read (line, '(70x,f21.4)') bias(ii, 4)
            bias(ii, 4) = bias(ii, 4)*VLIGHT*1d-9
          endif
        endif
      enddo
    endif
    read (lfn, '(a)', end=100) line
  enddo
300 close (lfn)

  return
100 write (*, '(2a)') '***ERROR(read_fcb): end of file ', trim(flnfcb)
  call exit(1)
200 write (*, '(2a)') '***ERROR(read_fcb): read file ', trim(flnfcb)
  call exit(1)
end
