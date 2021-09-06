!
!! read_bias.f90
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
!! Contributor: Jianghui Geng, Songfeng Yang, Xingyu Chen
!! 
!!
!!
!! purpose  : read fractional-cycle biases for ambiguity resolution in PPP
!! parameter:
!!    input : FCB -- fractional part of uncalibrated phase delays
!
subroutine read_bias(flnfcb, bias, fcbjd0, fcbjd1)
  implicit none
  include '../header/const.h'
!
!! local
  integer*4 i, j, lfn, jd0, jd1, iy, imon, id, ih, im, ierr, index_obs
  character*3 iprn, jprn
  character*3 prn_mat(MAXSAT)
  integer*4 ii, iprn_int
  real*8 sod0, sod1, is, nfcb, wfcb, l1(32), l2(32), p1(32), p2(32)
  integer*4 iyear0, idoy0, isod0, iy0, imon0, id0, ih0, imin0, is0, iyear1, idoy1, isod1, iy1, imon1, id1, ih1, imin1, is1
  real*8 bias(MAXSAT, 44), lamdw, lamdn, g, g1
  character*256 line
  character*20 flnfcb
  character*12 obs_prio_G,obs_prio_E,obs_prio_C,obs_prio_J,obs_type
  real*8 fcbjd0, fcbjd1
!
!! function used
  integer*4 get_valid_unit, modified_julday, pointer_string
  real*8 timdif
  real*8 f1, f2, a1, a2, bc, dw, dn, aw1, aw2, an1, an2
  obs_prio_G = 'NMYXLSCWP  '
  obs_prio_E = 'XCBQI'
  obs_prio_C = 'XQI'
  obs_prio_J = 'ZCXLS  '
  call prn_matbld(prn_mat)
!
!! open file
  lfn = get_valid_unit(500)
  open (lfn, file=flnfcb, status='old', iostat=ierr)
  if (ierr .ne. 0) then
    write (*, '(2a)') '###WARNING(read_fcb): open file ', trim(flnfcb)
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
        read (lfn, '(a)', iostat=ierr) line
        if (ierr .ne. 0) then
          backspace lfn
          exit
        endif
        if (index(line, 'OSB') .eq. 0) cycle
        !! read bias
        read (line, '(11x,a3,21x,2(i4,1x,i3,1x,i5,1x))', err=200) iprn,iyear0,idoy0,isod0,iyear1,idoy1,isod1
        call yeardoy2monthday(iyear0, idoy0, imon0, id0)
        jd0 = modified_julday(id0, imon0, iyear0)
        call yeardoy2monthday(iyear1, idoy1, imon1, id1)
        jd1 = modified_julday(id1, imon1, iyear1)
        if (jd1 + isod1/86400.d0 .le. fcbjd0 .or. jd0 + isod0/86400.d0 .ge. fcbjd1) cycle
        if (iprn .eq. '') cycle
        ii=pointer_string(MAXSAT,prn_mat,iprn)
        if (ii .eq. 0) cycle
        if (line(16:19) .ne. '') cycle
        read(line(26:28),*) obs_type
        if (iprn(1:1) .eq. 'G') then
          index_obs = index(obs_prio_G,obs_type(3:3))
          if(index_obs .eq. 0) cycle
          if(index(obs_type, 'L1') .ne. 0) then
            read (line, '(70x,f21.15)') bias(ii, index_obs)
            bias(ii, index_obs) = bias(ii, index_obs)*VLIGHT*1d-9
          elseif(index(obs_type, 'L2') .ne. 0) then
            read (line, '(70x,f21.15)') bias(ii, index_obs+11)
            bias(ii, index_obs+11) = bias(ii, index_obs+11)*VLIGHT*1d-9
          elseif(index(obs_type, 'C1') .ne. 0) then
            read (line, '(70x,f21.15)') bias(ii, index_obs+11*2)
            bias(ii, index_obs+11*2) = bias(ii, index_obs+11*2)*VLIGHT*1d-9
          elseif(index(obs_type, 'C2') .ne. 0) then
            read (line, '(70x,f21.15)') bias(ii, index_obs+11*3)
            bias(ii, index_obs+11*3) = bias(ii, index_obs+11*3)*VLIGHT*1d-9
          endif
        elseif (iprn(1:1) .eq. 'E') then
          index_obs = index(obs_prio_E,obs_type(3:3))
          if(index_obs .eq. 0) cycle
          if(index(obs_type, 'L1') .ne. 0) then
            read (line, '(70x,f21.15)') bias(ii, index_obs)
            bias(ii, index_obs) = bias(ii, index_obs)*VLIGHT*1d-9
          elseif(index(obs_type, 'L5') .ne. 0) then
            read (line, '(70x,f21.15)') bias(ii, index_obs+11)
            bias(ii, index_obs+11) = bias(ii, index_obs+11)*VLIGHT*1d-9
          elseif(index(obs_type, 'C1') .ne. 0) then
            read (line, '(70x,f21.15)') bias(ii, index_obs+11*2)
            bias(ii, index_obs+11*2) = bias(ii, index_obs+11*2)*VLIGHT*1d-9
          elseif(index(obs_type, 'C5') .ne. 0) then
            read (line, '(70x,f21.15)') bias(ii, index_obs+11*3)
            bias(ii, index_obs+11*3) = bias(ii, index_obs+11*3)*VLIGHT*1d-9
          endif
        elseif (iprn(1:1) .eq. 'C') then
          index_obs = index(obs_prio_C,obs_type(3:3))
          if(index_obs .eq. 0) cycle
          if(index(obs_type, 'L2') .ne. 0) then
            read (line, '(70x,f21.15)') bias(ii, index_obs)
            bias(ii, index_obs) = bias(ii, index_obs)*VLIGHT*1d-9
          elseif(index(obs_type, 'L6') .ne. 0) then
            read (line, '(70x,f21.15)') bias(ii, index_obs+11)
            bias(ii, index_obs+11) = bias(ii, index_obs+11)*VLIGHT*1d-9
          elseif(index(obs_type, 'C2') .ne. 0) then
            read (line, '(70x,f21.15)') bias(ii, index_obs+11*2)
            bias(ii, index_obs+11*2) = bias(ii, index_obs+11*2)*VLIGHT*1d-9
          elseif(index(obs_type, 'C6') .ne. 0) then
            read (line, '(70x,f21.15)') bias(ii, index_obs+11*3)
            bias(ii, index_obs+11*3) = bias(ii, index_obs+11*3)*VLIGHT*1d-9
          endif
          read (iprn(2:3), '(i2)') iprn_int
        elseif (iprn(1:1) .eq. 'J') then
          index_obs = index(obs_prio_J,obs_type(3:3))
          if(index_obs .eq. 0) cycle
          if(index(obs_type, 'L1') .ne. 0) then
            read (line, '(70x,f21.15)') bias(ii, index_obs)
            bias(ii, index_obs) = bias(ii, index_obs)*VLIGHT*1d-9
          elseif(index(obs_type, 'L2') .ne. 0) then
            read (line, '(70x,f21.15)') bias(ii, index_obs+11)
            bias(ii, index_obs+11) = bias(ii, index_obs+11)*VLIGHT*1d-9
          elseif(index(obs_type, 'C1') .ne. 0) then
            read (line, '(70x,f21.15)') bias(ii, index_obs+11*2)
            bias(ii, index_obs+11*2) = bias(ii, index_obs+11*2)*VLIGHT*1d-9
          elseif(index(obs_type, 'C2') .ne. 0) then
            read (line, '(70x,f21.15)') bias(ii, index_obs+11*3)
            bias(ii, index_obs+11*3) = bias(ii, index_obs+11*3)*VLIGHT*1d-9
          endif
        endif
      enddo
    endif
    read (lfn, '(a)', end=100) line
  enddo

100 if (index(line, '%=ENDBIA') .eq. 0) then 
    write (*, '(2a)') '***ERROR(read_fcb): end of file ', trim(flnfcb)
    call exit(1)
  else
    close (lfn)
    return
  endif
200 write (*, '(2a)') '***ERROR(read_fcb): read file ', trim(flnfcb)
  call exit(1)
end
