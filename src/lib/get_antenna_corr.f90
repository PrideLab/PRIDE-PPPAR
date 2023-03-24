!
!! get_ant_corr.f90
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
!! Contributor: Maorong Ge, Jianghui Geng, Songfeng Yang, Jihang Lin, Jing Zeng
!! 
!!
!!
!! purpose   : get antenna pointer in the whole table for receivers and satellites
!! parameter :
!!    input  : fjd_beg,fjd_end -- time span
!!             antnam,antnum   -- antenna name and serial number
!!    output : iptatx  -- pointer to global table
!!             enu     -- antenna offset
!
subroutine get_ant_ipt(fjd_beg, fjd_end, antnam, antnum, iptatx, enu, system, atx_snxcode)
  implicit none
  include '../header/const.h'
  include '../header/antatx.h'

  integer*4 iptatx, ript, sipt
  real*8 fjd_beg, fjd_end, zeni, azim, nadir, enu(1:*), var(2)
  character*20 antnam, antnum
  integer*4 sys1,sys2,frequency1,frequency2
  character*1 system
!
!! local
  integer*4 i, j, natx, izen, iazi
  real*8 rad2deg, alpha, zen, azi, nad, x1, x2
  type(antatx) AX(1 + MAXSAT), ATX
  character*60 atx_snxcode

  data natx/0/
  save natx, AX
  
  do i = 1, natx
    if (antnam .ne. AX(i)%antnam) cycle
    j = len_trim(antnum)
    if (j .ne. 0 .and. antnum(1:j) .ne. AX(i)%antnum(1:j)) cycle
    iptatx = i
    exit
  enddo
  if (iptatx .ne. 0) goto 5
!
!! if not in memory, read from file
  ATX%antnam = antnam
  ATX%antnum = antnum
  call rdatx(fjd_beg, fjd_end, ATX, atx_snxcode)
  natx = natx + 1
  AX(natx) = ATX
  antnam = ATX%antnam
  iptatx = natx
!
!! get antenna phase offset
  5 continue
  sys1=0
  sys2=0
  frequency1=0
  frequency2=0
  if(system .eq. 'G')then
    sys1=1
    sys2=1
    frequency1=1
    frequency2=2
  elseif(system .eq. 'R')then
    sys1=2
    sys2=2
    frequency1=1
    frequency2=2
    if(index(AX(iptatx)%sys_multi,'R') .eq. 0)then
      sys1=1
      frequency1=1
    endif
    if(index(AX(iptatx)%sys_multi2,'R') .eq. 0)then
      sys2=1
      frequency2=2
    endif
  elseif(system .eq. 'E')then
    sys1=3
    sys2=3
    frequency1=1
    frequency2=5
    if(index(AX(iptatx)%sys_multi,'E') .eq. 0)then
      sys1=1
      frequency1=1
    endif
    if(index(AX(iptatx)%sys_multi2,'E') .eq. 0)then
      sys2=1
      frequency2=2
    endif
  elseif(system .eq. 'C')then
    sys1=4
    sys2=4
    frequency1=2
    frequency2=6
    if(index(AX(iptatx)%sys_multi,'C') .eq. 0)then
      sys1=1
      frequency1=1
    endif
    if(index(AX(iptatx)%sys_multi2,'C') .eq. 0)then
      sys2=1
      frequency2=2
    endif
  elseif(system .eq. 'J')then
    sys1=5
    sys2=5
    frequency1=1
    frequency2=2
    if(index(AX(iptatx)%sys_multi,'J') .eq. 0)then
      sys1=1
      frequency1=1
    endif
    if(index(AX(iptatx)%sys_multi2,'J') .eq. 0)then
      sys2=1
      frequency2=2
    endif

  endif
  if(antnam(1:5).ne.'BLOCK' .and. antnam(1:7).ne.'GLONASS' .and. antnam(1:7).ne.'GALILEO' &
     .and. antnam(1:6).ne.'BEIDOU' .and. antnam(1:4).ne.'QZSS') then
    enu(1) = AX(iptatx)%neu(2, frequency1, sys1)
    enu(2) = AX(iptatx)%neu(1, frequency1, sys1)
    enu(3) = AX(iptatx)%neu(3, frequency1, sys1)
    enu(4) = AX(iptatx)%neu(2, frequency2, sys2)
    enu(5) = AX(iptatx)%neu(1, frequency2, sys2)
    enu(6) = AX(iptatx)%neu(3, frequency2, sys2)
  else
    do i = 1, 3
      enu(i) = AX(iptatx)%neu(i, frequency1, sys1)
      enu(i + 3) = AX(iptatx)%neu(i, frequency2, sys2)
    enddo
  endif

  return
!
!! purpose   : get antenna pcv for receivers and satellites
!! parameter :
!!    input  : ript,sipt -- pointer of receiver and satellite
!!             zeni  -- zenith angle in radian
!!             azim  -- azimuth angle in radian
!!             nadir -- nadir angle in radian
!!    output : var   -- phase center variation
!
  Entry get_ant_pcv(ript, sipt, zeni, azim, nadir, var, system) ! add multisystem-GREC
!
!! constant
  rad2deg = 180.d0/PI
  var = 0.d0
!
!! receiver pcv index
  if (ript .le. 0) goto 10
  sys1=0
  sys2=0
  frequency1=0
  frequency2=0
  if(system .eq. 'G')then
    sys1=1
    sys2=1
    frequency1=1
    frequency2=2
  elseif(system .eq. 'R')then
    sys1=2
    sys2=2
    frequency1=1
    frequency2=2
    if(index(AX(ript)%sys_multi,'R') .eq. 0)then
      sys1=1
      frequency1=1
    endif
    if(index(AX(ript)%sys_multi2,'R') .eq. 0)then
      sys2=1
      frequency2=2
    endif
  elseif(system .eq. 'E')then
    sys1=3
    sys2=3
    frequency1=1
    frequency2=5
    if(index(AX(ript)%sys_multi,'E') .eq. 0)then
      sys1=1
      frequency1=1
    endif
    if(index(AX(ript)%sys_multi2,'E') .eq. 0)then
      sys2=1
      frequency2=2
    endif
  elseif(system .eq. 'C')then
    sys1=4
    sys2=4
    frequency1=2
    frequency2=6
    if(index(AX(ript)%sys_multi,'C') .eq. 0)then
      sys1=1
      frequency1=1
    endif
    if(index(AX(ript)%sys_multi2,'C') .eq. 0)then
      sys2=1
      frequency2=2
    endif
  elseif(system .eq. 'J')then
    sys1=5
    sys2=5
    frequency1=1
    frequency2=2
    if(index(AX(ript)%sys_multi,'J') .eq. 0)then
      sys1=1
      frequency1=1
    endif
    if(index(AX(ript)%sys_multi2,'J') .eq. 0)then
      sys2=1
      frequency2=2
    endif
  endif
!
!! receiver pcv
  zen = zeni*rad2deg
  azi = azim*rad2deg
  if (zen .gt. AX(ript)%zen2) zen = AX(ript)%zen2
  if (zen .lt. AX(ript)%zen1) zen = AX(ript)%zen1
  if (azi .lt. 0.d0) azi = azi + 360.d0
!
!! azimuth dependent
  iazi = 0
  if (AX(ript)%dazi .ne. 0.d0) iazi = int(azi/AX(ript)%dazi) + 1
!
!! zenith dependent
  izen = int((zen - AX(ript)%zen1)/AX(ript)%dzen) + 1

  x1 = AX(ript)%pcv(izen, iazi, frequency1, sys1)
  x2 = AX(ript)%pcv(izen + 1, iazi, frequency1, sys1)
  if (iazi .ne. 0) then
    alpha = azi/AX(ript)%dazi - iazi + 1
    x1 = x1 + (AX(ript)%pcv(izen, iazi + 1, frequency1, sys1) - AX(ript)%pcv(izen, iazi, frequency1, sys1))*alpha
    x2 = x2 + (AX(ript)%pcv(izen + 1, iazi + 1, frequency1, sys1) - AX(ript)%pcv(izen + 1, iazi, frequency1, sys1))*alpha
  endif
  alpha = (zen - AX(ript)%zen1)/AX(ript)%dzen - izen + 1
  var(1) = var(1) + x1 + (x2 - x1)*alpha
  
  x1 = AX(ript)%pcv(izen, iazi, frequency2, sys2)
  x2 = AX(ript)%pcv(izen + 1, iazi, frequency2, sys2)
  if (iazi .ne. 0) then
    alpha = azi/AX(ript)%dazi - iazi + 1
    x1 = x1 + (AX(ript)%pcv(izen, iazi + 1, frequency2, sys2) - AX(ript)%pcv(izen, iazi, frequency2, sys2))*alpha
    x2 = x2 + (AX(ript)%pcv(izen + 1, iazi + 1, frequency2, sys2) - AX(ript)%pcv(izen + 1, iazi, frequency2, sys2))*alpha
  endif
  alpha = (zen - AX(ript)%zen1)/AX(ript)%dzen - izen + 1
  var(2) = var(2) + x1 + (x2 - x1)*alpha

!
!! satellite pcv index
 10 continue
  if (sipt .le. 0) goto 20
  sys1=0
  sys2=0
  frequency1=0
  frequency2=0
  if(system .eq. 'G')then
    sys1=1
    sys2=1
    frequency1=1
    frequency2=2
  elseif(system .eq. 'R')then
    sys1=2
    sys2=2
    frequency1=1
    frequency2=2
    if(index(AX(sipt)%sys_multi,'R') .eq. 0)then
      sys1=1
      frequency1=1
    endif
    if(index(AX(sipt)%sys_multi2,'R') .eq. 0)then
      sys2=1
      frequency2=2
    endif
  elseif(system .eq. 'E')then
    sys1=3
    sys2=3
    frequency1=1
    frequency2=5
    if(index(AX(sipt)%sys_multi,'E') .eq. 0)then
      sys1=1
      frequency1=1
    endif
    if(index(AX(sipt)%sys_multi2,'E') .eq. 0)then
      sys2=1
      frequency2=2
    endif
  elseif(system .eq. 'C')then
    sys1=4
    sys2=4
    frequency1=2
    frequency2=6
    if(index(AX(sipt)%sys_multi,'C') .eq. 0)then
      sys1=1
      frequency1=1
    endif
    if(index(AX(sipt)%sys_multi2,'C') .eq. 0)then
      sys2=1
      frequency2=2
    endif
  elseif(system .eq. 'J')then
    sys1=5
    sys2=5
    frequency1=1
    frequency2=2
    if(index(AX(sipt)%sys_multi,'J') .eq. 0)then
      sys1=1
      frequency1=1
    endif
    if(index(AX(sipt)%sys_multi2,'J') .eq. 0)then
      sys2=1
      frequency2=2
    endif
  endif
!
!! satellite pcv
  if (AX(sipt)%dzen .eq. 0.d0) return
  nad = nadir*rad2deg
  if (nad .gt. AX(sipt)%zen2) nad = AX(sipt)%zen2
  if (nad .lt. AX(sipt)%zen1) nad = AX(sipt)%zen1
!
!! only zenith dependent
  izen = int((nad - AX(sipt)%zen1)/AX(sipt)%dzen) + 1

  x1 = AX(sipt)%pcv(izen, 0, frequency1, sys1)
  x2 = AX(sipt)%pcv(izen + 1, 0, frequency1, sys1)
  alpha = (nad - AX(sipt)%zen1)/AX(sipt)%dzen - izen + 1
  var(1) = var(1) + x1 + (x2 - x1)*alpha
  
  x1 = AX(sipt)%pcv(izen, 0, frequency2, sys2)
  x2 = AX(sipt)%pcv(izen + 1, 0, frequency2, sys2)
  alpha = (nad - AX(sipt)%zen1)/AX(sipt)%dzen - izen + 1
  var(2) = var(2) + x1 + (x2 - x1)*alpha

 20 continue
  return
end subroutine
