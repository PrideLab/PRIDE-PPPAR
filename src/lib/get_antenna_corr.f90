!
!! get_ant_corr.f90
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
subroutine get_ant_ipt(fjd_beg, fjd_end, antnam, antnum, iptatx, enu, sys, atx_snxcode)
  implicit none
  include '../header/const.h'
  include '../header/antatx.h'

! common
  integer*4     idxfrq(MAXSYS, 2)
  common        idxfrq
! parameter
  real*8        fjd_beg, fjd_end
  character*20  antnam, antnum
  integer*4     iptatx
  real*8        enu(3, 2)
  character*1   sys
! local
  integer*4     i0, i, j, natx, izen, iazi
  integer*4     ifq1, ifq2
  integer*4     ript, sipt
  real*8        zeni, azim,azim_sat, nadir, pcv(2)
  real*8        rad2deg, alpha, zen, azi, nad, x1, x2
  type(antatx)  AX(1 + MAXSAT), ATX
  character*60  atx_snxcode
  character*3   cfq1, cfq2
!
!! function called
  integer*4     pointer_string

  data natx/0/
  save natx, AX

  iptatx = 0
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

  i0 = index(GNSS_PRIO, sys)
  if (pointer_string(MAXFRQ, AX(iptatx)%frq(:, i0), sys) .eq. 0) then
    i0 = index(GNSS_PRIO, 'G')
    ifq1 = 1
    ifq2 = 2
  end if
  write (cfq1, '(a1,i0.2)') GNSS_PRIO(i0:i0), idxfrq(i0, 1)
  write (cfq2, '(a1,i0.2)') GNSS_PRIO(i0:i0), idxfrq(i0, 2)
  ifq1 = pointer_string(MAXFRQ, AX(iptatx)%frq(:, i0), cfq1)
  ifq2 = pointer_string(MAXFRQ, AX(iptatx)%frq(:, i0), cfq2)

  if (antnam(1:5) .ne. 'BLOCK' .and. antnam(1:7) .ne. 'GLONASS' .and. antnam(1:7) .ne. 'GALILEO' &
     .and. antnam(1:6) .ne. 'BEIDOU' .and. antnam(1:4) .ne. 'QZSS') then
    !
    !! reciver antenna PCO
    if (ifq1 .eq. 0 .or. ifq2 .eq. 0) then
      i0 = index(GNSS_PRIO, 'G')
      ifq1 = 1
      ifq2 = 2
    end if
    enu(1, 1) = AX(iptatx)%neu(2, ifq1, i0)
    enu(2, 1) = AX(iptatx)%neu(1, ifq1, i0)
    enu(3, 1) = AX(iptatx)%neu(3, ifq1, i0)
    enu(1, 2) = AX(iptatx)%neu(2, ifq2, i0)
    enu(2, 2) = AX(iptatx)%neu(1, ifq2, i0)
    enu(3, 2) = AX(iptatx)%neu(3, ifq2, i0)
  else
    if (ifq2 .eq. 0) then
      ifq2 = 2
      if (cfq2 .eq. 'C05' .or. cfq2 .eq. 'C08') ifq2 = 7
    end if
    enu(1:3, 1) = AX(iptatx)%neu(1:3, ifq1, i0)
    enu(1:3, 2) = AX(iptatx)%neu(1:3, ifq2, i0)
  endif

  return
!
!! purpose   : get antenna pcv for receivers and satellites
!! parameter :
!!    input  : ript,sipt -- pointer of receiver and satellite
!!             zeni  -- zenith angle in radian
!!             azim  -- azimuth angle in radian
!!             nadir -- nadir angle in radian
!!    output : pcv   -- phase center variation
!
  Entry get_ant_pcv(ript, sipt, zeni, azim, azim_sat,  nadir, pcv, sys) ! add multisys-GREC

  rad2deg = 180.d0/PI
  pcv = 0.d0

  zen = zeni*rad2deg
  azi = azim*rad2deg
!
!! receiver frequency indecies
  i0 = index(GNSS_PRIO, sys)
  if (pointer_string(MAXFRQ, AX(ript)%frq(:, i0), sys) .eq. 0) then
    i0 = index(GNSS_PRIO, 'G')
    ifq1 = 1
    ifq2 = 2
  end if
  write (cfq1, '(a1,i0.2)') GNSS_PRIO(i0:i0), idxfrq(i0, 1)
  write (cfq2, '(a1,i0.2)') GNSS_PRIO(i0:i0), idxfrq(i0, 2)
  ifq1 = pointer_string(MAXFRQ, AX(ript)%frq(:, i0), cfq1)
  ifq2 = pointer_string(MAXFRQ, AX(ript)%frq(:, i0), cfq2)
  if (ifq1 .eq. 0 .or. ifq2 .eq. 0) then
    i0 = index(GNSS_PRIO, 'G')
    ifq1 = 1
    ifq2 = 2
  end if

!
!! receiver antenna PCVs
  if (zen .gt. AX(ript)%zen2) zen = AX(ript)%zen2
  if (zen .lt. AX(ript)%zen1) zen = AX(ript)%zen1
  if (azi .lt. 0.d0) azi = azi + 360.d0

  iazi = 0
  if (AX(ript)%dazi .ne. 0.d0) iazi = int(azi/AX(ript)%dazi) + 1

  izen = int((zen - AX(ript)%zen1)/AX(ript)%dzen) + 1

  x1 = AX(ript)%pcv(izen,   iazi, ifq1, i0)
  x2 = AX(ript)%pcv(izen+1, iazi, ifq1, i0)
  if (iazi .ne. 0) then
    alpha = azi/AX(ript)%dazi - iazi + 1
    x1 = x1 + (AX(ript)%pcv(izen, iazi + 1, ifq1, i0) - AX(ript)%pcv(izen, iazi, ifq1, i0))*alpha
    x2 = x2 + (AX(ript)%pcv(izen + 1, iazi + 1, ifq1, i0) - AX(ript)%pcv(izen + 1, iazi, ifq1, i0))*alpha
  endif
  alpha = (zen - AX(ript)%zen1)/AX(ript)%dzen - izen + 1
  pcv(1) = pcv(1) + x1 + (x2 - x1)*alpha

  x1 = AX(ript)%pcv(izen,   iazi, ifq2, i0)
  x2 = AX(ript)%pcv(izen+1, iazi, ifq2, i0)
  if (iazi .ne. 0) then
    alpha = azi/AX(ript)%dazi - iazi + 1
    x1 = x1 + (AX(ript)%pcv(izen, iazi + 1, ifq2, i0) - AX(ript)%pcv(izen, iazi, ifq2, i0))*alpha
    x2 = x2 + (AX(ript)%pcv(izen + 1, iazi + 1, ifq2, i0) - AX(ript)%pcv(izen + 1, iazi, ifq2, i0))*alpha
  endif
  alpha = (zen - AX(ript)%zen1)/AX(ript)%dzen - izen + 1
  pcv(2) = pcv(2) + x1 + (x2 - x1)*alpha

!
!! satellite frequency indecies
  i0 = index(GNSS_PRIO, sys)
  if (pointer_string(MAXFRQ, AX(sipt)%frq(:, i0), sys) .eq. 0) then
    i0 = index(GNSS_PRIO, 'G')
    ifq1 = 1
    ifq2 = 2
  end if
  write (cfq1, '(a1,i0.2)') GNSS_PRIO(i0:i0), idxfrq(i0, 1)
  write (cfq2, '(a1,i0.2)') GNSS_PRIO(i0:i0), idxfrq(i0, 2)
  ifq1 = pointer_string(MAXFRQ, AX(sipt)%frq(:, i0), cfq1)
  ifq2 = pointer_string(MAXFRQ, AX(sipt)%frq(:, i0), cfq2)
  if (ifq1 .eq. 0 .or. ifq2 .eq. 0) then
    i0 = index(GNSS_PRIO, 'G')
    ifq1 = 1
    ifq2 = 2
  end if

!
!! satellite antenna PCVs
  if (AX(sipt)%dzen .eq. 0.d0) return
  nad = nadir*rad2deg
  azi = azim_sat*rad2deg
  if (nad .gt. AX(sipt)%zen2) nad = AX(sipt)%zen2
  if (nad .lt. AX(sipt)%zen1) nad = AX(sipt)%zen1

  iazi=0
  if (AX(sipt)%dazi .ne. 0.d0) iazi = int(azi/AX(sipt)%dazi) + 1  

  izen = int((nad - AX(sipt)%zen1)/AX(sipt)%dzen) + 1

  x1 = AX(sipt)%pcv(izen,   iazi, ifq1, i0)
  x2 = AX(sipt)%pcv(izen+1, iazi, ifq1, i0)
  if (iazi.ne.0) then
    alpha = azi/AX(sipt)%dazi - iazi +1
    x1 = x1 + (AX(sipt)%pcv(izen, iazi + 1, ifq1, i0) - AX(sipt)%pcv(izen, iazi, ifq1, i0))*alpha
    x2 = x2 + (AX(sipt)%pcv(izen + 1, iazi + 1, ifq1, i0) - AX(sipt)%pcv(izen + 1, iazi, ifq1, i0))*alpha
  endif
  alpha = (nad - AX(sipt)%zen1)/AX(sipt)%dzen - izen + 1
  pcv(1) = pcv(1) + x1 + (x2 - x1)*alpha
  
  x1 = AX(sipt)%pcv(izen,   iazi, ifq2, i0)
  x2 = AX(sipt)%pcv(izen+1, iazi, ifq2, i0)
  if (iazi .ne.0 ) then
    alpha = azi/AX(sipt)%dazi - iazi + 1
    x1 = x1 + (AX(sipt)%pcv(izen, iazi + 1, ifq2, i0) - AX(sipt)%pcv(izen, iazi, ifq2, i0))*alpha
    x2 = x2 + (AX(sipt)%pcv(izen + 1, iazi + 1, ifq2, i0) - AX(sipt)%pcv(izen + 1, iazi, ifq2, i0))*alpha
  endif
  alpha = (nad - AX(sipt)%zen1)/AX(sipt)%dzen - izen + 1
  pcv(2) = pcv(2) + x1 + (x2 - x1)*alpha

  return
end subroutine
