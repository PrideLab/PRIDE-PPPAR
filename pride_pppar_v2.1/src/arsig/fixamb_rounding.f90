!
!! fixamb_rounding.f90
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
!! purpose  : fix ambiguity at one station
!! parameter:
!!    input : ACF  -- fractional cycle biases
!!    output: AS   -- ambiguity station struct
!!            ASD  -- single-difference struct
!
subroutine fixamb_rounding(ACF, AS, ASD)
  implicit none
  include '../header/const.h'
  include 'ambssd.h'
  include 'arscfg.h'
  include 'ambsit.h'

  type(arscfg) ACF
  type(ambsit) AS
  type(ambssd) ASD(1:*)
!
!! local
  logical*1 ldep
  integer*4 i, j, k, l, jd, ik, il, isat, jsat, isd, nwnfx, nwlfx, ntot, ndef
  integer*4 nwnfx_G, nwlfx_G, nwnfx_R, nwlfx_R, nwnfx_E, nwlfx_E, nwnfx_C, nwlfx_C, nwnfx_3, nwlfx_3, nwnfx_J, nwlfx_J
  real*8 prob, alpha, val, sodb, sode, fwl, swl
  character*5 fixed
  type(ambssd) SDX, SD(MAXSD_SIT)

  nwnfx = 0
  nwlfx = 0
  !! Statistics of each system fixing rate
  nwnfx_G = 0
  nwlfx_G = 0
  nwnfx_R = 0
  nwlfx_R = 0
  nwnfx_E = 0
  nwlfx_E = 0
  nwnfx_C = 0
  nwlfx_C = 0
  nwnfx_3 = 0
  nwlfx_3 = 0
  nwnfx_J = 0
  nwlfx_J = 0
  do isd = 1, AS%nsd
    fixed = ' '
    ASD(isd)%dec = 0.d0
!
!! pointer to satellites
    isat = AS%isat(ASD(isd)%ipt(1))
    jsat = AS%isat(ASD(isd)%ipt(2))
!
!! try to fix
    call bdeci(ASD(isd)%rwl, ASD(isd)%swl, 1, ACF%wl_maxdev, ACF%wl_maxsig, prob, alpha)
    if (prob .gt. 0.d0 .and. alpha .gt. ACF%wl_alpha) then
      ASD(isd)%id = 1               ! wide-lane fixed
      fixed = 'WL'
      nwlfx = nwlfx + 1
      !! Statistics of each system fixing rate
      if(ASD(isd)%sys .eq. 'G')then
        nwlfx_G = nwlfx_G + 1
      elseif(ASD(isd)%sys .eq. 'E')then
        nwlfx_E = nwlfx_E + 1
      elseif(ASD(isd)%sys .eq. 'C')then
        nwlfx_C = nwlfx_C + 1
      elseif(ASD(isd)%sys .eq. '3')then
        nwlfx_3 = nwlfx_3 + 1
      elseif(ASD(isd)%sys .eq. 'J')then
        nwlfx_J = nwlfx_J + 1
      endif
!
!! compute narrow-lane ambiguities
      if(ASD(isd)%sys .eq. 'G')then
        ASD(isd)%rnl = (FREQ1_G+FREQ2_G)/FREQ1_G*(AS%xamb(ASD(isd)%ipt(1)) - AS%xamb(ASD(isd)%ipt(2))) &
                       - FREQ2_G/(FREQ1_G-FREQ2_G)*nint(ASD(isd)%rwl)
      elseif(ASD(isd)%sys .eq. 'E')then
        ASD(isd)%rnl = (FREQ1_E+FREQ2_E)/FREQ1_E*(AS%xamb(ASD(isd)%ipt(1)) - AS%xamb(ASD(isd)%ipt(2))) &
                       - FREQ2_E/(FREQ1_E-FREQ2_E)*nint(ASD(isd)%rwl)
      elseif(ASD(isd)%sys .eq. 'C' .or. ASD(isd)%sys .eq. '3')then
        ASD(isd)%rnl = (FREQ1_C+FREQ2_C)/FREQ1_C*(AS%xamb(ASD(isd)%ipt(1)) - AS%xamb(ASD(isd)%ipt(2))) &
                       - FREQ2_C/(FREQ1_C-FREQ2_C)*nint(ASD(isd)%rwl)
      elseif(ASD(isd)%sys .eq. 'J')then
        ASD(isd)%rnl = (FREQ1_J+FREQ2_J)/FREQ1_J*(AS%xamb(ASD(isd)%ipt(1)) - AS%xamb(ASD(isd)%ipt(2))) &
                       - FREQ2_J/(FREQ1_J-FREQ2_J)*nint(ASD(isd)%rwl)
      endif
      ASD(isd)%snl = 0.05d0 ! presumed statistics, see Ge et al. (2005)
!
!! try to fix
      call bdeci(ASD(isd)%rnl, ASD(isd)%snl, 1, ACF%nl_maxdev, ACF%nl_maxsig, prob, alpha)
      ASD(isd)%dec = alpha
      if (prob .gt. 0.d0 .and. alpha .gt. ACF%nl_alpha) then
        ASD(isd)%id = 0                     ! narrow-lane fixed already
        fixed = 'WL_NL'
        nwnfx = nwnfx + 1
        !! Statistics of each system fixing rate
        if(ASD(isd)%sys .eq. 'G')then
          nwnfx_G = nwnfx_G + 1
        elseif(ASD(isd)%sys .eq. 'E')then
          nwnfx_E = nwnfx_E + 1
        elseif(ASD(isd)%sys .eq. 'C')then
          nwnfx_C = nwnfx_C + 1
        elseif(ASD(isd)%sys .eq. '3')then
          nwnfx_3 = nwnfx_3 + 1
        elseif(ASD(isd)%sys .eq. 'J')then
          nwnfx_J = nwnfx_J + 1
        endif
      endif
!
!! output
      call timinc(ACF%jd0, ACF%sod0, (ASD(isd)%iepc(1) - 1)*ACF%dintv, jd, sodb)
      call timinc(ACF%jd0, ACF%sod0, (ASD(isd)%iepc(2) - 1)*ACF%dintv, jd, sode)
      write (*, '(i5,1x,a4,2a3,2(f14.3,f8.3),2f10.1,1x,a5)') isd, AS%name, ACF%prn(isat), ACF%prn(jsat), &
        ASD(isd)%rwl, ASD(isd)%swl, ASD(isd)%rnl, ASD(isd)%snl, sodb, sode, fixed
    endif
  enddo
  !! Statistics of each system fixing rate
  write(*,'(a28,4x,3a8,a7,2x,a7)') 'Integer rounding: ','#NL','#WL','#TOT','WL/TOT','NL/TOT'
  if (nwlfx_G .ne. 0) then
    write (*, '(a,a4,3i8,2(f6.1,a1,2x))') 'G_Wide/Narrow-lane FR(all): ', AS%name, nwnfx_G, nwlfx_G, AS%nsd_G, &
      nwlfx_G*1.d2/AS%nsd_G, '%', nwnfx_G*1.d2/AS%nsd_G, '%'
  endif
  if (nwlfx_R .ne. 0) then
    write (*, '(a,a4,3i8,2(f6.1,a1,2x))') 'R_Wide/Narrow-lane FR(all): ', AS%name, nwnfx_R, nwlfx_R, AS%nsd_R, &
      nwlfx_R*1.d2/AS%nsd_R, '%', nwnfx_R*1.d2/AS%nsd_R, '%'
  endif
  if (nwlfx_E .ne. 0) then
    write (*, '(a,a4,3i8,2(f6.1,a1,2x))') 'E_Wide/Narrow-lane FR(all): ', AS%name, nwnfx_E, nwlfx_E, AS%nsd_E, &
      nwlfx_E*1.d2/AS%nsd_E, '%', nwnfx_E*1.d2/AS%nsd_E, '%'
  endif
  if (nwlfx_C .ne. 0) then
    write (*, '(a,a4,3i8,2(f6.1,a1,2x))') 'C2Wide/Narrow-lane FR(all): ', AS%name, nwnfx_C, nwlfx_C, AS%nsd_C, &
      nwlfx_C*1.d2/AS%nsd_C, '%', nwnfx_C*1.d2/AS%nsd_C, '%'
  endif
  if (nwlfx_3 .ne. 0) then
    write (*, '(a,a4,3i8,2(f6.1,a1,2x))') 'C3Wide/Narrow-lane FR(all): ', AS%name, nwnfx_3, nwlfx_3, AS%nsd_3, &
      nwlfx_3*1.d2/AS%nsd_3, '%', nwnfx_3*1.d2/AS%nsd_3, '%'
  endif
  if (nwlfx_J .ne. 0) then
    write (*, '(a,a4,3i8,2(f6.1,a1,2x))') 'J_Wide/Narrow-lane FR(all): ', AS%name, nwnfx_J, nwlfx_J, AS%nsd_J, &
      nwlfx_J*1.d2/AS%nsd_J, '%', nwnfx_J*1.d2/AS%nsd_J, '%'
  endif
  if (nwlfx .ne. 0) then
    write (*, '(a,a4,3i8,2(f6.1,a1,2x))') 'mGWide/Narrow-lane FR(all): ', AS%name, nwnfx, nwlfx, AS%nsd, &
        nwlfx*1.d2/AS%nsd, '%', nwnfx*1.d2/AS%nsd, '%'
  endif
!
!! sort ASD
  do k = 1, 2
    do i = 1, AS%nsd - 1
      if (k .eq. 1 .and. ASD(i)%id .ne. 0) cycle
      if (k .eq. 2 .and. ASD(i)%id .ne. 1) cycle
      do j = i + 1, AS%nsd
        if (k .eq. 1 .and. ASD(j)%id .ne. 0) cycle
        if (k .eq. 2 .and. ASD(j)%id .ne. 1) cycle
        if (ASD(i)%dec .lt. ASD(j)%dec) then
          SDX = ASD(i)
          ASD(i) = ASD(j)
          ASD(j) = SDX
        endif
      enddo
    enddo
  enddo
!
!! select independent SD for this station
  ndef = 0
  do k = 1, 3
    do i = 1, AS%nsd
      if (k .eq. 1 .and. ASD(i)%id .ne. 0) cycle
      if (k .eq. 2 .and. ASD(i)%id .ne. 1) cycle
      if (k .eq. 3 .and. ASD(i)%id .ne. 2) cycle
      call check_amb_depend(AS%now, ndef, ASD(i)%ipt, ldep)
      if (ldep) cycle
      if (k .le. 2) SD(ndef) = ASD(i)
    enddo
    if (k .eq. 1) nwnfx = ndef
    if (k .eq. 2) nwlfx = ndef
    if (k .eq. 3) ntot = ndef
  enddo
  call check_amb_depend(0, 0, (/1, 1/), ldep)
  write (*, '(a,a4,3i8,2(f6.1,a1,2x))') 'mGWide/Narrow-lane FR(ind): ', AS%name, nwnfx, nwlfx, ntot, &
      nwlfx*1.d2/ntot, '%', nwnfx*1.d2/ntot, '%'
!! only save widelane fixed ones
  AS%nsd = nwlfx
  do i = 1, AS%nsd
    ASD(i) = SD(i)
  enddo
  return
end
