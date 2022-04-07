!
!! fixamb_rounding.f90
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
!! Contributor: Maorong Ge, Jianghui Geng, Songfeng Yang, Jihang Lin
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
  real*8 prob, alpha, val, sodb, sode, fwl, swl
  character*5 fixed
  type(ambssd) SDX, SD(MAXSD_SIT)

  nwnfx = 0
  nwlfx = 0
  ACF%ntot_G = 0
  ACF%ntot_E = 0
  ACF%ntot_C = 0
  ACF%ntot_3 = 0
  ACF%ntot_J = 0
  !! Statistics of each system fixing rate
  ACF%nwnfx_G = 0
  ACF%nwlfx_G = 0
  ACF%nwnfx_E = 0
  ACF%nwlfx_E = 0
  ACF%nwnfx_C = 0
  ACF%nwlfx_C = 0
  ACF%nwnfx_3 = 0
  ACF%nwlfx_3 = 0
  ACF%nwnfx_J = 0
  ACF%nwlfx_J = 0
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
        ACF%nwlfx_G = ACF%nwlfx_G + 1
      elseif(ASD(isd)%sys .eq. 'E')then
        ACF%nwlfx_E = ACF%nwlfx_E + 1
      elseif(ASD(isd)%sys .eq. 'C')then
        ACF%nwlfx_C = ACF%nwlfx_C + 1
      elseif(ASD(isd)%sys .eq. '3')then
        ACF%nwlfx_3 = ACF%nwlfx_3 + 1
      elseif(ASD(isd)%sys .eq. 'J')then
        ACF%nwlfx_J = ACF%nwlfx_J + 1
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
          ACF%nwnfx_G = ACF%nwnfx_G + 1
        elseif(ASD(isd)%sys .eq. 'E')then
          ACF%nwnfx_E = ACF%nwnfx_E + 1
        elseif(ASD(isd)%sys .eq. 'C')then
          ACF%nwnfx_C = ACF%nwnfx_C + 1
        elseif(ASD(isd)%sys .eq. '3')then
          ACF%nwnfx_3 = ACF%nwnfx_3 + 1
        elseif(ASD(isd)%sys .eq. 'J')then
          ACF%nwnfx_J = ACF%nwnfx_J + 1
        endif
      endif
!
!! output
      call timinc(ACF%jd0, ACF%sod0, (ASD(isd)%iepc(1) - 1)*ACF%dintv, jd, sodb)
      sodb = sodb + (jd-ACF%jd0)*86400.d0
      call timinc(ACF%jd0, ACF%sod0, (ASD(isd)%iepc(2) - 1)*ACF%dintv, jd, sode)
      sode = sode + (jd-ACF%jd0)*86400.d0
      write (*, '(i5,2(1x,a3),2(f14.3,f8.3),2f10.1,1x,a5)') isd, ACF%prn(isat), ACF%prn(jsat), &
        ASD(isd)%rwl, ASD(isd)%swl, ASD(isd)%rnl, ASD(isd)%snl, sodb, sode, fixed
    endif
  enddo
  !! Statistics of each system fixing rate
  write(*,'(a28,3a8,a7,2x,a7)') 'Integer rounding: ','#NL','#WL','#TOT','WL/TOT','NL/WL'
  ACF%ntot_G=AS%nsd_G
  if (ACF%ntot_G .ne. 0) then
    write (*, '(a,3i8,2(f6.1,a1,2x))') 'G_Wide/Narrow-lane FR(all): ', ACF%nwnfx_G, ACF%nwlfx_G, ACF%ntot_G, &
      ACF%nwlfx_G*1.d2/ACF%ntot_G, '%', ACF%nwnfx_G*1.d2/ACF%nwlfx_G, '%'
  endif
  ACF%ntot_E=AS%nsd_E
  if (ACF%ntot_E .ne. 0) then
    write (*, '(a,3i8,2(f6.1,a1,2x))') 'E_Wide/Narrow-lane FR(all): ', ACF%nwnfx_E, ACF%nwlfx_E, ACF%ntot_E, &
      ACF%nwlfx_E*1.d2/ACF%ntot_E, '%', ACF%nwnfx_E*1.d2/ACF%nwlfx_E, '%'
  endif
  ACF%ntot_C=AS%nsd_C
  if (ACF%ntot_C .ne. 0) then
    write (*, '(a,3i8,2(f6.1,a1,2x))') 'C2Wide/Narrow-lane FR(all): ', ACF%nwnfx_C, ACF%nwlfx_C, ACF%ntot_C, &
      ACF%nwlfx_C*1.d2/ACF%ntot_C, '%', ACF%nwnfx_C*1.d2/ACF%nwlfx_C, '%'
  endif
  ACF%ntot_3=AS%nsd_3
  if (ACF%ntot_3 .ne. 0) then
    write (*, '(a,3i8,2(f6.1,a1,2x))') 'C3Wide/Narrow-lane FR(all): ', ACF%nwnfx_3, ACF%nwlfx_3, ACF%ntot_3, &
      ACF%nwlfx_3*1.d2/ACF%ntot_3, '%', ACF%nwnfx_3*1.d2/ACF%nwlfx_3, '%'
  endif
  ACF%ntot_J=AS%nsd_J
  if (ACF%ntot_J .ne. 0) then
    write (*, '(a,3i8,2(f6.1,a1,2x))') 'J_Wide/Narrow-lane FR(all): ', ACF%nwnfx_J, ACF%nwlfx_J, ACF%ntot_J, &
      ACF%nwlfx_J*1.d2/ACF%ntot_J, '%', ACF%nwnfx_J*1.d2/ACF%nwlfx_J, '%'
  endif
  if (nwlfx .ne. 0) then
    write (*, '(a,3i8,2(f6.1,a1,2x))') 'mGWide/Narrow-lane FR(all): ', nwnfx, nwlfx, AS%nsd, &
        nwlfx*1.d2/AS%nsd, '%', nwnfx*1.d2/nwlfx, '%'
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
  if (ACF%lsearch) then
    ACF%ntot_G = 0
    ACF%ntot_E = 0
    ACF%ntot_C = 0
    ACF%ntot_3 = 0
    ACF%ntot_J = 0
    ACF%nwlfx_G = 0
    ACF%nwlfx_E = 0
    ACF%nwlfx_C = 0
    ACF%nwlfx_3 = 0
    ACF%nwlfx_J = 0
  endif
  ndef = 0
  do k = 1, 3
    do i = 1, AS%nsd
      if (k .eq. 1 .and. ASD(i)%id .ne. 0) cycle
      if (k .eq. 2 .and. ASD(i)%id .ne. 1) cycle
      if (k .eq. 3 .and. ASD(i)%id .ne. 2) cycle
      call check_amb_depend(AS%now, ndef, ASD(i)%ipt, ldep)
      if (ldep) cycle
      if (ACF%lsearch) then
        if(ASD(i)%sys.eq.'G') then
          ACF%ntot_G = ACF%ntot_G + 1
          if(k.le.2) ACF%nwlfx_G = ACF%nwlfx_G + 1
        else if(ASD(i)%sys.eq.'E') then
          ACF%ntot_E = ACF%ntot_E + 1
          if(k.le.2) ACF%nwlfx_E = ACF%nwlfx_E + 1
        else if(ASD(i)%sys.eq.'C') then
          ACF%ntot_C = ACF%ntot_C + 1
          if(k.le.2) ACF%nwlfx_C = ACF%nwlfx_C + 1
        else if(ASD(i)%sys.eq.'3') then
          ACF%ntot_3 = ACF%ntot_3 + 1
          if(k.le.2) ACF%nwlfx_3 = ACF%nwlfx_3 + 1
        else if(ASD(i)%sys.eq.'J') then
          ACF%ntot_J = ACF%ntot_J + 1
          if(k.le.2) ACF%nwlfx_J = ACF%nwlfx_J + 1
        endif
      endif
      if (k .le. 2) SD(ndef) = ASD(i)
    enddo
    if (k .eq. 1) nwnfx = ndef
    if (k .eq. 2) nwlfx = ndef
    if (k .eq. 3) ntot = ndef
  enddo
  call check_amb_depend(0, 0, (/1, 1/), ldep)
  write (*, '(a,3i8,2(f6.1,a1,2x))') 'mGWide/Narrow-lane FR(ind): ', nwnfx, nwlfx, ntot, &
      nwlfx*1.d2/ntot, '%', nwnfx*1.d2/nwlfx, '%'
  ACF%ntot_ind=ntot
  ACF%nwl_ind =nwlfx
  ACF%nwn_ind =nwnfx
!! only save widelane fixed ones
  AS%nsd = nwlfx
  do i = 1, AS%nsd
    ASD(i) = SD(i)
  enddo
  return
end
