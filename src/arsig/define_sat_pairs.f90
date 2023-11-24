!
!! define_sat_pairs.f90
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
!! Contributor: Maorong Ge, Jianghui Geng
!! 
!!
!!
!! purpose  : define satellite pairs per station
!! parameter:
!!    input : ACF -- fractional cycle biases
!!            AS  -- station ambiguity
!!    output: ASD -- single-difference ambiguities
!
subroutine define_sat_pairs(ACF, AS, ASD)
  implicit none
  include '../header/const.h'
  include 'ambssd.h'
  include 'arscfg.h'
  include 'ambsit.h'

  type(arscfg) ACF
  type(ambsit) AS
  type(ambssd) ASD(MAXSD_SIT)
!
!! local
  integer*4 j, k, isat, jsat, im, in, is
  integer*4 iprn_int,jprn_int
!
!! function used
  integer*4 pointer_int, pointer_string

  AS%nsd   = 0
  AS%nsd_G = 0
  AS%nsd_R = 0
  AS%nsd_E = 0
  AS%nsd_C = 0
  AS%nsd_3 = 0
  AS%nsd_J = 0

  do im=1,AS%now-1
    if(pointer_string(ACF%fcbnprn, ACF%fcbprn, ACF%prn(AS%isat(im))) .eq. 0) cycle
    do in=im+1,AS%now
      if(pointer_string(ACF%fcbnprn, ACF%fcbprn, ACF%prn(AS%isat(in))) .eq. 0) cycle
      if(ACF%prn(AS%isat(im))(1:1) .ne. ACF%prn(AS%isat(in))(1:1)) cycle
      if(ACF%prn(AS%isat(im))(1:1).eq.'C') then
        read(ACF%prn(AS%isat(im))(2:3),'(i2)') j
        read(ACF%prn(AS%isat(in))(2:3),'(i2)') k
        if(j.lt.18.and.k.gt.18.or.j.gt.18.and.k.lt.18) cycle
      endif
      j = max(AS%iepc(1, im), AS%iepc(1, in))
      k = min(AS%iepc(2, im), AS%iepc(2, in))
      if((k - j)*ACF%dintv .le. ACF%minsec_common) cycle
      AS%nsd = AS%nsd + 1
      if(AS%nsd .gt. MAXSD_SIT) then
        write (*, '(a)') '***ERROR(define_sat_pairs): too many single difference ambiguities'
        call exit(1)
      endif
      ASD(AS%nsd)%id = 2                         ! wide-lane not fixed
      ASD(AS%nsd)%ipt(1) = im                    ! pointer to un-differenced ambiguity in the station
      ASD(AS%nsd)%ipt(2) = in
      ASD(AS%nsd)%rwl = AS%xrwl(im) - AS%xrwl(in)
      ASD(AS%nsd)%swl = dsqrt(AS%xswl(im)**2 + AS%xswl(in)**2)
      ASD(AS%nsd)%iepc(1) = j
      ASD(AS%nsd)%iepc(2) = k
      if(ACF%prn(AS%isat(im))(1:1).ne.'C') then
        ASD(AS%nsd)%sys = ACF%prn(AS%isat(im))(1:1)
        if(ASD(AS%nsd)%sys.eq.'G') AS%nsd_G = AS%nsd_G+1
        if(ASD(AS%nsd)%sys.eq.'E') AS%nsd_E = AS%nsd_E+1
        if(ASD(AS%nsd)%sys.eq.'J') AS%nsd_J = AS%nsd_J+1
      else
        read(ACF%prn(AS%isat(im))(2:3),'(i2)') j
        if(j.lt.18) then
          ASD(AS%nsd)%sys = 'C'
          AS%nsd_C = AS%nsd_C+1
        else
          ASD(AS%nsd)%sys = '3'
          AS%nsd_3 = AS%nsd_3+1
        endif
      endif
    enddo
  enddo

  return
end
