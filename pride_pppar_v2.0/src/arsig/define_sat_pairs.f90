!
!! define_sat_pairs.f90
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
!! Contributor: Maorong Ge, Jianghui Geng, Songfeng Yang
!! 
!!
!!
!! purpose  : define satellite pairs per station
!! parameter:
!!    input : FCB -- fractional cycle biases
!!            AS  -- station ambiguity
!!    output: ASD -- single-difference ambiguities
!
subroutine define_sat_pairs(FCB, AS, ASD)
  implicit none
  include '../header/const.h'
  include '../header/difamb.h'
  include 'abfcb.h'
  include 'ambsit.h'

  type(abfcb) FCB
  type(ambsit) AS
  type(difamb) ASD(MAXSD_SIT)
!
!! local
  integer*4 j, k, isat, jsat, im, in, is
  integer*4 iprn_int,jprn_int
!
!! function used
  integer*4 pointer_int

  AS%nsd = 0
  AS%nsd_G = 0
  AS%nsd_R = 0
  AS%nsd_E = 0
  AS%nsd_C = 0
  AS%nsd_3 = 0
  AS%nsd_J = 0
  ASD%sys=''
  do isat = 1, FCB%nprn - 1
    im = 0
    do while (im .lt. AS%now)
      is = pointer_int(AS%now - im, AS%isat(im + 1:), isat)
      im = is + im
      if (is .eq. 0) exit
      do jsat = isat + 1, FCB%nprn
        in = 0
        do while (in .lt. AS%now)
          is = pointer_int(AS%now - in, AS%isat(in + 1:), jsat)
          in = is + in
          if (is .eq. 0) exit
          if(FCB%prn(isat)(1:1) .eq. 'G' .and. FCB%prn(jsat)(1:1) .eq. 'G')then
            j = max(AS%iepc(1, im), AS%iepc(1, in))
            k = min(AS%iepc(2, im), AS%iepc(2, in))
            if ((k - j)*FCB%dintv .le. FCB%minsec_common) cycle
!            
!! single difference ambiguity between satellites
            AS%nsd = AS%nsd + 1
            AS%nsd_G = AS%nsd_G + 1
            if (AS%nsd .gt. MAXSD_SIT) then
              write (*, '(a,a4)') '***ERROR(define_sat_pairs): too many single difference ambiguities ', AS%name
              call exit(1)
            endif
            ASD(AS%nsd)%id = 2                         ! wide-lane not fixed
            ASD(AS%nsd)%ipt(1) = im                    ! pointer to un-differenced ambiguity in the station
            ASD(AS%nsd)%ipt(2) = in
            ASD(AS%nsd)%ipt(3:4) = 0
            ASD(AS%nsd)%rwl = AS%xrwl(im) - AS%xrwl(in)
            ASD(AS%nsd)%swl = dsqrt(AS%xswl(im)**2 + AS%xswl(in)**2)
            ASD(AS%nsd)%iepc(1) = j
            ASD(AS%nsd)%iepc(2) = k
            ASD(AS%nsd)%rsv = 0.d0
            ASD(AS%nsd)%sys = 'G'
          elseif(FCB%prn(isat)(1:1) .eq. 'E' .and. FCB%prn(jsat)(1:1) .eq. 'E')then
            j = max(AS%iepc(1, im), AS%iepc(1, in))
            k = min(AS%iepc(2, im), AS%iepc(2, in))
            if ((k - j)*FCB%dintv .le. FCB%minsec_common) cycle
!            
!! single difference ambiguity between satellites
            AS%nsd = AS%nsd + 1
            AS%nsd_E = AS%nsd_E + 1
            if (AS%nsd .gt. MAXSD_SIT) then
              write (*, '(a,a4)') '***ERROR(define_sat_pairs): too many single difference ambiguities ', AS%name
              call exit(1)
            endif
            ASD(AS%nsd)%id = 2                         ! wide-lane not fixed
            ASD(AS%nsd)%ipt(1) = im                    ! pointer to un-differenced ambiguity in the station
            ASD(AS%nsd)%ipt(2) = in
            ASD(AS%nsd)%ipt(3:4) = 0
            ASD(AS%nsd)%rwl = AS%xrwl(im) - AS%xrwl(in)
            ASD(AS%nsd)%swl = dsqrt(AS%xswl(im)**2 + AS%xswl(in)**2)
            ASD(AS%nsd)%iepc(1) = j
            ASD(AS%nsd)%iepc(2) = k
            ASD(AS%nsd)%rsv = 0.d0
            ASD(AS%nsd)%sys = 'E'
          elseif(FCB%prn(isat)(1:1) .eq. 'C' .and. FCB%prn(jsat)(1:1) .eq. 'C')then
            read(FCB%prn(isat)(2:3),'(i2)')iprn_int
            read(FCB%prn(jsat)(2:3),'(i2)')jprn_int
            if(iprn_int .le. 17 .and. jprn_int .le. 17)then
              j = max(AS%iepc(1, im), AS%iepc(1, in))
              k = min(AS%iepc(2, im), AS%iepc(2, in))
              if ((k - j)*FCB%dintv .le. FCB%minsec_common) cycle
!            
!! single difference ambiguity between satellites
              AS%nsd = AS%nsd + 1
              AS%nsd_C = AS%nsd_C + 1
              if (AS%nsd .gt. MAXSD_SIT) then
                write (*, '(a,a4)') '***ERROR(define_sat_pairs): too many single difference ambiguities ', AS%name
                call exit(1)
              endif
              ASD(AS%nsd)%id = 2                         ! wide-lane not fixed
              ASD(AS%nsd)%ipt(1) = im                    ! pointer to un-differenced ambiguity in the station
              ASD(AS%nsd)%ipt(2) = in
              ASD(AS%nsd)%ipt(3:4) = 0
              ASD(AS%nsd)%rwl = AS%xrwl(im) - AS%xrwl(in)
              ASD(AS%nsd)%swl = dsqrt(AS%xswl(im)**2 + AS%xswl(in)**2)
              ASD(AS%nsd)%iepc(1) = j
              ASD(AS%nsd)%iepc(2) = k
              ASD(AS%nsd)%rsv = 0.d0
              ASD(AS%nsd)%sys = 'C'
            elseif(iprn_int .ge. 18 .and. jprn_int .ge. 18)then
              j = max(AS%iepc(1, im), AS%iepc(1, in))
              k = min(AS%iepc(2, im), AS%iepc(2, in))
              if ((k - j)*FCB%dintv .le. FCB%minsec_common) cycle
!            
!! single difference ambiguity between satellites
              AS%nsd = AS%nsd + 1
              AS%nsd_3 = AS%nsd_3 + 1
              if (AS%nsd .gt. MAXSD_SIT) then
                write (*, '(a,a4)') '***ERROR(define_sat_pairs): too many single difference ambiguities ', AS%name
                call exit(1)
              endif
              ASD(AS%nsd)%id = 2                         ! wide-lane not fixed
              ASD(AS%nsd)%ipt(1) = im                    ! pointer to un-differenced ambiguity in the station
              ASD(AS%nsd)%ipt(2) = in
              ASD(AS%nsd)%ipt(3:4) = 0
              ASD(AS%nsd)%rwl = AS%xrwl(im) - AS%xrwl(in)
              ASD(AS%nsd)%swl = dsqrt(AS%xswl(im)**2 + AS%xswl(in)**2)
              ASD(AS%nsd)%iepc(1) = j
              ASD(AS%nsd)%iepc(2) = k
              ASD(AS%nsd)%rsv = 0.d0
              ASD(AS%nsd)%sys = '3'
            endif
          elseif(FCB%prn(isat)(1:1) .eq. 'J' .and. FCB%prn(jsat)(1:1) .eq. 'J')then
            j = max(AS%iepc(1, im), AS%iepc(1, in))
            k = min(AS%iepc(2, im), AS%iepc(2, in))
            if ((k - j)*FCB%dintv .le. FCB%minsec_common) cycle
!            
!! single difference ambiguity between satellites
            AS%nsd = AS%nsd + 1
            AS%nsd_J = AS%nsd_J + 1
            if (AS%nsd .gt. MAXSD_SIT) then
              write (*, '(a,a4)') '***ERROR(define_sat_pairs): too many single difference ambiguities ', AS%name
              call exit(1)
            endif
            ASD(AS%nsd)%id = 2                         ! wide-lane not fixed
            ASD(AS%nsd)%ipt(1) = im                    ! pointer to un-differenced ambiguity in the station
            ASD(AS%nsd)%ipt(2) = in
            ASD(AS%nsd)%ipt(3:4) = 0
            ASD(AS%nsd)%rwl = AS%xrwl(im) - AS%xrwl(in)
            ASD(AS%nsd)%swl = dsqrt(AS%xswl(im)**2 + AS%xswl(in)**2)
            ASD(AS%nsd)%iepc(1) = j
            ASD(AS%nsd)%iepc(2) = k
            ASD(AS%nsd)%rsv = 0.d0
            ASD(AS%nsd)%sys = 'J'
          endif
        enddo
      enddo
    enddo
  enddo

  return
end
