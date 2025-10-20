!
!! write_ambcon.f90
!!
!!    Copyright (C) 2018 by Wuhan University
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
!! Contributor: Maorong Ge, Jianghui Geng, Songfeng Yang
!! 
!!
!!
!! purpose  : write ambiguity constraint
!! parameter:
!!    input : indp -- # of independent SD ambiguities
!!            SD   -- single-difference ambiguities
!!            ACF  -- fractional part of initial phases
!!            AS   -- ambiguity station
!!    output:
!
subroutine write_ambcon(indp, SD, ACF, AS)
  implicit none
  include '../header/const.h'
  include 'ambssd.h'
  include 'arscfg.h'
  include 'ambsit.h'

  integer*4 indp
  type(arscfg) ACF
  type(ambsit) AS
  type(ambssd) SD(1:*)
!
!! local
  integer*4 i, j, lfn, jd, iy(2), imon(2), id(2), ih(2), im(2)
  real*8 is(2), sod
  type(ambssd) XSD
  integer*4 isat, jsat
!
!! function used
  integer*4 get_valid_unit
!
!! re-arrange SD
  do i = 1, indp - 1
    do j = i + 1, indp
      if (SD(i)%iepc(2) .gt. SD(j)%iepc(2)) then
        XSD = SD(i)
        SD(i) = SD(j)
        SD(j) = XSD
      endif
    enddo
  enddo
!
!! write constraint header
  lfn = get_valid_unit(10)
  open (lfn, file=ACF%flncon)
  write (lfn, '(a38,2x,a4,16x,a)') 'Single-Difference Ambiguity Constraint', AS%name, 'COMMENT'
  write (lfn, '(4x,a2,54x,a)') 'SD', 'TYPE OF CONSTRAINT'
  if(ACF%ntot_G.gt.0) write (lfn, '(3i6,42x,a)') ACF%ntot_G,ACF%nwlfx_G,ACF%nwnfx_G,'G   AMB FIXING (T/W/N)'
  if(ACF%ntot_E.gt.0) write (lfn, '(3i6,42x,a)') ACF%ntot_E,ACF%nwlfx_E,ACF%nwnfx_E,'E   AMB FIXING (T/W/N)'
  if(ACF%ntot_C.gt.0) write (lfn, '(3i6,42x,a)') ACF%ntot_C,ACF%nwlfx_C,ACF%nwnfx_C,'C2  AMB FIXING (T/W/N)'
  if(ACF%ntot_3.gt.0) write (lfn, '(3i6,42x,a)') ACF%ntot_3,ACF%nwlfx_3,ACF%nwnfx_3,'C3  AMB FIXING (T/W/N)'
  if(ACF%ntot_J.gt.0) write (lfn, '(3i6,42x,a)') ACF%ntot_J,ACF%nwlfx_J,ACF%nwnfx_J,'J   AMB FIXING (T/W/N)'
  if(ACF%lsearch) then
    write (lfn, '(3i6,4x,a8,30x,a)') ACF%ntot_ind,ACF%nwl_ind,ACF%nwn_ind,'LAMBDA  ','IND AMB FIXING (T/W/N)'
    if(ACF%features(1) .gt. 1.1) then
        write(lfn,'(i12, 3f12.4,12x,a)') nint(ACF%features(1)), (ACF%features(i), i=2,4), 'NCAD DISMIN RATIO PENAL'
        write(lfn,'(4f12.4,12x,a)')      (ACF%features(i), i=5,8), 'ADOP WRATIO PRATIO NOBS'
    endif
  else
    write (lfn, '(3i6,4x,a8,30x,a)') ACF%ntot_ind,ACF%nwl_ind,ACF%nwn_ind,'ROUNDING','IND AMB FIXING (T/W/N)'
  endif

  if(indp .ne. 0) then
     write (lfn, '(4x, a6,50x,a)') 'SUCESS','AR RESULT'
  else
     write (lfn, '(4x, a6,50x,a)') 'FAIL  ','AR RESULT'
  endif

  write (lfn, '(60x,a)') 'END OF HEADER'
!
!! write each constraint
  do i = 1, indp
    do j = 1, 2
      call timinc(ACF%jd0, ACF%sod0, ACF%dintv*(SD(i)%iepc(j) - 1), jd, sod)
      call mjd2date(jd, sod, iy(j), imon(j), id(j), ih(j), im(j), is(j))
    enddo
    write (lfn, '(1x,2(a3,1x),2(i4,4i3,f10.6,1x),2i13)') ACF%prn(SD(i)%ipt(1)), ACF%prn(SD(i)%ipt(2)), &
      iy(1), imon(1), id(1), ih(1), im(1), is(1), iy(2), imon(2), id(2), ih(2), im(2), is(2), nint(SD(i)%rwl), &
      nint(SD(i)%rnl)
  enddo
  close (lfn)

  return
end
