!
!! lsq_add_newamb.f90
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
!! purpose   : add the new ambiguities to the estimator. Newly occured or detected
!!             is indicated by OB%flag.ne.0. Parameter table must be updated and
!!             so do the estimator.
!! parameter :
!!    input  : jd,sod -- time tag
!!             LCF    -- LSQ struct
!!             OB     -- observation struct
!!             NM,PM  -- normal matrix & PAR table
!
subroutine lsq_add_newamb(jd, sod, name, LCF, OB, NM, PM)
  implicit none
  include '../header/const.h'
  include '../header/rnxobs.h'
  include 'lsqcfg.h'
  include 'lsq.h'

  integer*4 jd
  real*8 sod
  character*4 name
  type(lsqcfg) LCF
  type(rnxobr) OB
  type(norm) NM
  type(prmt) PM(1:*)
!
!! local
  integer*4 i, j, iy, imon, id, ih, im, isat, ipar, nbias, ierr
  real*8 sec, ambini, rwini
!
!! function called
  integer*4 pointer_string

  ipar = pointer_string(OB%npar, OB%pname, 'AMBC')

  nbias = 0
  do isat = 1, LCF%nprn
    if (LCF%prn(isat) .eq. '' .or. OB%obs(isat, 1) .eq. 0.d0) cycle
    if (OB%flag(isat) .eq. 1) then
!
!! check whether inserted
      i = OB%ltog(ipar, isat)
      if (i .ne. 0 .and. dabs(PM(i)%ptbeg - OB%lifamb(isat, 1))*86400.d0 .lt. MAXWND .and. &
          dabs(PM(i)%ptend - OB%lifamb(isat, 2))*86400.d0 .lt. MAXWND) cycle
!
!! add new ambiguity
      nbias = nbias + 1
      if (NM%imtx + nbias + 1 .gt. NM%nmtx) then
        write (*, '(a,i5,f15.7)') '***ERROR(lsq_add_newamb): normal matrix too small ', jd, sod
        call exit(1)
      endif
      if (NM%ipm + nbias .gt. NM%npm) then
        write (*, '(a,i5,f15.7)') '***ERROR(lsq_add_newamb): parameter array too small ', jd, sod
        call exit(1)
      endif
      PM(NM%ipm + nbias)%pname = 'AMBC'
      PM(NM%ipm + nbias)%ptype = 'S'
      PM(NM%ipm + nbias)%ipt = NM%imtx + nbias
      PM(NM%ipm + nbias)%iobs = 0
      PM(NM%ipm + nbias)%pcode(1) = 1
      PM(NM%ipm + nbias)%pcode(2) = isat
      PM(NM%ipm + nbias)%ptime(1:2) = 0.d0
      PM(NM%ipm + nbias)%xini = 0.d0
      PM(NM%ipm + nbias)%xcor = 0.d0
      PM(NM%ipm + nbias)%xrms = 0.d0
      PM(NM%ipm + nbias)%xrwl = 0.d0
      PM(NM%ipm + nbias)%xswl = 0.d0
      PM(NM%ipm + nbias)%mele = 0.d0       ! mean elevation
      PM(NM%ipm + nbias)%rw = 0.d0         ! weight of widelane
      PM(NM%ipm + nbias)%zw = 0.d0         ! initial wide-lane value
      OB%ltog(ipar, isat) = NM%ipm + nbias
      NM%iptp(NM%imtx + nbias) = NM%ipm + nbias
!
!! set living time for ambiguity parameters
      if (OB%lifamb(isat, 2) .ne. 0.d0) then
        if ((OB%lifamb(isat, 1) - LCF%jd0)*86400.d0 - LCF%sod0 .lt. 0.d0) OB%lifamb(isat, 1) = LCF%jd0 + LCF%sod0/86400.d0
        if ((OB%lifamb(isat, 2) - LCF%jd1)*86400.d0 - LCF%sod1 .gt. 0.d0) OB%lifamb(isat, 2) = LCF%jd1 + LCF%sod1/86400.d0
        if (OB%lifamb(isat, 1) .gt. OB%lifamb(isat, 2)) then
          write (*, '(a,2f12.5)') '***ERROR(lsq_add_newamb): active time ', OB%lifamb(isat, 1), OB%lifamb(isat, 2)
          call exit(1)
        endif
        PM(NM%ipm + nbias)%ptbeg = OB%lifamb(isat, 1)
        PM(NM%ipm + nbias)%ptend = OB%lifamb(isat, 2)
      else
        call mjd2date(jd, sod, iy, imon, id, ih, im, sec)
        write (*, '(a,a3,1x,a,i5,4i3,f11.7)') '***ERROR(lsq_add_newamb): living ', OB%prn(isat), &
                                                  name//' at', iy, imon, id, ih, im, sec
        call exit(1)
      endif
    endif
  enddo
  if (nbias .eq. 0) return
!
!! shift the rightsite to the last column
  do i = 1, NM%imtx + nbias
    if (i .gt. NM%imtx) then
      NM%norx(i, NM%imtx + 1 + nbias) = 0.d0
    else
      NM%norx(i, NM%imtx + 1 + nbias) = NM%norx(i, NM%imtx + 1)
    endif
  enddo
!
!! clean the part for the new parameters
  do j = 1, nbias
    do i = 1, NM%imtx + j
      if (i .lt. NM%imtx + j) then
        NM%norx(i, NM%imtx + j) = 0.d0
      else if (i .eq. NM%imtx + j) then
        NM%norx(i, i) = 1.d-8
      endif
    enddo
  enddo

  NM%ipm = NM%ipm + nbias
  NM%imtx = NM%imtx + nbias
  NM%ns = NM%ns + nbias

  return
end
