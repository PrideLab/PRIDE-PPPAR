!
!! lsq_add_ambcon.f90
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
!! Contributor: Maorong Ge, Jianghui Geng, Songfeng Yang, Jihang Lin
!!
!!
!!
!! purpose  : add ambiguity resolution constraint
!! parameter:
!!    input : jd,sod -- epoch time
!!            LCF    -- lsq control struct
!!            SITE   -- station struct
!!    output: NM,PM  -- normal matrix & parameter array
!
subroutine lsq_add_ambcon(lfncid, lfnobs, jd, sod, LCF, SITE, NM, PM, OB)
  implicit none
  include '../header/const.h'
  include '../header/station.h'
  include '../header/rnxobs.h'
  include 'lsqcfg.h'
  include 'lsq.h'

! common
  integer*4     idxfrq(MAXSYS, 2)
  common        idxfrq
! parameter
  integer*4     jd, lfncid, lfnobs
  real*8        sod
  type(lsqcfg)  LCF
  type(station) SITE
  type(norm)    NM
  type(prmt)    PM(1:*)
  type(rnxobr)  OB
! local
  logical*1     lfirst, lexist
  integer*4     i0, i, j, ipar, ierr, lfn
  integer*4     jdw, jdx, iy(2), imon(2), id(2), ih(2), im(2)
  integer*4     ir, ic, ix(2), ip(2)
  real*8        is(2), iwl, rnl, sodw, sodx, ini(2), bc, op(2), weight
  real*8        f1(MAXSYS), f2(MAXSYS), rg(MAXSYS), c1(MAXSYS), c2(MAXSYS)
  character*3   isat, jsat
  character*2   ctyp
  character*4   cname
  character*256 line
! function called
  integer*4     get_valid_unit
  integer*4     pointer_string
  integer*4     modified_julday
  real*8        timdif

  data lfirst, lexist/.true., .false./, op/1.d0, -1.d0/, weight/1.d10/
  save lfirst, lexist, lfn, f1, f2, rg, c1 ,c2, ctyp, weight

  if (lfirst) then
    lfirst = .false.
    inquire (file=LCF%flncon, exist=lexist)
    if (lexist) then
      lfn = get_valid_unit(10)
      open (lfn, file=LCF%flncon, status='old', iostat=ierr)
      if (ierr .ne. 0) then
        write (*, '(2a)') '***ERROR(lsq_add_ambcon): open file ', trim(LCF%flncon)
        call exit(1)
      end if
!
!! read header
      read (lfn, '(a)') line
      read (line(41:), '(a4)') cname
      do while (index(line, 'END OF HEADER') .eq. 0)
        if (index(line, 'TYPE OF CONSTRAINT') .ne. 0) then
          ctyp = line(5:6)
        end if
        read (lfn, '(a)') line
      end do
      if (cname .ne. SITE%name .or. ctyp .ne. 'SD') then
        lexist = .false.
        close (lfn)
        return
      end if
      write (*, '(a)') '%%%MESSAGE(lsq_add_ambcon): ambiguity constraint found '
      do i0 = 1, MAXSYS
      !! frequency
        f1(i0) = FREQ_SYS(idxfrq(i0, 1), i0)
        if (f1(i0) .eq. 0.d0) goto 250
        f2(i0) = FREQ_SYS(idxfrq(i0, 2), i0)
        if (f2(i0) .eq. 0.d0) goto 250
      !! ratio and coeff
        rg(i0) = f1(i0)/f2(i0)
        c1(i0) = rg(i0)/(rg(i0)**2 - 1.d0)
        c2(i0) = 1.d0/(1.d0 + 1.d0/rg(i0))
      end do
    end if
  end if

  if (.not. lexist) return

10 continue
  read (lfn, '(a)', end=100, err=100) line
  if (line(1:1) .ne. ' ') goto 10
  read (line, '(1x,2(a3,1x),2(i4,4i3,f10.6,1x),2f13.0)', err=200) isat, jsat, &
    iy(1), imon(1), id(1), ih(1), im(1), is(1), iy(2), imon(2), id(2), ih(2), im(2), is(2), iwl, rnl
  if (isat(1:1) .eq. 'R') goto 10
  i0 = index(GNSS_PRIO, isat(1:1))
!
!! time comaprison
  jdw = modified_julday(id(1), imon(1), iy(1))
  sodw = ih(1)*3600.d0 + im(1)*60.d0 + is(1)
  jdx = modified_julday(id(2), imon(2), iy(2))
  sodx = ih(2)*3600.d0 + im(2)*60.d0 + is(2)
!
!! if time is outside required duration
  if (timdif(LCF%jd0, LCF%sod0, jdw, sodw) .gt. 0.d0) then
    jdw = LCF%jd0
    sodw = LCF%sod0
  end if
  if (timdif(jdx, sodx, LCF%jd1, LCF%sod1) .gt. 0.d0) then
    jdx = LCF%jd1
    sodx = LCF%sod1
  end if
  if (timdif(jdx, sodx, jdw, sodw) .le. 0.d0) goto 10
!
!! check whether add new constraints
  if (timdif(jdx, sodx, jd, sod) .gt. MAXWND) then
    backspace lfn
    return
  end if
!
!! search corresponding one-way ambiguities in normal euation
  ix = 0
  do i = NM%nc + NM%np + 1, NM%imtx
    ipar = NM%iptp(i)
    if (PM(ipar)%pname(1:4) .ne. 'AMBC') cycle
    if (LCF%prn(PM(ipar)%psat) .ne. isat .and. LCF%prn(PM(ipar)%psat) .ne. jsat) cycle
    if ((PM(ipar)%ptbeg - jdw)*86400.d0 - sodw .gt. MAXWND .or. (PM(ipar)%ptend - jdx)*86400.d0 - sodx .lt. -MAXWND) cycle
    if (LCF%prn(PM(ipar)%psat) .eq. isat) then
      ix(1) = PM(ipar)%ipt
      ip(1) = ipar
      ini(1) = PM(ipar)%xini
    else if (LCF%prn(PM(ipar)%psat) .eq. jsat) then
      ix(2) = PM(ipar)%ipt
      ip(2) = ipar
      ini(2) = PM(ipar)%xini
    end if
    if (ix(1) .ne. 0 .and. ix(2) .ne. 0) exit
  end do
  if (ix(1) .eq. 0 .or. ix(2) .eq. 0) then
    write (*, '(a,1x,2(a3,1x),2(i4,4i3,f10.6,1x))') '***ERROR(lsq_add_ambcon): ambiguity not found ', &
      isat, jsat, iy(1), imon(1), id(1), ih(1), im(1), is(1), iy(2), imon(2), id(2), ih(2), im(2), is(2)
    call exit(1)
  end if
!
!! add constraint to normal equation
  bc = iwl*c1(i0)/f1(i0)*VLIGHT - ini(1) + rnl*c2(i0)/f1(i0)*VLIGHT + ini(2)
  do i = 1, 2
    if (ix(i) .eq. 0) exit
    do j = i, 2
      if (ix(j) .eq. 0) exit
      ir = min(ix(i), ix(j))
      ic = max(ix(i), ix(j))
      NM%norx(ir, ic) = NM%norx(ir, ic) + op(i)*weight*op(j)
    end do
    ir = ix(i)
    ic = NM%imtx + 1
    NM%norx(ir, ic) = NM%norx(ir, ic) + op(i)*weight*bc
  end do
  if (lfnobs .ne. 0) then
    write (lfncid) 'cn'
    i = pointer_string(LCF%nprn, LCF%prn, isat)
    j = pointer_string(LCF%nprn, LCF%prn, jsat)
    write (lfnobs) jdw, sodw, jdx, sodx, i, j, ip(1:2), bc, weight
  end if
  if (isat(1:1) .eq. 'G') then
    LCF%nconG = LCF%nconG + 1
  elseif (isat(1:1) .eq. 'E') then
    LCF%nconE = LCF%nconE + 1
  elseif (isat(1:1) .eq. 'C') then
    read (isat, '(1x,i2)') i
    if (i .le. 18) then
      LCF%nconC2 = LCF%nconC2 + 1
    else
      LCF%nconC3 = LCF%nconC3 + 1
    end if
  elseif (isat(1:1) .eq. 'J') then
    LCF%nconJ = LCF%nconJ + 1
  end if
  NM%ltpl = NM%ltpl + weight*bc**2
  NM%nobs = NM%nobs + 1
  goto 10

100 continue
  return

200 continue
  write (*, '(a,/,a)') '***ERROR(lsq_add_ambcon): read file ', trim(line)
  call exit(1)

250 continue
  write (*, '(a,5(1x,a,2i1))') '***ERROR(read_rinex_file): invalid frequency number:', &
    'G', idxfrq(index(GNSS_PRIO, 'G'), 1:2), &
    'R', idxfrq(index(GNSS_PRIO, 'R'), 1:2), &
    'E', idxfrq(index(GNSS_PRIO, 'E'), 1:2), &
    'C', idxfrq(index(GNSS_PRIO, 'C'), 1:2), &
    'J', idxfrq(index(GNSS_PRIO, 'J'), 1:2)
  call exit(1)
end subroutine
