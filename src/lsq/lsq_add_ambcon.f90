!
!! lsq_add_ambcon.f90
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

  integer*4 jd, lfncid, lfnobs
  real*8 sod
  type(lsqcfg) LCF
  type(station) SITE
  type(norm) NM
  type(prmt) PM(1:*)
  type(rnxobr) OB
!
!! local
  logical*1 lfirst, lexist
  integer*4 i, j, ir, ic, lfn, jdw, jdx, iy(2), imon(2), id(2), ih(2), im(2), ipar, ix(2), ip(2), ierr
  character*3 isat,jsat
  real*8 is(2), iwl, rnl, sodw, sodx, ini(2), op(2), weight
  real*8 bc, c_G(2), c_E(2), c_C(2), c_J(2)
  character*2 ctyp
  character*4 cname
  character*256 line
!
!! function called
  integer*4 get_valid_unit, pointer_string, modified_julday
  real*8 timdif

  data lfirst, lexist/.true., .false./, op/1.d0, -1.d0/, weight/1.d10/
  save lfirst, lexist, lfn, c_G, c_E, c_C, c_J, ctyp, weight

  if (lfirst) then
    lfirst = .false.
    inquire (file=LCF%flncon, exist=lexist)
    if (lexist) then
      lfn = get_valid_unit(10)
      open (lfn, file=LCF%flncon, status='old', iostat=ierr)
      if (ierr .ne. 0) then
        write (*, '(2a)') '***ERROR(lsq_add_ambcon): open file ', trim(LCF%flncon)
        call exit(1)
      endif
!
!! read header
      read (lfn, '(a)') line
      read (line(41:),'(a4)') cname
      do while (index(line, 'END OF HEADER') .eq. 0)
        if (index(line, 'TYPE OF CONSTRAINT') .ne. 0) then
          ctyp = line(5:6)
        endif
        read (lfn, '(a)') line
      enddo
      if (cname.ne.SITE%name .or. ctyp.ne.'SD') then
        lexist=.false.
        close(lfn)
        return
      endif
      write (*, '(a)') '%%%MESSAGE(lsq_add_ambcon): ambiguity constraint found '
      ! G
      bc=FREQ1_G/FREQ2_G
      c_G(1)=bc/(bc**2-1.d0)
      c_G(2)=FREQ1_G/(FREQ1_G+FREQ2_G)
      ! E
      bc=FREQ1_E/FREQ2_E
      c_E(1)=bc/(bc**2-1.d0)
      c_E(2)=FREQ1_E/(FREQ1_E+FREQ2_E)
      ! C
      bc=FREQ1_C/FREQ2_C
      c_C(1)=bc/(bc**2-1.d0)
      c_C(2)=FREQ1_C/(FREQ1_C+FREQ2_C)
      ! J
      bc=FREQ1_J/FREQ2_J
      c_J(1)=bc/(bc**2-1.d0)
      c_J(2)=FREQ1_J/(FREQ1_J+FREQ2_J)
    endif
  endif
  if (.not. lexist) return

10 read (lfn, '(a)', end=100, err=100) line
  if(line(1:1).ne.' ') goto 10
  read (line, '(1x,2(a3,1x),2(i4,4i3,f10.6,1x),2f13.0)', err=200) isat, jsat, &
    iy(1), imon(1), id(1), ih(1), im(1), is(1), iy(2), imon(2), id(2), ih(2), im(2), is(2), iwl, rnl
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
  endif
  if (timdif(jdx, sodx, LCF%jd1, LCF%sod1) .gt. 0.d0) then
    jdx = LCF%jd1
    sodx = LCF%sod1
  endif
  if (timdif(jdx, sodx, jdw, sodw) .le. 0.d0) goto 10
!
!! check whether add new constraints
  if (timdif(jdx, sodx, jd, sod) .gt. MAXWND) then
    backspace lfn
    return
  endif
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
    endif
    if (ix(1) .ne. 0 .and. ix(2) .ne. 0) exit
  enddo
  if (ix(1) .eq. 0 .or. ix(2) .eq. 0) then
    write (*, '(a,1x,2a3,2(1x,i4,4i3,f10.6))') '***ERROR(lsq_add_ambcon): ambiguity not found ', &
      isat, jsat, iy(1), imon(1), id(1), ih(1), im(1), is(1), iy(2), imon(2), id(2), ih(2), im(2), is(2)
    call exit(1)
  endif
!
!! add constraint to normal equation
  if(isat(1:1) .eq. 'G')then
    bc=c_G(1)*iwl/(FREQ1_G/VLIGHT)+c_G(2)*rnl/(FREQ1_G/VLIGHT)-ini(1)+ini(2)
  elseif(isat(1:1) .eq. 'E')then
    bc=c_E(1)*iwl/(FREQ1_E/VLIGHT)+c_E(2)*rnl/(FREQ1_E/VLIGHT)-ini(1)+ini(2)
  elseif(isat(1:1) .eq. 'C')then
    bc=c_C(1)*iwl/(FREQ1_C/VLIGHT)+c_C(2)*rnl/(FREQ1_C/VLIGHT)-ini(1)+ini(2)
  elseif(isat(1:1) .eq. 'J')then
    bc=c_J(1)*iwl/(FREQ1_J/VLIGHT)+c_J(2)*rnl/(FREQ1_J/VLIGHT)-ini(1)+ini(2)
  end if
  do i = 1, 2
    if (ix(i) .eq. 0) exit
    do j = i, 2
      if (ix(j) .eq. 0) exit
      ir = min(ix(i), ix(j))
      ic = max(ix(i), ix(j))
      NM%norx(ir, ic) = NM%norx(ir, ic) + op(i)*op(j)*weight
    enddo
    ir = ix(i)
    ic = NM%imtx + 1
    if(isat(1:1) .eq. 'G')then
      NM%norx(ir,ic)=NM%norx(ir,ic)+op(i)*weight*bc
    elseif(isat(1:1) .eq. 'E')then
      NM%norx(ir,ic)=NM%norx(ir,ic)+op(i)*weight*bc
    elseif(isat(1:1) .eq. 'C')then
      NM%norx(ir,ic)=NM%norx(ir,ic)+op(i)*weight*bc
    elseif(isat(1:1) .eq. 'J')then
      NM%norx(ir,ic)=NM%norx(ir,ic)+op(i)*weight*bc
    end if
  enddo

  if (lfnobs.ne.0) then
    write(lfncid) 'cn'
    i=pointer_string(LCF%nprn,LCF%prn,isat)
    j=pointer_string(LCF%nprn,LCF%prn,jsat)
    write(lfnobs) jdw,sodw,jdx,sodx,i,j,ip(1:2),bc,weight
  endif

  if(isat(1:1) .eq. 'G')then
    NM%ltpl=NM%ltpl+weight*bc**2
    LCF%nconG = LCF%nconG + 1
  elseif(isat(1:1) .eq. 'E')then
    NM%ltpl=NM%ltpl+weight*bc**2
    LCF%nconE = LCF%nconE + 1
  elseif(isat(1:1) .eq. 'C')then
    NM%ltpl=NM%ltpl+weight*bc**2
    read(isat,'(1x,i2)') i
    if(i.le.18) then
      LCF%nconC2 = LCF%nconC2 + 1
    else
      LCF%nconC3 = LCF%nconC3 + 1
    endif
  elseif(isat(1:1) .eq. 'J')then
    NM%ltpl=NM%ltpl+weight*bc**2
    LCF%nconJ = LCF%nconJ + 1
  end if
  NM%nobs = NM%nobs + 1

  goto 10

100 return
200 write (*, '(a,/,a)') '***ERROR(lsq_add_ambcon): read file ', trim(line)
  call exit(1)
end
