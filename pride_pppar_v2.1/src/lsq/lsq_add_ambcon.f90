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
subroutine lsq_add_ambcon(jd, sod, LCF, SITE, NM, PM, OB)
  implicit none
  include '../header/const.h'
  include '../header/station.h'
  include '../header/rnxobs.h'
  include 'lsqcfg.h'
  include 'lsq.h'

  integer*4 jd
  real*8 sod
  type(lsqcfg) LCF
  type(station) SITE
  type(norm) NM
  type(prmt) PM(1:*)
  type(rnxobr) OB
!
!! local
  logical*1 lfirst, lexist
  integer*4 i, j, ir, ic, lfn, jdw, jdx, iy(2), imon(2), id(2), ih(2), im(2), ipar, ix(2), ierr
  character*3 isat,jsat
  real*8 is(2), iwl, rnl, sodw, sodx, ini(4), op(4)
  character*2 ctyp
  character*4 cname, dname
  character*256 line
!
!! function called
  integer*4 get_valid_unit, modified_julday
  real*8 timdif

  real*8 bc_G,c_G(2)
  real*8 bc_E,c_E(2)
  real*8 bc_C,c_C(2)
  real*8 bc_J,c_J(2)
  character*3 tmp_prn

  data lfirst, lexist/.true., .false./, op/1.d0, -1.d0, -1.d0, 1.d0/
  save lfirst, lexist, lfn, c_G, c_E, c_C, c_J, ctyp

  if (lfirst) then
    lfirst = .false.
    inquire (file=LCF%flncon, exist=lexist)
    if (lexist) then
      write (*, '(a)') '%%%MESSAGE(lsq_add_ambcon): ambiguity constraint found '
      lfn = get_valid_unit(10)
      open (lfn, file=LCF%flncon, status='old', iostat=ierr)
      if (ierr .ne. 0) then
        write (*, '(2a)') '***ERROR(lsq_add_ambcon): open file ', trim(LCF%flncon)
        call exit(1)
      endif
      ! G
      bc_G=FREQ1_G/FREQ2_G
      c_G(1)=bc_G/(bc_G**2-1.d0)
      c_G(2)=FREQ1_G/(FREQ1_G+FREQ2_G)
      ! E
      bc_E=FREQ1_E/FREQ2_E
      c_E(1)=bc_E/(bc_E**2-1.d0)
      c_E(2)=FREQ1_E/(FREQ1_E+FREQ2_E)
      ! C
      bc_C=FREQ1_C/FREQ2_C
      c_C(1)=bc_C/(bc_C**2-1.d0)
      c_C(2)=FREQ1_C/(FREQ1_C+FREQ2_C)
      ! J
      bc_J=FREQ1_J/FREQ2_J
      c_J(1)=bc_J/(bc_J**2-1.d0)
      c_J(2)=FREQ1_J/(FREQ1_J+FREQ2_J)
!
!! read header
      read (lfn, '(a)') line
      do while (index(line, 'END OF HEADER') .eq. 0)
        if (index(line, 'TYPE OF CONSTRAINT') .ne. 0) then
          ctyp = line(5:6)
        endif
        read (lfn, '(a)') line
      enddo
    endif
  endif
  if (.not. lexist) return

10 read (lfn, '(a)', end=100, err=100) line
  if (ctyp .eq. 'SD') then
    dname = ' '
    read (line, '(a4,1x,2(a3,1x),2(i4,4i3,f10.6,1x),2f13.0)', err=200) cname, isat, jsat, &
      iy(1), imon(1), id(1), ih(1), im(1), is(1), iy(2), imon(2), id(2), ih(2), im(2), is(2), iwl, rnl
  endif
!
!! if station missing ...
  j = 0
  if (SITE%name .eq. cname) j = j + 1
  if (ctyp .eq. 'SD' .and. j .lt. 1) goto 10
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
    if (SITE%name .ne. cname) cycle
    if (LCF%prn(PM(ipar)%psat) .ne. isat .and. LCF%prn(PM(ipar)%psat) .ne. jsat) cycle
    if ((PM(ipar)%ptbeg - jdw)*86400.d0 - sodw .gt. MAXWND .or. (PM(ipar)%ptend - jdx)*86400.d0 - sodx .lt. -MAXWND) cycle
    if (SITE%name .eq. cname .and. LCF%prn(PM(ipar)%psat) .eq. isat) then
      tmp_prn=isat
      ix(1) = PM(ipar)%ipt
      ini(1) = PM(ipar)%xini
    else if (SITE%name .eq. cname .and. LCF%prn(PM(ipar)%psat) .eq. jsat) then
      tmp_prn=jsat
      ix(2) = PM(ipar)%ipt
      ini(2) = PM(ipar)%xini
    endif
    if (ctyp .eq. 'SD' .and. ix(1) .ne. 0 .and. ix(2) .ne. 0) exit
  enddo
  if (ix(1) .eq. 0 .or. ix(2) .eq. 0) then
    write (*, '(a,a4,1x,a4,2a3,2(1x,i4,4i3,f10.6))') '***ERROR(lsq_add_ambcon): ambiguity not found ', cname, dname, &
      isat, jsat, iy(1), imon(1), id(1), ih(1), im(1), is(1), iy(2), imon(2), id(2), ih(2), im(2), is(2)
    call exit(1)
  endif
!
!! add constraint to normal equation
  if (ctyp .eq. 'SD') then
    if(tmp_prn(1:1) .eq. 'G')then
      bc_G=c_G(1)*iwl/(FREQ1_G/VLIGHT)+c_G(2)*rnl/(FREQ1_G/VLIGHT)-ini(1)+ini(2)
    elseif(tmp_prn(1:1) .eq. 'E')then
      bc_E=c_E(1)*iwl/(FREQ1_E/VLIGHT)+c_E(2)*rnl/(FREQ1_E/VLIGHT)-ini(1)+ini(2)
    elseif(tmp_prn(1:1) .eq. 'C')then
      bc_C=c_C(1)*iwl/(FREQ1_C/VLIGHT)+c_C(2)*rnl/(FREQ1_C/VLIGHT)-ini(1)+ini(2)
    elseif(tmp_prn(1:1) .eq. 'J')then
      bc_J=c_J(1)*iwl/(FREQ1_J/VLIGHT)+c_J(2)*rnl/(FREQ1_J/VLIGHT)-ini(1)+ini(2)
    end if
  endif
  do i = 1, 2
    if (ix(i) .eq. 0) exit
    do j = i, 2
      if (ix(j) .eq. 0) exit
      ir = min(ix(i), ix(j))
      ic = max(ix(i), ix(j))
      NM%norx(ir, ic) = NM%norx(ir, ic) + op(i)*op(j)*1.d10
    enddo
    ir = ix(i)
    ic = NM%imtx + 1
    if(tmp_prn(1:1) .eq. 'G')then
      NM%norx(ir,ic)=NM%norx(ir,ic)+op(i)*1.d10*bc_G
    elseif(tmp_prn(1:1) .eq. 'E')then
      NM%norx(ir,ic)=NM%norx(ir,ic)+op(i)*1.d10*bc_E
    elseif(tmp_prn(1:1) .eq. 'C')then
      NM%norx(ir,ic)=NM%norx(ir,ic)+op(i)*1.d10*bc_C
    elseif(tmp_prn(1:1) .eq. 'J')then
      NM%norx(ir,ic)=NM%norx(ir,ic)+op(i)*1.d10*bc_J
    end if
  enddo
  if(tmp_prn(1:1) .eq. 'G')then
    NM%ltpl=NM%ltpl+1.d10*bc_G**2
    LCF%nconG = LCF%nconG + 1
  elseif(tmp_prn(1:1) .eq. 'E')then
    NM%ltpl=NM%ltpl+1.d10*bc_E**2
    LCF%nconE = LCF%nconE + 1
  elseif(tmp_prn(1:1) .eq. 'C')then
    NM%ltpl=NM%ltpl+1.d10*bc_C**2
    read(tmp_prn,'(1x,i2)') i
    if(i.le.18) then
      LCF%nconC2 = LCF%nconC2 + 1
    else
      LCF%nconC3 = LCF%nconC3 + 1
    endif
  elseif(tmp_prn(1:1) .eq. 'J')then
    NM%ltpl=NM%ltpl+1.d10*bc_J**2
    LCF%nconJ = LCF%nconJ + 1
  end if
  NM%nobs = NM%nobs + 1

  goto 10

100 return
200 write (*, '(a,/,a)') '***ERROR(lsq_add_ambcon): read file ', trim(line)
  call exit(1)
end
