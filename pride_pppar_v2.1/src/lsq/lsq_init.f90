!
!! lsq_init.f90
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
!! purpose   : initialization of the least squares estimator
!! parameters:
!!             LCF   -- LSQ control struct
!!             SITE  -- station struct
!!             OB    -- observation struct
!!             NM,PM -- normal matrix & PAR table
!
subroutine lsq_init(LCF, SITE, SAT, OB, NM, PM)
  implicit none
  include '../header/const.h'
  include '../header/station.h'
  include '../header/satellite.h'
  include '../header/rnxobs.h'
  include 'lsqcfg.h'
  include 'lsq.h'

  type(lsqcfg) LCF
  type(station) SITE
  type(satellite) SAT(MAXSAT)
  type(rnxobr) OB
  type(norm) NM
  type(prmt) PM(1:*)
!
!! local
  integer*4 isat, ipar, irow
  integer*4 i, j, k, ic, ip, kpar, minut
  character*4 sitename
  character*3 xyz, prnname
  character*2 htg
  real*8 val

  data xyz, htg/'XYZ', 'CS'/
!
!! initialization
  do ipar = 1, NM%imtx + 1
    if (ipar .ne. NM%imtx + 1) then
      PM(ipar)%ptype = ' '
      PM(ipar)%ipt = 0
      PM(ipar)%iobs = 0
      PM(ipar)%iobs_G = 0
      PM(ipar)%iobs_R = 0
      PM(ipar)%iobs_E = 0
      PM(ipar)%iobs_C = 0
      PM(ipar)%iobs_3 = 0
      PM(ipar)%iobs_J = 0
      PM(ipar)%ptime(1) = LCF%jd0 + LCF%sod0/86400.d0
      PM(ipar)%ptime(2) = LCF%jd1 + LCF%sod1/86400.d0
      PM(ipar)%ptbeg = 0.d0
      PM(ipar)%ptend = 0.d0
      PM(ipar)%map = 0.d0
      PM(ipar)%rw = 0.d0
      PM(ipar)%xcor = 0.d0
      PM(ipar)%xsig = 0.d0
      PM(ipar)%xrms = 0.d0
      PM(ipar)%xrwl = 0.d0      ! Especially Melbourne-Wubbena
      PM(ipar)%xswl = 0.d0
      PM(ipar)%mele = 0.d0
      PM(ipar)%zw = 0.d0        ! initial value of  wide-lane
    endif
    do irow = 1, NM%imtx
      NM%norx(irow, ipar) = 0.d0
    enddo
  enddo

  do ipar = 1, MAXPAR_STA
    do isat = 1, MAXSAT
      OB%ltog(ipar, isat) = 0
    enddo
  enddo
  do isat = 1, MAXSAT
    OB%delay(isat) = 0.d0
    OB%flag(isat) = 0
  enddo
  ic = 0
  ip = 0
  minut = 0
!
!! station depended
  OB%npar = 0
  if (SITE%skd(1:1) .eq. 'S' .or. SITE%skd(1:1) .eq. 'F') then
    do i = 1, 3
      ic = ic + 1
      ipar = ic
      PM(ipar)%pname = 'STAP'//xyz(i:i)
      PM(ipar)%ptype = 'C'
      PM(ipar)%ipt = ipar
      PM(ipar)%psat = 0
      PM(ipar)%xini = SITE%x(i)*1.d3       ! unit: m
      OB%npar = OB%npar + 1
      OB%pname(OB%npar) = PM(ipar)%pname
      OB%ltog(OB%npar, 1) = ipar
      NM%norx(ipar, ipar) = 1.d0/SITE%dx0(i)**2
      NM%iptp(ipar) = ipar
    enddo
  else if (SITE%skd(1:1) .eq. 'K') then
    do i = 1, 3
      ip = ip + 1
      ipar = NM%nc + ip
      PM(ipar)%pname = 'STAP'//xyz(i:i)
      PM(ipar)%ptype = 'P'
      PM(ipar)%ipt = ipar
      PM(ipar)%psat = 0
      PM(ipar)%xini = 1.d3       ! unit: m
      PM(ipar)%ptime(2) = LCF%jd0 + LCF%sod0/86400.d0
      PM(ipar)%map = 0.d0
      PM(ipar)%rw = 1.d0/SITE%dx0(i)
      OB%npar = OB%npar + 1
      OB%pname(OB%npar) = PM(ipar)%pname
      OB%ltog(OB%npar, 1) = ipar
      NM%norx(ipar, ipar) = 1.d0/SITE%dx0(i)**2
      NM%iptp(ipar) = ipar
    enddo
  else
    write (*, '(2a)') '***ERROR(lsq_init): unknown site type ', SITE%name//' '//SITE%skd
    call exit(1)
  endif
!
!! receiver clock offset
! G
  if(LCF%sys(1:1) .eq. 'G')then
    ip = ip + 1
    ipar = NM%nc + ip
    PM(ipar)%pname = 'RECCLK_G'
    PM(ipar)%ptype = 'P'
    PM(ipar)%ipt = ipar
    PM(ipar)%psat = 0
    PM(ipar)%ptime(2) = LCF%jd0 + LCF%sod0/86400.d0
    PM(ipar)%map = 0.d0
    PM(ipar)%rw = 1.d0/SITE%dclk0_G
    PM(ipar)%xini = 0.d0
    OB%npar = OB%npar + 1
    OB%pname(OB%npar) = PM(ipar)%pname
    OB%ltog(OB%npar, 1) = ipar
    NM%norx(ipar, ipar) = 1.d0/SITE%dclk0_G**2
    NM%iptp(ipar) = ipar
  endif
! R
  if(LCF%sys(2:2) .eq. 'R')then
    ip = ip + 1
    ipar = NM%nc + ip
    PM(ipar)%pname = 'RECCLK_R'
    PM(ipar)%ptype = 'P' 
    PM(ipar)%ipt = ipar
    PM(ipar)%psat = 0 
    PM(ipar)%ptime(2) = LCF%jd0 + LCF%sod0/86400.d0
    PM(ipar)%map = 0.d0
    PM(ipar)%rw = 1.d0/SITE%dclk0_R
    PM(ipar)%xini = 0.d0
    OB%npar = OB%npar + 1 
    OB%pname(OB%npar) = PM(ipar)%pname
    OB%ltog(OB%npar, 1) = ipar
    NM%norx(ipar, ipar) = 1.d0/SITE%dclk0_R**2
    NM%iptp(ipar) = ipar
  endif
! E
  if(LCF%sys(3:3) .eq. 'E')then
    ip = ip + 1
    ipar = NM%nc + ip
    PM(ipar)%pname = 'RECCLK_E'
    PM(ipar)%ptype = 'P' 
    PM(ipar)%ipt = ipar
    PM(ipar)%psat = 0 
    PM(ipar)%ptime(2) = LCF%jd0 + LCF%sod0/86400.d0
    PM(ipar)%map = 0.d0
    PM(ipar)%rw = 1.d0/SITE%dclk0_E
    PM(ipar)%xini = 0.d0
    OB%npar = OB%npar + 1 
    OB%pname(OB%npar) = PM(ipar)%pname
    OB%ltog(OB%npar, 1) = ipar
    NM%norx(ipar, ipar) = 1.d0/SITE%dclk0_E**2
    NM%iptp(ipar) = ipar
  endif
! C
  if(LCF%sys(4:4) .eq. 'C')then
    ip = ip + 1
    ipar = NM%nc + ip
    PM(ipar)%pname = 'RECCLK_C'
    PM(ipar)%ptype = 'P' 
    PM(ipar)%ipt = ipar
    PM(ipar)%psat = 0 
    PM(ipar)%ptime(2) = LCF%jd0 + LCF%sod0/86400.d0
    PM(ipar)%map = 0.d0
    PM(ipar)%rw = 1.d0/SITE%dclk0_C
    PM(ipar)%xini = 0.d0
    OB%npar = OB%npar + 1 
    OB%pname(OB%npar) = PM(ipar)%pname
    OB%ltog(OB%npar, 1) = ipar
    NM%norx(ipar, ipar) = 1.d0/SITE%dclk0_C**2
    NM%iptp(ipar) = ipar
  endif
! 3
  if(LCF%sys(5:5) .eq. '3')then
    ip = ip + 1
    ipar = NM%nc + ip
    PM(ipar)%pname = 'RECCLK_3'
    PM(ipar)%ptype = 'P' 
    PM(ipar)%ipt = ipar
    PM(ipar)%psat = 0 
    PM(ipar)%ptime(2) = LCF%jd0 + LCF%sod0/86400.d0
    PM(ipar)%map = 0.d0
    PM(ipar)%rw = 1.d0/SITE%dclk0_3
    PM(ipar)%xini = 0.d0
    OB%npar = OB%npar + 1 
    OB%pname(OB%npar) = PM(ipar)%pname
    OB%ltog(OB%npar, 1) = ipar
    NM%norx(ipar, ipar) = 1.d0/SITE%dclk0_3**2
    NM%iptp(ipar) = ipar
  endif
! J
  if(LCF%sys(6:6) .eq. 'J')then
    ip = ip + 1
    ipar = NM%nc + ip
    PM(ipar)%pname = 'RECCLK_J'
    PM(ipar)%ptype = 'P' 
    PM(ipar)%ipt = ipar
    PM(ipar)%psat = 0 
    PM(ipar)%ptime(2) = LCF%jd0 + LCF%sod0/86400.d0
    PM(ipar)%map = 0.d0
    PM(ipar)%rw = 1.d0/SITE%dclk0_J
    PM(ipar)%xini = 0.d0
    OB%npar = OB%npar + 1 
    OB%pname(OB%npar) = PM(ipar)%pname
    OB%ltog(OB%npar, 1) = ipar
    NM%norx(ipar, ipar) = 1.d0/SITE%dclk0_J**2
    NM%iptp(ipar) = ipar
  endif
!
!! zenith atmospheric delay
  if (LCF%ztdmod(1:3) .eq. 'PWC') then
    k = index(LCF%ztdmod, ':')
    read (LCF%ztdmod(k + 1:), *) minut
    kpar = nint(((LCF%jd1 - LCF%jd0)*1440.d0 + (LCF%sod1 - LCF%sod0)/60.d0)/minut)
  endif
  ip = ip + 1
  ipar = NM%nc + ip
  PM(ipar)%pname = 'ZTD'//trim(LCF%ztdmod)
  PM(ipar)%ptype = 'P'
  PM(ipar)%ipt = ipar
  PM(ipar)%psat = 0
  if (LCF%ztdmod(1:3) .eq. 'STO') then
    PM(ipar)%ptime(2) = LCF%jd0 + LCF%sod0/86400.d0
  else if (LCF%ztdmod(1:3) .eq. 'PWC' .and. kpar .gt. 1) then
    PM(ipar)%ptime(2) = LCF%jd0 + LCF%sod0/86400.d0 + minut/1440.d0
  endif
  PM(ipar)%map = 1.d0
  if (LCF%ztdmod(1:3) .eq. 'STO') then
    PM(ipar)%rw = 1.d0/(SITE%qztd*dsqrt(LCF%dintv/3600.d0))
  else if (LCF%ztdmod(1:3) .eq. 'PWC') then
    PM(ipar)%rw = 1.d0/(SITE%qztd*dsqrt(minut/60.d0))
  endif
  PM(ipar)%xini = 0.d0
  OB%npar = OB%npar + 1
  OB%pname(OB%npar) = PM(ipar)%pname
  OB%ltog(OB%npar, 1) = ipar
  NM%norx(ipar, ipar) = 1.d0/SITE%dztd0**2
  NM%iptp(ipar) = ipar
!
!! horizontal troposphere gradients
  if (LCF%htgmod(1:3) .eq. 'PWC') then
    k = index(LCF%htgmod, ':')
    read (LCF%htgmod(k + 1:), *) minut
    kpar = nint(((LCF%jd1 - LCF%jd0)*1440.d0 + (LCF%sod1 - LCF%sod0)/60.d0)/minut)
  endif
  if (LCF%htgmod(1:3) .ne. 'NON') then
    do i = 1, 2
      ip = ip + 1
      ipar = NM%nc + ip
      PM(ipar)%pname = 'HTG'//htg(i:i)//trim(LCF%htgmod)
      PM(ipar)%ptype = 'P'
      PM(ipar)%ipt = ipar
      PM(ipar)%psat = 0
      if (LCF%htgmod(1:3) .eq. 'STO') then
        PM(ipar)%ptime(2) = LCF%jd0 + LCF%sod0/86400.d0
      else if (LCF%htgmod(1:3) .eq. 'PWC' .and. kpar .gt. 1) then
        PM(ipar)%ptime(2) = LCF%jd0 + LCF%sod0/86400.d0 + minut/1440.d0
      endif
      PM(ipar)%map = 1.d0
      if (LCF%htgmod(1:3) .eq. 'STO') then
        PM(ipar)%rw = 1.d0/(SITE%qhtg*dsqrt(LCF%dintv/3600.d0))
      else if (LCF%htgmod(1:3) .eq. 'PWC') then
        PM(ipar)%rw = 1.d0/(SITE%qhtg*dsqrt(minut/720.d0))        !m/sqrt(12 hour)
      endif
      PM(ipar)%xini = 0.d0
      OB%npar = OB%npar + 1
      OB%pname(OB%npar) = PM(ipar)%pname
      OB%ltog(OB%npar, 1) = ipar
      NM%norx(ipar, ipar) = 1.d0/SITE%dhtg0**2
      NM%iptp(ipar) = ipar
    enddo
  endif
  NM%ipm = NM%imtx
!
!! copy common parameter index
  do ipar = 1, OB%npar
    do isat = 2, LCF%nprn
      OB%ltog(ipar, isat) = OB%ltog(ipar, 1)
    enddo
  enddo
!
!! ambiguity (dynamic)
  OB%npar = OB%npar + 1
  OB%pname(OB%npar) = 'AMBC'
  do isat = 1, LCF%nprn
    OB%ltog(OB%npar, isat) = 0
    OB%lifamb(isat, 1) = 0.d0
    OB%lifamb(isat, 2) = 0.d0
  enddo
!
!! check consistence
  if (ic .ne. NM%nc .or. ip .ne. NM%np .or. ic + ip .ne. NM%imtx) then
    write (*, '(a)') '***ERROR(lsq_init): parameter number not consistent with that from lsq_cnt_prmt'
    write (*, '(2(a,2i5,/))') ' nc, ic ', NM%nc, ic, ' np, ip ', NM%np, ip
    call exit(1)
  endif
!
!! output to screen
  write (*, '(a)') ' ALL PARAMETERS '
  do ipar = 1, NM%imtx
    val = NM%norx(ipar, ipar)
    sitename = SITE%name
    if (PM(ipar)%psat .eq. 0) then
      prnname = '   '
    else
      write (prnname, '(i3)') PM(ipar)%psat
    endif
    write (*, '(2i4,1x,a4,1x,a3,1x,a15)') ipar, PM(ipar)%psat, sitename, prnname, PM(ipar)%pname
  enddo

  return
end
