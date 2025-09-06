

!! lsq_init.f90
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
!! Contributor: Maorong Ge, Jianghui Geng, Songfeng Yang, Jihang Lin, Jing Zeng
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
subroutine lsq_init(LCF, SITE, OB, NM, PM)
  implicit none
  include '../header/const.h'
  include '../header/station.h'
  include '../header/rnxobs.h'
  include 'lsqcfg.h'
  include 'lsq.h'

  type(lsqcfg) LCF
  type(station) SITE
  type(rnxobr) OB
  type(norm) NM
  type(prmt) PM(1:*)
!
!! local
  integer*4 isat, ipar, irow
  integer*4 i, j, k, ic, ip, kpar, sec, refrck 
  character*4 sitename
  character*3 xyz, prnname
  character*2 htg, tmpsys
  character*6 sys
  real*8 val

  data xyz, htg, sys/'XYZ', 'CS', 'GREC3J'/
!
!! initialization
  do ipar = 1, NM%imtx + 1
    if (ipar .ne. NM%imtx + 1) then
      PM(ipar)%ptype = ' '
      PM(ipar)%ipt = 0
      PM(ipar)%iepo = 0
      PM(ipar)%iobs = 0
      PM(ipar)%iobs_G = 0
      PM(ipar)%iobs_R = 0
      PM(ipar)%iobs_E = 0
      PM(ipar)%iobs_C = 0
      PM(ipar)%iobs_3 = 0
      PM(ipar)%iobs_J = 0
      PM(ipar)%ptime(1) = LCF%jd0 + LCF%sod0/864.d2
      PM(ipar)%ptime(2) = LCF%jd1 + LCF%sod1/864.d2
      PM(ipar)%ptbeg = 0.d0
      PM(ipar)%ptend = 0.d0
      PM(ipar)%map = 0.d0
      PM(ipar)%rw = 0.d0
      PM(ipar)%xcor = 0.d0
      PM(ipar)%xsig = 0.d0
      PM(ipar)%xrms = 0.d0
      PM(ipar)%xrwl = 0.d0      ! especially Melbourne-Wubbena
      PM(ipar)%xswl = 0.d0
      PM(ipar)%mele = 0.d0
      PM(ipar)%zw = 0.d0        ! initial value of wide-lane
    end if
    do irow = 1, NM%imtx
      NM%norx(irow, ipar) = 0.d0
    end do
  end do

  do ipar = 1, MAXPAR_STA
    do isat = 1, MAXSAT
      OB%ltog(ipar, isat) = 0
    end do
  end do
  do isat = 1, MAXSAT
    OB%delay(isat) = 0.d0
    OB%flag(isat) = 0
  end do
  ic = 0
  ip = 0
  sec = 0
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
      PM(ipar)%xini = SITE%x(i)*1.d3
      OB%npar = OB%npar + 1
      OB%pname(OB%npar) = PM(ipar)%pname
      OB%ltog(OB%npar, 1) = ipar
      NM%norx(ipar, ipar) = 1.d0/SITE%dx0(i)**2
      NM%iptp(ipar) = ipar
    end do
  else if (SITE%skd(1:1) .eq. 'P') then
    do i = 1, 3
      ip = ip + 1
      ipar = NM%nc + ip
      write (PM(ipar)%pname, '(a,i0)') 'STAP'//xyz(i:i)//':', SITE%pospd
      PM(ipar)%ptype = 'P'
      PM(ipar)%ipt = ipar
      PM(ipar)%psat = 0
      PM(ipar)%xini = SITE%x(i)*1.d3
      PM(ipar)%ptime(1) = LCF%jd0 + LCF%sod0/864.d2 
      PM(ipar)%ptime(2) = LCF%jd0 + LCF%sod0/864.d2 + SITE%pospd/864.d2 
      PM(ipar)%map = 1.d0
      PM(ipar)%rw = 1.d0/(SITE%rx0*dsqrt(SITE%pospd*1.d0))
      OB%npar = OB%npar + 1
      OB%pname(OB%npar) = PM(ipar)%pname
      OB%ltog(OB%npar, 1) = ipar
      NM%norx(ipar, ipar) = 1.d0/SITE%dx0(i)**2
      NM%iptp(ipar) = ipar
    end do
  else if (SITE%skd(1:1) .eq. 'K'.or. SITE%skd(1:1) .eq. 'L') then
    do i = 1, 3
      ip = ip + 1
      ipar = NM%nc + ip
      PM(ipar)%pname = 'STAP'//xyz(i:i)
      PM(ipar)%ptype = 'P'
      PM(ipar)%ipt = ipar
      PM(ipar)%psat = 0
      PM(ipar)%xini = 1.d3
      PM(ipar)%ptime(1) = LCF%jd0 + LCF%sod0/864.d2 
      PM(ipar)%ptime(2) = LCF%jd0 + LCF%sod0/864.d2 
      PM(ipar)%map = 0.d0
      PM(ipar)%rw = 1.d0/SITE%dx0(i)
      OB%npar = OB%npar + 1
      OB%pname(OB%npar) = PM(ipar)%pname
      OB%ltog(OB%npar, 1) = ipar
      NM%norx(ipar, ipar) = 1.d0/SITE%dx0(i)**2
      NM%iptp(ipar) = ipar
    end do
  else
    write (*, '(2a,1x,a)') '***ERROR(lsq_init): unknown positioning mode ', SITE%name, SITE%skd
    call exit(1)
  end if
!
!! receiver clock offset
  refrck = 0; tmpsys = ''
  do i = 1, 6
    if (index(LCF%sys, sys(i:i)) .ne. 0) then
      if (index(LCF%isbsys, sys(i:i)) .eq. 0)then
        if(refrck.eq.0) refrck = i 
        ip = ip + 1
        ipar = NM%nc + ip
        PM(ipar)%pname = 'RECCLK_'//sys(i:i)//'_'//trim(LCF%rckmod)
        PM(ipar)%ptype = 'P'
        PM(ipar)%ipt = ipar
        PM(ipar)%psat = 0
        if (LCF%rckmod(1:3) .eq. 'STO') then
          PM(ipar)%map = 1.d0
          PM(ipar)%ptime(2) = LCF%jd0 + LCF%sod0/864.d2
          PM(ipar)%rw = 1.d0/(SITE%qrck*dsqrt(LCF%dintv/3600.d0))
        else if (LCF%rckmod(1:3) .eq. 'WNO') then
          PM(ipar)%map = 0.d0
          PM(ipar)%ptime(2) = LCF%jd0 + LCF%sod0/864.d2
          PM(ipar)%rw = 1.d0/SITE%dclk0
        end if
        PM(ipar)%xini = 0.d0
        OB%npar = OB%npar + 1
        OB%pname(OB%npar) = PM(ipar)%pname
        do isat = 1, LCF%nprn
          tmpsys(1:1) = LCF%prn(isat) (1:1)
          if (LCF%prn(isat) (1:1) .eq. 'C') then
            read (LCF%prn(isat) (2:3), '(i2)') k
            if (k .gt. 17) tmpsys(1:1) = '3'
          endif
          if (tmpsys(1:1) .eq. sys(i:i)) then
            OB%ltog(OB%npar, isat) = ipar
          else if (index(LCF%isbsys, tmpsys(1:1)) .ne. 0 .and. i.eq.refrck)then
            OB%ltog(OB%npar, isat) = ipar 
          endif
        end do
        NM%norx(ipar, ipar) = 1.d0/SITE%dclk0**2
        NM%iptp(ipar) = ipar
      else
        ip = ip + 1
        ipar = NM%nc + ip
        write (PM(ipar)%pname, '(a,i0)') 'RECISB_'//sys(i:i)//':', 86400
        PM(ipar)%ptype = 'P'
        PM(ipar)%ipt = ipar
        PM(ipar)%psat = 0
        PM(ipar)%map = 0.d0
        PM(ipar)%ptime(2) = LCF%jd0 + LCF%sod0/864.d2 + 864.d2/864.d2
        PM(ipar)%rw = 1.d0/SITE%dclk0
        PM(ipar)%xini = 0.d0
        OB%npar = OB%npar + 1
        OB%pname(OB%npar) = PM(ipar)%pname
        do isat = 1, LCF%nprn
          if (LCF%prn(isat) (1:1) .eq. 'C') then
            if (sys(i:i) .eq. 'C') then
              read (LCF%prn(isat) (2:3), '(i2)') k
              if (k .le. 17) OB%ltog(OB%npar, isat) = ipar
            else if (sys(i:i) .eq. '3') then
              read (LCF%prn(isat)(2:3), '(i2)') k
              if (k .gt. 17) OB%ltog(OB%npar, isat) = ipar
            end if
            cycle
          end if
          if (LCF%prn(isat)(1:1) .eq. sys(i:i)) OB%ltog(OB%npar, isat) = ipar
        end do
        NM%norx(ipar, ipar) = 1.d0/SITE%dclk0**2
        NM%iptp(ipar) = ipar
          
      endif
    end if
  end do

!
!! disable atmosphere models for LEO satellites
  if (SITE%skd(1:1) .eq. 'L') goto 50

!
!! zenith atmospheric delay
  if (LCF%ztdmod(1:3) .ne. 'NON') then
    if (LCF%ztdmod(1:3) .eq. 'PWC') then
      k = index(LCF%ztdmod, ':')
      read (LCF%ztdmod(k + 1:), *) sec 
      sec = sec * 60
      kpar = nint(((LCF%jd1 - LCF%jd0)*864.d2 + (LCF%sod1 - LCF%sod0))/sec)
    end if
    ip = ip + 1
    ipar = NM%nc + ip
    if (LCF%ztdmod(1:3) .eq. 'STO') then
      PM(ipar)%pname = 'ZTDSTO'
    else if (LCF%ztdmod(1:3) .eq. 'PWC') then
      write (PM(ipar)%pname, '(a,i0)') 'ZTDPWC:', sec
    end if
    PM(ipar)%ptype = 'P'
    PM(ipar)%ipt = ipar
    PM(ipar)%psat = 0
    if (LCF%ztdmod(1:3) .eq. 'STO') then
      PM(ipar)%ptime(2) = LCF%jd0 + LCF%sod0/864.d2
    else if (LCF%ztdmod(1:3) .eq. 'PWC' .and. kpar .gt. 1) then
      PM(ipar)%ptime(2) = LCF%jd0 + LCF%sod0/864.d2 + sec/864.d2
    end if
    PM(ipar)%map = 1.d0
    if (LCF%ztdmod(1:3) .eq. 'STO') then
      PM(ipar)%rw = 1.d0/(SITE%qztd*dsqrt(LCF%dintv/3600.d0))
    else if (LCF%ztdmod(1:3) .eq. 'PWC') then
      PM(ipar)%rw = 1.d0/(SITE%qztd*dsqrt(sec/3600.d0))
    end if
    PM(ipar)%xini = 0.d0
    OB%npar = OB%npar + 1
    OB%pname(OB%npar) = PM(ipar)%pname
    OB%ltog(OB%npar, 1) = ipar
    NM%norx(ipar, ipar) = 1.d0/SITE%dztd0**2
    NM%iptp(ipar) = ipar
  end if

!
!! horizontal troposphere gradients
  if (LCF%htgmod(1:3) .ne. 'NON') then
    if (LCF%htgmod(1:3) .eq. 'PWC') then
      k = index(LCF%htgmod, ':')
      read (LCF%htgmod(k + 1:), *) sec
      sec = sec * 60 
      kpar = nint(((LCF%jd1 - LCF%jd0)*864.d2 + (LCF%sod1 - LCF%sod0))/sec)
    end if
    do i = 1, 2
      ip = ip + 1
      ipar = NM%nc + ip
      if (LCF%htgmod(1:3) .eq. 'STO') then
        PM(ipar)%pname = 'HTG'//htg(i:i)//'STO'
      else if (LCF%htgmod(1:3) .eq. 'PWC') then
        write (PM(ipar)%pname, '(a,i0)') 'HTG'//htg(i:i)//'PWC:', sec
      end if
      PM(ipar)%ptype = 'P'
      PM(ipar)%ipt = ipar
      PM(ipar)%psat = 0
      if (LCF%htgmod(1:3) .eq. 'STO') then
        PM(ipar)%ptime(2) = LCF%jd0 + LCF%sod0/864.d2
      else if (LCF%htgmod(1:3) .eq. 'PWC' .and. kpar .gt. 1) then
        PM(ipar)%ptime(2) = LCF%jd0 + LCF%sod0/864.d2 + sec/864.d2
      end if
      PM(ipar)%map = 1.d0
      if (LCF%htgmod(1:3) .eq. 'STO') then
        PM(ipar)%rw = 1.d0/(SITE%qhtg*dsqrt(LCF%dintv/3600.d0))
      else if (LCF%htgmod(1:3) .eq. 'PWC') then
        PM(ipar)%rw = 1.d0/(SITE%qhtg*dsqrt(sec/43200.d0))
      end if
      PM(ipar)%xini = 0.d0
      OB%npar = OB%npar + 1
      OB%pname(OB%npar) = PM(ipar)%pname
      OB%ltog(OB%npar, 1) = ipar
      NM%norx(ipar, ipar) = 1.d0/SITE%dhtg0**2
      NM%iptp(ipar) = ipar
    end do
  end if

50 continue
  NM%ipm = NM%imtx

!
!! copy common parameter index
  do ipar = 1, OB%npar
    if (OB%pname(ipar) (1:6) .ne. 'RECCLK' .and. OB%pname(ipar) (1:6) .ne. 'RECISB') then
      do isat = 2, LCF%nprn
        OB%ltog(ipar, isat) = OB%ltog(ipar, 1)
      end do
    end if
  end do
!
!! ambiguity (dynamic)
  OB%npar = OB%npar + 1
  OB%pname(OB%npar) = 'AMBC'
  do isat = 1, LCF%nprn
    OB%ltog(OB%npar, isat) = 0
    OB%lifamb(isat, 1) = 0.d0
    OB%lifamb(isat, 2) = 0.d0
  end do
!
!! check consistence
  if (ic .ne. NM%nc .or. ip .ne. NM%np .or. ic + ip .ne. NM%imtx) then
    write (*, '(a)') '***ERROR(lsq_init): parameter number not consistent with that from lsq_cnt_prmt'
    write (*, '(2(a,2i5,/))') ' nc, ic ', NM%nc, ic, ' np, ip ', NM%np, ip
    call exit(1)
  end if
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
    end if
    write (*, '(2i4,1x,a4,1x,a3,1x,a15)') ipar, PM(ipar)%psat, sitename, prnname, PM(ipar)%pname
  end do

  return
end subroutine
