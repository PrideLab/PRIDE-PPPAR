!
!! write_diag_rpt.f90
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
!! purpose  : write rinex health diagnose report
!! parameter:
!!    input : flnrhd   -- diagnose file
!!            nepo     -- # of epochs
!!            nsat     -- # of satellites
!!            jd0,ti   -- time tag
!!            flagall  -- flags
!!            interval -- sampling interval
!
subroutine write_diag_rpt(flnrhd, nepo, nsat, stanam, jd0, ti, session_length, flagall, interval)
  implicit none
  include '../header/const.h'

  integer*4 nepo, nsat, jd0, flagall(nepo, MAXSAT)
  real*8 ti(1:*), interval, session_length
  character*(*) flnrhd,stanam
!
!! local
  integer*4 i, j, ns, lfn, iepo, isat, ie, jd, iy, imon, id, ih, im
  integer*4 maxamb, totamb, remobs, avaobs, xamb(nepo), flag(MAXSAT, nepo)
  real*8 sec, sod
  character*256 string, str(MAXSAT)
  character*3 prn0(MAXSAT)
!
!! function called
  logical*1 istrue
  integer*4 get_valid_unit, find_ambd
  
  call prn_matbld(prn0)
!
!! initialization
  maxamb = 0
  totamb = 0
  remobs = 0
  avaobs = 0
  flag = transpose(flagall)
  do i = 1, nepo
    xamb(i) = 0
  enddo
!
!! statistics about rinex health
  do isat = 1, nsat
    do iepo = 1, nepo
      if (istrue(flagall(iepo, isat), 'nodata')) cycle
      if (istrue(flagall(iepo, isat), 'ok')) then
        avaobs = avaobs + 1
        if (istrue(flagall(iepo, isat), 'amb')) then
          totamb = totamb + 1
          j = find_ambd(nepo, flagall(1, isat), iepo)
          do i = iepo, j
            xamb(i) = xamb(i) + 1
          enddo
        endif
      else
        remobs = remobs + 1
      endif
    enddo
  enddo
  maxamb = maxval(xamb, 1)
!
!! write header
  lfn = get_valid_unit(10)
  open (lfn, file=flnrhd)
  write (lfn, '(a20,6x,a4,30x,a)') 'RINEX health Logfile', stanam, 'COMMENT'
  write (lfn, '(f10.2,50x,a)') interval, 'INTERVAL'
  write (lfn, '(3i10,30x,a)') maxamb, totamb, 0, 'AMB MAX/TOT/NEW'
  write (lfn, '(3i10,30x,a)') avaobs, remobs, 0, 'EPO AVA/REM/NEW'
  call mjd2date(jd0, ti(1), iy, imon, id, ih, im, sec)
  write (lfn, '(i5,4i3,f11.7,f12.2,20x,a)') iy, imon, id, ih, im, sec, session_length, 'RES TIME BEG/LEN'
  write (lfn, '(60x,a)') 'END OF HEADER'
!
!! epoch diagnose
  do iepo = 1, nepo
    ns = 0
    do isat = 1, nsat
      if (istrue(flag(isat, iepo), 'nodata') .or. istrue(flag(isat, iepo), 'good')) cycle
      if (.not. istrue(flag(isat, iepo), 'ok')) then
        ns = ns + 1
        string = ' '
        if (istrue(flag(isat, iepo), 'lowele')) then
          string = 'LOWELEVATION'
        else if (istrue(flag(isat, iepo), 'shrt')) then
          string = 'SHORTPIECE'
        else if (istrue(flag(isat, iepo), 'lwbad')) then
          string = 'BADWIDELANE'
        else if (istrue(flag(isat, iepo), 'lgbad')) then
          string = 'BADIONOSPHERE'
        else if (istrue(flag(isat, iepo), 'lccheck')) then
          string = 'CANNOTCHECKLC'
        else if (istrue(flag(isat, iepo), 'no4')) then
          string = 'LESSTHAN4OBS'
        else if (istrue(flag(isat, iepo), 'pcbad')) then
          string = 'BADRANGE'
        else if (istrue(flag(isat, iepo), 'pc1ms')) then
          string = 'BAD1MS'
        endif
        write (str(ns), '(a3,57x,a)') prn0(isat), 'DEL_'//trim(string) ! add multisystem-GREC
      else if (istrue(flag(isat, iepo), 'amb')) then
        ns = ns + 1
        string = ' '
        if (istrue(flag(isat, iepo), 'lli')) then
          string = 'FLAGINRINEX'
        else if (istrue(flag(isat, iepo), 'bigsd')) then
          string = 'BIGLCJUMP'
        else if (istrue(flag(isat, iepo), 'gap')) then
          string = 'BIGGAP'
        else if (istrue(flag(isat, iepo), 'lwjump')) then
          string = 'WIDELANEJUMP'
        else if (istrue(flag(isat, iepo), 'lwconn')) then
          string = 'IONO.FAILED'
        else if (istrue(flag(isat, iepo), 'lgjump')) then
          string = 'IONO.JUMP'
        endif
        ie = find_ambd(nepo, flagall(1, isat), iepo)
        jd = jd0 + int(ti(ie)/86400.d0)
        sod = ti(ie) - (jd - jd0)*86400.d0
        call mjd2date(jd, sod, iy, imon, id, ih, im, sec)
        write (str(ns), '(a3,28x,i5,4i3,f11.7,1x,a)') prn0(isat), iy, imon, id, ih, im, sec, &
          'AMB_'//trim(string)
      endif
    enddo
    if (ns .eq. 0) cycle
    jd = jd0 + int(ti(iepo)/86400.d0)
    sod = ti(iepo) - (jd - jd0)*86400.d0
    call mjd2date(jd, sod, iy, imon, id, ih, im, sec)
    write (lfn, '(a3,i5,4i3,f11.7)') 'TIM', iy, imon, id, ih, im, sec
    do i = 1, ns
      write (lfn, '(a)') str(i) (1:len_trim(str(i)))
    enddo
  enddo

  close (lfn)
  return
end
