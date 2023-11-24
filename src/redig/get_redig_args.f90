!
!! get_redig_args.f90
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
!! motified by: Songfeng Yang : add multisystem-GREC
!!                           Email:sfyang@whu.edu.cn
!!
!!
!! purpose  : get arguments of redig
!! parameter:
!!            input : nepo -- # of epochs in residual file
!!            output: RCF  -- redig configure struct
!
subroutine get_redig_args(nepo, RCF)
  implicit none
  include '../header/const.h'
  include 'rescfg.h'

!
!! common
  integer*4 idxfrq(MAXSYS, 2)
  common idxfrq
!
!! parameter
  integer*4 nepo
  type(rescfg) RCF
!
!! local
  integer*4 i, j, nargs, isat
  integer*4 iy, imon, id, ih, im, is, iy_ses, imon_ses, id_ses, ih_ses, im_ses, is_ses, ierr
  integer*4 jy, jmon, jday, jdoy
  character*3 prn(MAXSAT), temprn, prn_mat(MAXSAT)
  real*8 sec, sec_ses, seslen, seslen_ses
  character*20 resfil, sttfil
  character*256 line
!
!! function called
  integer*4 get_valid_unit, modified_julday, pointer_string
  character*4 lower_string

  RCF%lfnres = 0
  RCF%lfnstt = 0
  RCF%jump = 0.d0
  RCF%pjump = 0.d0
  RCF%nsht = 0
  RCF%pcedit = .false.
  RCF%amb_at_dbd = 'DEFAULT'
  call prn_matbld(prn_mat)
!
!! read arguments
  nargs = iargc()
  if (nargs .eq. 0) then
    write (*, '(a)') 'Usage: redig resfil [-jmp jump -sht nsht -pce pcedit]'
    call exit(4)
  end if
  call getarg(1, resfil)
  i = 2
  do while (i .le. nargs)
    call getarg(i, line)
    i = i + 1
    if (line(1:4) .eq. '-jmp') then
      call getarg(i, line)
      read (line, *, err=100) RCF%jump
      RCF%pjump = RCF%jump/10.d0
    else if (line(1:4) .eq. '-sht') then
      call getarg(i, line)
      read (line, *, err=100) RCF%nsht
    else if (line(1:4) .eq. '-pce') then
      RCF%pcedit = .true.
    end if
    i = i + 1
  end do
!
!! read header of residuals
  RCF%lfnres = get_valid_unit(10)
  open (RCF%lfnres, file=resfil, status='old', iostat=ierr)
  if (ierr .ne. 0) then
    write (*, '(2a)') '***ERROR(get_redig_args): open file ', trim(resfil)
    call exit(1)
  end if
  line = ' '
  RCF%snam = ' '
  isat = 1
  do i = 1, MAXSAT
    RCF%prn(i) = ''
  end do
  prn = ''
  do while (index(line, 'END OF HEADER') .eq. 0)
    read (RCF%lfnres, '(a)') line
    if (index(line, '# OF SAT') .ne. 0) then
      read (line, '(i10)') RCF%nprn
    else if (index(line, 'STATION') .ne. 0) then
      read (line, '(a4)') RCF%snam
    else if (index(line, 'SATELLITE LIST') .ne. 0) then
      read (line, '(15(a3,1x))') prn(isat:isat + 14)
      do i = isat, isat + 15
        if (prn(i) .eq. '') then
          exit
        else
          RCF%prn(i) = prn(i)
        end if
      end do
      isat = i
    else if (index(line, 'INT / OBS TYPE') .ne. 0) then
      read (line, '(f10.2,10x,a11)') RCF%dintv, RCF%obstyp
    else if (index(line, 'RES TIME BEG/LEN') .ne. 0) then
      read (line, *) iy, imon, id, ih, im, sec, seslen
      RCF%jd0 = modified_julday(id, imon, iy)
      RCF%sod0 = ih*3600.d0 + im*60.d0 + sec
      nepo = nint(seslen/RCF%dintv) + 1
    else if (index(line, 'SYS / FREQUENCY BAND') .ne. 0) then
      do i = 1, MAXSYS
        if (line(1:3) .eq. GNSS_NAME_LEN3(i)) then
          do j = 1, MAXFRQ
            if (line(11:13) .eq. FREQ_NAME_SYS(j, i)) idxfrq(i, 1) = j
            if (line(16:18) .eq. FREQ_NAME_SYS(j, i)) idxfrq(i, 2) = j
          end do
          exit
        end if
      end do
    else if (index(line, 'AMB AT DAT-BOUNDARY') .ne. 0) then
      read (line, *) RCF%amb_at_dbd
    end if
  end do
!
!! read year/month/day from file name instead of records
  jy = 0
  jmon = 0
  jday = 0

  i = index(resfil, 'res_')
  read (resfil(i:), '(4xi4i3)', err=200) jy, jdoy
  call yeardoy2monthday(jy, jdoy, jmon, jday)

200 continue
  if (jy .le. 0 .or. jmon .le. 0 .or. jmon .gt. 12 .or. jday .le. 0 .or. jday .gt. 31) then
    jy = iy
    jmon = imon
    jday = id
  end if
!
!! status file
  call file_name(.false., 'stt', ' ', jy, jmon, jday, ih, sttfil)
  RCF%lfnstt = get_valid_unit(10)
  open (RCF%lfnstt, file=sttfil)
!
!! rinex health diagnose file
  call file_name(.false., 'log', 'SNAM='//lower_string(RCF%snam), jy, jmon, jday, ih, RCF%flnrhd)
  RCF%lfnrhd = get_valid_unit(10)
  open (RCF%lfnrhd, file=RCF%flnrhd, status='old', action='read', iostat=ierr)
  if (ierr .ne. 0) then
    write (*, '(2a)') '***ERROR(get_redig_args): open file ', trim(RCF%flnrhd)
    call exit(1)
  end if
!
!! update prn according to log
  call log_sat(RCF, nepo)
  ! sort prn
  do i = 1, RCF%nprn - 1
    do j = i + 1, RCF%nprn
      id = pointer_string(MAXSAT, prn_mat, RCF%prn(i))
      ih = pointer_string(MAXSAT, prn_mat, RCF%prn(j))
      if (id .gt. ih) then
        temprn = RCF%prn(i)
        RCF%prn(i) = RCF%prn(j)
        RCF%prn(j) = temprn
      end if
    end do
  end do

  return
100 write (*, '(2a)') '***ERROR(get_redig_args): read arguments ', trim(line)
  call exit(1)
end subroutine
