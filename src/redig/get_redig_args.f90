!
!! get_redig_args.f90
!!
!!    Copyright (C) 2022 by Wuhan University
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
!! Contributor: Maorong Ge, Jianghui Geng, Songfeng Yang, Jihang Lin
!! 
!!
!!
!! motified by: Songfeng Yang : add multisystem-GREC
!!                           Email:sfyang@whu.edu.cn
!!
!!
!! purpose  : get arguments of redig
!! parameter:
!!    input : nepo -- # of epochs in residual file
!!    output: RCF  -- redig configure struct
!
subroutine get_redig_args(nepo, RCF)
  implicit none
  include '../header/const.h'
  include 'rescfg.h'

  integer*4 nepo
  type(rescfg) RCF
!
!! local
  integer*4 i, j, nargs, isat
  integer*4 iy, imon, id, ih, im, is, iy_ses, imon_ses, id_ses, ih_ses, im_ses, is_ses, ierr
  integer*4 jy, jmon, jday, jdoy
  character*3 prn(MAXSAT),temprn,prn_mat(MAXSAT)
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
  RCF%nsht = 0
  call prn_matbld(prn_mat)
!
!! read arguments
  nargs = iargc()
  if (nargs .eq. 0) then
    write (*, '(a)') 'Usage: redig resfil [-jmp jump -sht nsht]'
    call exit(4)
  endif
  call getarg(1, resfil)
  i = 2
  do while (i .le. nargs)
    call getarg(i, line)
    i = i + 1
    if (line(1:4) .eq. '-jmp') then
      call getarg(i, line)
      read (line, *, err=100) RCF%jump
    else if (line(1:4) .eq. '-sht') then
      call getarg(i, line)
      read (line, *, err=100) RCF%nsht
    endif
    i = i + 1
  enddo
!
!! read header of residuals
  RCF%lfnres = get_valid_unit(10)
  open (RCF%lfnres, file=resfil, status='old', iostat=ierr)
  if (ierr .ne. 0) then
    write (*, '(2a)') '***ERROR(get_redig_args): open file ', trim(resfil)
    call exit(1)
  endif
  line = ' '
  RCF%snam = ' '
  isat = 1
  do i = 1, MAXSAT
    RCF%prn(i) = ''
  enddo
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
        endif
      enddo
      isat = i
    else if (index(line, 'INT / OBS TYPE') .ne. 0) then
      read (line, '(f10.2,10x,a11)') RCF%dintv, RCF%obstyp
    else if (index(line, 'RES TIME BEG/LEN') .ne. 0) then
      read (line, *) iy, imon, id, ih, im, sec, seslen
      RCF%jd0 = modified_julday(id, imon, iy)
      RCF%sod0 = ih*3600.d0 + im*60.d0 + sec
      nepo = nint(seslen/RCF%dintv) + 1
    endif
  enddo
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
  endif
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
  endif
!
!! update prn according to log
  call log_sat(RCF,nepo)
  ! sort prn
  do i=1,RCF%nprn-1
    do j=i+1,RCF%nprn
      id=pointer_string(MAXSAT, prn_mat, RCF%prn(i))
      ih=pointer_string(MAXSAT, prn_mat, RCF%prn(j))
      if(id .gt. ih) then
        temprn=RCF%prn(i)
        RCF%prn(i)=RCF%prn(j)
        RCF%prn(j)=temprn
      endif
    enddo
  enddo

  return
100 write (*, '(2a)') '***ERROR(get_redig_args): read arguments ', trim(line)
  call exit(1)
end
