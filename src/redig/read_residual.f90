!
!! read_residual.f90
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
!! purpose  : read residuals
!! parameter:
!!    input : nepo -- number of epochs
!!            RCF  -- anares control
!!    output: resi -- residuals
!!            flag -- flag info
!!            trsi -- time tag for each epoch
!
subroutine read_residual(nepo, resi, flag, trsi, RCF)
  implicit none
  include '../header/const.h'
  include 'rescfg.h'
  include 'data_flag.h'

  integer*4 nepo
  type(rescfg) RCF
  integer*4 flag(nepo, RCF%nprn)
  real*8 resi(nepo, RCF%nprn)
  character*27 trsi(nepo)
!
!! local
  integer*4 i, j, jdx, iy, imon, id, ih, im, iepo, isat, ierr
  real*8 sec, sodx, phs, rag, elev, azim, rmstbl(0:MAXSAT)
  character*3 iprn
  character*1024 line, line1
!
!! function called
  logical*1 istrue
  integer*4 modified_julday
  real*8 timdif
  integer*4 pointer_string
!
!! initialization
  do isat = 1, RCF%nprn
    do iepo = 1, nepo
      resi(iepo, isat) = 0.d0
      flag(iepo, isat) = NODATA
    enddo
    rmstbl(isat) = 0.d0
  enddo
  do iepo = 1, nepo
    trsi(iepo) = ' '
  enddo
!
!! read residuals
  do while (.true.)
    read (RCF%lfnres, '(a)', iostat=ierr) line
    if (ierr .ne. 0) exit
    if (line(1:3) .ne. 'TIM') cycle
    read (line, '(3x,i5,4i3,f11.7)', err=100) iy, imon, id, ih, im, sec
    jdx = modified_julday(id, imon, iy)
    sodx = ih*3600.d0 + im*60.d0 + sec
    iepo = nint(timdif(jdx, sodx, RCF%jd0, RCF%sod0)/RCF%dintv) + 1
    trsi(iepo) = line(5:31)
    if (iepo .gt. nepo) then
      write (*, '(a)') '***ERROR(read_residual): number of epochs exceeded '
      call exit(1)
    endif
    do while (.true.)
      read (RCF%lfnres, '(a)', iostat=ierr) line
      if (ierr .ne. 0) exit
      if (line(1:3) .eq. 'TIM') then
        backspace RCF%lfnres
        exit
      elseif (line(1:3).eq.'CST') then
        cycle
      endif
      read (line, '(a3)', iostat=ierr) iprn
      isat=pointer_string(RCF%nprn,RCF%prn,iprn)
      read (line, '(3x,f10.3,42x,i3,f8.3,f9.3)', iostat=ierr) phs, flag(iepo, isat), elev, azim
      if (ierr .ne. 0) then
        flag(iepo, isat) = 4
        phs = 0.d0
      endif
      if(index(RCF%obstyp,'LC').ne.0) resi(iepo,isat)=phs*1.d3
    enddo
  enddo
  close (RCF%lfnres)
  RCF%lfnres = 0

  if (RCF%lfnstt .eq. 0) return
!
!! reports about residuals
  rmstbl(0) = 0.d0
  j = 0
  do isat = 1, RCF%nprn
    i = 0
    do iepo = 1, nepo
      if (istrue(flag(iepo, isat), 'OK')) then
        i = i + 1
        rmstbl(isat) = rmstbl(isat) + resi(iepo, isat)**2
      endif
    enddo
    j = j + i
    rmstbl(0) = rmstbl(0) + rmstbl(isat)
    if (i .gt. 0) rmstbl(isat) = dsqrt(rmstbl(isat)/i)
  enddo
  if (j .gt. 0) rmstbl(0) = dsqrt(rmstbl(0)/j)
  if (index(RCF%obstyp, 'LC') .ne. 0) then
    write (RCF%lfnstt, '(a)') '+RMS OF RESIDUALS---PHASE(MM)'
  else
    write (RCF%lfnstt, '(a)') '+RMS OF RESIDUALS---RANGE(DM)'
  endif
  write (RCF%lfnstt, '(a4,1x,a4,190(1x,a3))') 'NAME', 'SUMM', (RCF%prn(i), i=1, RCF%nprn)
  write (RCF%lfnstt, '(a4,1x,i4,190(1x,i3))') RCF%snam, (nint(rmstbl(i)), i=0, RCF%nprn)
  write (RCF%lfnstt, '(a4,1x,a4,190(1x,a3))') 'NAME', 'SUMM', (RCF%prn(i), i=1, RCF%nprn)
  if (index(RCF%obstyp, 'LC') .ne. 0) then
    write (RCF%lfnstt, '(a)') '-RMS OF RESIDUALS---PHASE(MM)'
  else
    write (RCF%lfnstt, '(a)') '-RMS OF RESIDUALS---RANGE(DM)'
  endif
!
!! residual time series
  if (index(RCF%obstyp, 'LC') .ne. 0) then
    write (RCF%lfnstt, '(a)') '+TIME SERIES OF RESIDUALS---PHASE(MM)'
  else
    write (RCF%lfnstt, '(a)') '+TIME SERIES OF RESIDUALS---RANGE(DM)'
  endif
  write (RCF%lfnstt, '(2x,a4,190(1x,a3))') RCF%snam, (RCF%prn(i), i=1, RCF%nprn)
  do iepo = 1, nepo
    line = ' '
    write (line, '(i6)') iepo
    j = 6
    do isat = 1, RCF%nprn
      if (istrue(flag(iepo, isat), 'OK')) then
        write (line(j + 1:), '(a1,i3)') ' ', nint(resi(iepo, isat))
      else
        write (line(j + 1:), '(a4)') '    '
      endif
      j = j + 4
    enddo
    write (RCF%lfnstt, '(a)') trim(line)
  enddo
  write (RCF%lfnstt, '(2x,a4,190(1x,a3))') RCF%snam, (RCF%prn(i), i=1, RCF%nprn)
  if (index(RCF%obstyp, 'LC') .ne. 0) then
    write (RCF%lfnstt, '(a)') '-TIME SERIES OF RESIDUALS---PHASE(MM)'
  else
    write (RCF%lfnstt, '(a)') '-TIME SERIES OF RESIDUALS---RANGE(DM)'
  endif
  close (RCF%lfnstt)
  RCF%lfnstt = 0
  write (*, '(a)') '%%%+RMS OF RESIDUALS---PHASE(MM)'
  isat = 1
  line = ' '
  line1 = ' '
  j = 0
  do isat = 1, RCF%nprn
    if (rmstbl(isat) .gt. 0 .and. RCF%prn(isat)(1:1) .eq. 'G') then
      write (line(j + 1:), '(a3,a1)') RCF%prn(isat), ' '
      write (line1(j + 1:), '(i3,a1)') nint(rmstbl(isat)), ' '
      j = j + 4
    endif
  enddo
  if (line .ne. '') then
    write (*, '(a)') trim(line)
    write (*, '(a)') trim(line1)
  endif
  isat = 1
  line = ' '
  line1 = ' '
  j = 0
  do isat = 1, RCF%nprn
    if (rmstbl(isat) .gt. 0 .and. RCF%prn(isat)(1:1) .eq. 'R') then
      write (line(j + 1:), '(a3,a1)') RCF%prn(isat), ' '
      write (line1(j + 1:), '(i3,a1)') nint(rmstbl(isat)), ' '
      j = j + 4
    endif
  enddo
  if (line .ne. '') then
    write (*, '(a)') trim(line)
    write (*, '(a)') trim(line1)
  endif
  isat = 1
  line = ' '
  line1 = ' '
  j = 0
  do isat = 1, RCF%nprn
    if (rmstbl(isat) .gt. 0 .and. RCF%prn(isat)(1:1) .eq. 'E') then
      write (line(j + 1:), '(a3,a1)') RCF%prn(isat), ' '
      write (line1(j + 1:), '(i3,a1)') nint(rmstbl(isat)), ' '
      j = j + 4
    endif
  enddo
  if (line .ne. '') then
    write (*, '(a)') trim(line)
    write (*, '(a)') trim(line1)
  endif
  isat = 1
  line = ' '
  line1 = ' '
  j = 0
  do isat = 1, RCF%nprn
    if (rmstbl(isat) .gt. 0 .and. RCF%prn(isat)(1:1) .eq. 'C') then
      write (line(j + 1:), '(a3,a1)') RCF%prn(isat), ' '
      write (line1(j + 1:), '(i3,a1)') nint(rmstbl(isat)), ' '
      j = j + 4
    endif
  enddo
  if (line .ne. '') then
    write (*, '(a)') trim(line)
    write (*, '(a)') trim(line1)
  endif
  isat = 1
  line = ' '
  line1 = ' '
  j = 0
  do isat = 1, RCF%nprn
    if (rmstbl(isat) .gt. 0 .and. RCF%prn(isat)(1:1) .eq. 'J') then
      write (line(j + 1:), '(a3,a1)') RCF%prn(isat), ' '
      write (line1(j + 1:), '(i3,a1)') nint(rmstbl(isat)), ' '
      j = j + 4
    endif
  enddo
  if (line .ne. '') then
    write (*, '(a)') trim(line)
    write (*, '(a)') trim(line1)
  endif
  write (*, '(a)') '%%%-RMS OF RESIDUALS---PHASE(MM)'

  return
100 write (*, '(2a)') '***ERROR(read_residual): read file ', trim(line)
  call exit(1)
end
