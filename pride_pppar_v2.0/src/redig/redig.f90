!
!! redig.f90
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
!! residual diagnose
!
program redig
  implicit none
  include '../header/const.h'
  include 'rescfg.h'

  logical*1 again
  integer*4 i, j, k, nepo, iepo, isat, jdb, jdx, iy, imon, id, ih, im, ierr, suitepo
  character*3 iprn
  integer*4 avaobs, remobs, maxamb, totamb, xamb(MAXEPO), newam, newrem, rtot, atot
  integer*4, pointer :: flag(:, :)
  real*8 tb, tx, sec, intv
  real*8, pointer :: resi(:, :)
  character*20, pointer :: rhd(:, :)
  character*27 trsi(MAXEPO)
  character*256 line, str(MAXSAT)
  type(rescfg) RCF
!
!! function called
  logical*1 istrue
  integer*4 modified_julday, pointer_string, find_ambd
  real*8 timdif
!
!! instruction
  write (*, '(a)') '++++++++++++++++++++++++++++++++++++'
  write (*, '(a)') 'DIAGNOSE CARRIER-PHASE RESIDUALS'
  write (*, '(a)') '++++++++++++++++++++++++++++++++++++'
!
!! read configure arguments
  call get_redig_args(nepo, RCF)
!
!! allocate memory
  suitepo = nint(RCF%seslen/RCF%dintv)+1
  allocate (rhd(MAXSAT, suitepo), stat=ierr)
  if (ierr .ne. 0) then
    write (*, '(a)') '***ERROR(redig): memory allocation rhd'
    call exit(1)
  endif
  allocate (resi(nepo, RCF%nprn), stat=ierr)
  if (ierr .ne. 0) then
    write (*, '(a)') '***ERROR(redig): memory allocation resi'
    call exit(1)
  endif
  allocate (flag(nepo, RCF%nprn), stat=ierr)
  if (ierr .ne. 0) then
    write (*, '(a)') '***ERROR(redig): memory allocation flag'
    call exit(1)
  endif
!
!! read residuals
  write (*, '(a)') '$$$MESSAGE: READING RESIDUALS ...'
  call read_residual(nepo, resi, flag, trsi, RCF)
!
!! edit residuals
  if (RCF%jump .ne. 0.d0) then
    write (*, '(a)') '$$$MESSAGE: EDITING RESIDUALS ...'
    write (*, '(a,a4)') '%%% SCREENING FOR SITE ', RCF%snam
    do isat = 1, RCF%nprn
      write (*, '(a,a3)') ' $$ FOR PRN ', RCF%prn(isat)
!
!! remove huge residuals
      if (RCF%xres .ne. 0.d0) call delet_huge(RCF%lupd, RCF%xres, nepo, flag(1, isat), resi(1, isat), trsi)
!
!! remove short before checking jump
      call remov_shrt(RCF%lupd, RCF%nsht, nepo, flag(1, isat), trsi)
!
!! check jump
      again = .true.
      do while (again)
        call check_slip(RCF%lupd, RCF%jump, nepo, flag(1, isat), resi(1, isat), trsi, again)
      enddo
!
!! remove short after checking jump
      call remov_shrt(RCF%lupd, RCF%nsht, nepo, flag(1, isat), trsi)
    enddo
  endif
!
!! update rinex diagnose file
  rtot = 0; atot = 0
  write (*, '(a)') '$$$MESSAGE: UPDATING RHD FILES ...'
!
!! initialization
  jdb = 0
  tb = 0.d0
  do iepo = 1, suitepo
    xamb(iepo) = 0
    do isat = 1, MAXSAT
      rhd(isat, iepo) = ' '
    enddo
  enddo
!
!! read original health diagnose info
  line = ' '
  do while (line(1:3) .ne. 'TIM')
    read (RCF%lfnrhd, '(a)') line
    if (index(line, 'INT AMB/DEL') .ne. 0) read (line(11:), *) intv
  enddo
  read (line, '(3x,i5,4i3,f11.7)') iy, imon, id, ih, im, sec
  jdb = modified_julday(id, imon, iy)
  tb = ih*3600.d0 + im*60.d0 + sec
  backspace RCF%lfnrhd
  do while (.true.)
    read (RCF%lfnrhd, '(a)', iostat=ierr) line
    if (ierr .ne. 0) exit
    if (line(1:3) .ne. 'TIM') cycle
    read (line, '(3x,i5,4i3,f11.7)') iy, imon, id, ih, im, sec
    jdx = modified_julday(id, imon, iy)
    tx = ih*3600.d0 + im*60.d0 + sec
    iepo = nint(timdif(jdx, tx, jdb, tb)/intv) + 1
    if (iepo .gt. suitepo) goto 100
    do while (.true.)
      read (RCF%lfnrhd, '(a)', iostat=ierr) line
      if (ierr .ne. 0) exit
      if (line(1:3) .eq. 'TIM') then
        backspace RCF%lfnrhd
        exit
      endif
      if (line(61:63) .eq. 'AMB') cycle
      read (line, *) iprn
      isat = pointer_string(RCF%nprn, RCF%prn, iprn)
      if (isat .eq. 0) cycle
      rhd(isat, iepo) = line(61:len_trim(line))
    enddo
  enddo
!
!! set ambiguity flags in rhd
  newam = 0
  totamb = 0
  avaobs = 0
  remobs = 0
  newrem = 0
  do iepo = 1, nepo
    call timinc(RCF%jd0, RCF%sod0, (iepo - 1)*RCF%dintv, jdx, tx)
    i = nint(timdif(jdx, tx, jdb, tb)/intv) + 1
    if (i .gt. suitepo) goto 100
    if (i .lt. 1) cycle
    do isat = 1, RCF%nprn
      if (istrue(flag(iepo, isat), 'AMB')) then
        if (istrue(flag(iepo, isat), 'NEWAMB')) newam = newam + 1
        k = find_ambd(nepo, flag(1, isat), iepo)
        write (rhd(isat, i), '(a3,i7.7)') 'AMB', k
        do j = iepo, k
          xamb(j) = xamb(j) + 1
        enddo
      else if (istrue(flag(iepo, isat), 'GOOD')) then
        avaobs = avaobs + 1
      else if (.not. istrue(flag(iepo, isat), 'NODATA')) then
        rhd(isat, i) = 'DEL'
        newrem = newrem + 1
      endif
      if (rhd(isat, i) (1:3) .eq. 'AMB') then
        totamb = totamb + 1
      else if (rhd(isat, i) (1:3) .eq. 'DEL') then
        remobs = remobs + 1
      endif
    enddo
  enddo
  maxamb = maxval(xamb, 1)
  rtot = rtot + newrem
  atot = atot + newam
!
!! write new rhd
  rewind RCF%lfnrhd
  write (RCF%lfnrhd, '(a21,5x,a4,30x,a)') 'Rinex Health Diagnose', RCF%snam, 'COMMENT'
  write (RCF%lfnrhd, '(2f10.2,40x,a)') RCF%dintv, intv, 'INT AMB/DEL'
  write (RCF%lfnrhd, '(3i10,30x,a)') maxamb, totamb, newam, 'AMB MAX/TOT/NEW'
  write (RCF%lfnrhd, '(3i10,30x,a)') avaobs, remobs, newrem, 'EPO AVA/REM/NEW'
  write (RCF%lfnrhd, '(60x,a)') 'END OF HEADER'
  do iepo = 1, suitepo
    i = 0
    do isat = 1, RCF%nprn
      if (rhd(isat, iepo) (1:3) .eq. 'AMB') then
        i = i + 1
        read (rhd(isat, iepo), '(3x,i7)', err=200) k
        call timinc(RCF%jd0, RCF%sod0, (k - 1)*RCF%dintv, jdx, tx)
        call mjd2date(jdx, tx, iy, imon, id, ih, im, sec)
        write (str(i), '(a3,28x,i5,4i3,f11.7,1x,a)') RCF%prn(isat), iy, imon, id, ih, im, sec, 'AMB'
      else if (rhd(isat, iepo) (1:3) .eq. 'DEL') then
        i = i + 1
        write (str(i), '(a3,57x,a)') RCF%prn(isat), rhd(isat, iepo) (1:len_trim(rhd(isat, iepo)))
      endif
    enddo
    if (i .eq. 0) cycle
    tx = tb + (iepo - 1)*intv
    jdx = jdb + int(tx/86400.d0)
    tx = tx - (jdx - jdb)*86400.d0
    call mjd2date(jdx, tx, iy, imon, id, ih, im, sec)
    write (RCF%lfnrhd, '(a3,i5,4i3,f11.7)') 'TIM', iy, imon, id, ih, im, sec
    do j = 1, i
      write (RCF%lfnrhd, '(a)') str(j) (1:len_trim(str(j)))
    enddo
  enddo
  close (RCF%lfnrhd)
  RCF%lfnrhd = 0
!
!! total statistics
  write (*, '(a)') '$$$MESSAGE: TOTAL STATISTICS '
  write (*, '(a,i10,/,a,i10)') 'NEWLY REMOVED: ', rtot, 'NEWLY AMBIGUT: ', atot
!
!! release memory
50 deallocate (resi)
  deallocate (flag)
  deallocate (rhd)

  stop
100 write (*, '(a)') '***ERROR(redig): MAXEPO exceeded '
  call exit(1)
200 write (*, '(a)') '***ERROR(redig): read index '
  call exit(1)
end
