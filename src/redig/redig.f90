!
!! redig.f90
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
!! residual diagnose
!
program redig
  implicit none
  include '../header/const.h'
  include 'rescfg.h'

!
!! common
  integer*4 idxfrq(MAXSYS, 2)
  common idxfrq
!
!! parameter
  logical*1 again
  integer*4 i, j, k, nepo, iepo, isat, jdb, iy, imon, id, ih, im, ierr
  character*3 iprn
  integer*4 avaobs, remobs, maxamb, totamb, oldam, newam, oldrem, newrem
  integer*4, pointer :: xamb(:), flag(:, :)
  real*8 tb, sec, intv
  real*8, pointer :: resi(:, :), resp(:, :)
  character*30, pointer :: rhd(:, :)
  character*27, pointer :: trsi(:)
  character*256 line, str(MAXSAT)
  type(rescfg) RCF
!
!! function called
  logical*1 istrue
  integer*4 modified_julday, pointer_string, get_valid_unit, find_ambd
  real*8 timdif, cutemod
!
!! read configure arguments
  call get_redig_args(nepo, RCF)
!
!! allocate memory
  allocate (xamb(nepo), stat=ierr)
  if (ierr .ne. 0) then
    write (*, '(a)') '***ERROR(redig): memory allocation xamb'
    call exit(1)
  end if
  allocate (rhd(RCF%nprn, nepo), stat=ierr)
  if (ierr .ne. 0) then
    write (*, '(a)') '***ERROR(redig): memory allocation rhd'
    call exit(1)
  end if
  allocate (resi(nepo, RCF%nprn), stat=ierr)
  if (ierr .ne. 0) then
    write (*, '(a)') '***ERROR(redig): memory allocation resi'
    call exit(1)
  end if
  allocate (flag(nepo, RCF%nprn), stat=ierr)
  if (ierr .ne. 0) then
    write (*, '(a)') '***ERROR(redig): memory allocation flag'
    call exit(1)
  end if
  allocate (trsi(nepo), stat=ierr)
  if (ierr .ne. 0) then
    write (*, '(a)') '***ERROR(redig): memory allocation trsi'
    call exit(1)
  end if
  if (RCF%pcedit) then
    allocate (resp(nepo, RCF%nprn), stat=ierr)
    if (ierr .ne. 0) then
      write (*, '(a)') '***ERROR(redig): memory allocation resp'
      call exit(1)
    end if
  end if
!
!! read residuals
  write (*, '(a)') '$$$MESSAGE: READING RESIDUALS ...'
  call read_residual(nepo, resi, resp, flag, trsi, RCF)
!
!! edit residuals
  if (RCF%jump .ne. 0.d0) then
    write (*, '(a)') '$$$MESSAGE: EDITING RESIDUALS ...'
    write (*, '(a,a4)') '%%% SCREENING FOR SITE ', RCF%snam
    do isat = 1, RCF%nprn
      write (*, '(a,a3)') ' $$ FOR PRN ', RCF%prn(isat)
!
!! remove short before checking jump
      call remov_shrt(RCF%nsht, nepo, flag(1, isat), trsi)
      if (RCF%pcedit) call check_pcres(RCF%pjump, nepo, flag(1, isat), resp(1, isat), trsi)
!
!! check jump
      again = .true.
      do while (again)
        call check_slip(RCF, nepo, flag(1, isat), resi(1, isat), trsi, again)
      end do
      if (RCF%pcedit) call check_arcslip(RCF%jump, nepo, flag(1, isat), resi(1, isat), trsi, again)
!
!! remove short after checking jump
      call remov_shrt(RCF%nsht, nepo, flag(1, isat), trsi)
    end do
  else
    goto 50
  end if
!
!! update rinex diagnose file
  write (*, '(a)') '$$$MESSAGE: UPDATING LOG FILES ...'
!
!! initialization
  jdb = 0
  tb = 0.d0
  do iepo = 1, nepo
    xamb(iepo) = 0
    do isat = 1, RCF%nprn
      rhd(isat, iepo) = ' '
    end do
  end do
!
!! read original health diagnose info
  oldam = 0
  oldrem = 0
  line = ' '
  do while (.true.)
    read (RCF%lfnrhd, '(a)', iostat=ierr) line
    if (ierr .ne. 0) exit
    if (index(line, 'AMB MAX/TOT/NEW') .ne. 0) read (line(21:), *) oldam
    if (index(line, 'EPO AVA/REM/NEW') .ne. 0) read (line(21:), *) oldrem
    if (index(line, 'SYS / FREQUENCY BAND') .ne. 0) then
      do i = 1, MAXSYS
        if (line(1:3) .eq. GNSS_NAME_LEN3(i)) then
          if (line(11:13) .ne. FREQ_NAME_SYS(idxfrq(i, 1), i) .or. &
              line(16:18) .ne. FREQ_NAME_SYS(idxfrq(i, 2), i)) then
            write (*, '(a)') '***ERROR(redig): inconsistent frequency band'
            write (*, '(a,a3,7x,2(a3,2x))') 'resfil: ', GNSS_NAME_LEN3(i), &
              FREQ_NAME_SYS(idxfrq(i, 1), i), FREQ_NAME_SYS(idxfrq(i, 2), i)
            write (*, '(a,a3,7x,2(a3,2x))') 'logfil: ', line(1:3), line(11:13), line(16:18)
            call exit(1)
          end if
        end if
      end do
    end if
    if (line(1:3) .ne. 'TIM') cycle
    read (line, '(3x,i5,4i3,f11.7)') iy, imon, id, ih, im, sec
    jdb = modified_julday(id, imon, iy)
    tb = 3600.d0*ih + 60.d0*im + sec
    intv = timdif(jdb, tb, RCF%jd0, RCF%sod0)
    if (intv .lt. 0.d0 .or. cutemod(intv, RCF%dintv) .ne. 0.d0) cycle
    iepo = nint(intv/RCF%dintv) + 1
    if (iepo .gt. nepo) then
      write (*, '(a,i5,f10.2)') '%%%MESSAGE(redig): log file truncated at ', jdb, tb
      exit
    end if
    do while (.true.)
      read (RCF%lfnrhd, '(a)', iostat=ierr) line
      if (ierr .ne. 0) exit
      if (line(1:3) .eq. 'TIM') then
        backspace RCF%lfnrhd
        exit
      end if
      if (line(61:63) .eq. 'AMB') cycle
      read (line, *) iprn
      isat = pointer_string(RCF%nprn, RCF%prn, iprn)
      if (isat .eq. 0) cycle
      rhd(isat, iepo) = line(61:len_trim(line))
    end do
  end do
  close (RCF%lfnrhd)
!
!! set ambiguity flags in log
  newam = 0
  totamb = 0
  avaobs = 0
  remobs = 0
  newrem = 0
  do iepo = 1, nepo
    do isat = 1, RCF%nprn
      if (istrue(flag(iepo, isat), 'AMB')) then
        if (istrue(flag(iepo, isat), 'NEWAMB')) newam = newam + 1
        k = find_ambd(nepo, flag(1, isat), iepo)
        write (rhd(isat, iepo), '(a3,i7.7)', iostat=ierr) 'AMB', k
        if (ierr .ne. 0) then
          write (*, '(a,i7)') '***ERROR(redig): write rhd amb ', k
          call exit(1)
        end if
        do j = iepo, k
          xamb(j) = xamb(j) + 1
        end do
      else if (istrue(flag(iepo, isat), 'GOOD')) then
        avaobs = avaobs + 1
      else if (istrue(flag(iepo, isat), 'DELLOW')) then
        rhd(isat, iepo) = 'DEL_LOWELEVATION'
        newrem = newrem + 1
      else if (istrue(flag(iepo, isat), 'DELORB')) then
        rhd(isat, iepo) = 'DEL_NOORBIT'
        newrem = newrem + 1
      else if (istrue(flag(iepo, isat), 'DELCLK')) then
        rhd(isat, iepo) = 'DEL_NOCLOCK'
        newrem = newrem + 1
      else if (istrue(flag(iepo, isat), 'DELBIA')) then
        rhd(isat, iepo) = 'DEL_NOBIA'
        newrem = newrem + 1
      else if (istrue(flag(iepo, isat), 'DELSHT')) then
        rhd(isat, iepo) = 'DEL_SHORTPIECE'
        newrem = newrem + 1
      else if (.not. istrue(flag(iepo, isat), 'NODATA')) then
        rhd(isat, iepo) = 'DEL_OUTLIER'
        newrem = newrem + 1
      end if
      if (rhd(isat, iepo) (1:3) .eq. 'AMB') then
        totamb = totamb + 1
      else if (rhd(isat, iepo) (1:3) .eq. 'DEL') then
        remobs = remobs + 1
      end if
    end do
  end do
  maxamb = maxval(xamb, 1)
!
!! write new log
  if (newam .gt. 0 .or. newrem .gt. 0 .or. newam .ne. oldam .or. newrem .ne. oldrem) then
    RCF%lfnrhd = get_valid_unit(10)
    open (RCF%lfnrhd, file=RCF%flnrhd, status='old', iostat=ierr)
    write (RCF%lfnrhd, '(a20,6x,a4,30x,a)') 'RINEX health Logfile', RCF%snam, 'COMMENT'
    write (RCF%lfnrhd, '(f10.2,50x,a)') RCF%dintv, 'INTERVAL'
    write (RCF%lfnrhd, '(a15,45x,a)') RCF%amb_at_dbd, 'AMB AT DAY-BOUNDARY'
    write (RCF%lfnrhd, '(3i10,30x,a)') maxamb, totamb, newam, 'AMB MAX/TOT/NEW'
    write (RCF%lfnrhd, '(3i10,30x,a)') avaobs, remobs, newrem, 'EPO AVA/REM/NEW'
    call mjd2date(RCF%jd0, RCF%sod0, iy, imon, id, ih, im, sec)
    write (RCF%lfnrhd, '(i5,4i3,f11.7,f12.2,20x,a)') iy, imon, id, ih, im, sec, (nepo - 1)*RCF%dintv, 'RES TIME BEG/LEN'
    do i = 1, MAXSYS
      write (RCF%lfnrhd, '(a3,7x,2(a3,2x),40x,a)') GNSS_NAME_LEN3(i:i), &
        FREQ_NAME_SYS(idxfrq(i, 1), i), FREQ_NAME_SYS(idxfrq(i, 2), i), 'SYS / FREQUENCY BAND'
    end do
    write (RCF%lfnrhd, '(60x,a)') 'END OF HEADER'
    do iepo = 1, nepo
      i = 0
      do isat = 1, RCF%nprn
        if (rhd(isat, iepo) (1:3) .eq. 'AMB') then
          i = i + 1
          read (rhd(isat, iepo), '(3x,i7)', iostat=ierr) k
          if (ierr .ne. 0) then
            write (*, '(a)') '***ERROR(redig): read index '
            call exit(1)
          end if
          call timinc(RCF%jd0, RCF%sod0, (k - 1)*RCF%dintv, jdb, tb)
          call mjd2date(jdb, tb, iy, imon, id, ih, im, sec)
          str(i) = ' '
          write (str(i), '(a3,28x,i5,4i3,f11.7,1x,a)') RCF%prn(isat), iy, imon, id, ih, im, sec, 'AMB'
        else if (rhd(isat, iepo) (1:3) .eq. 'DEL') then
          i = i + 1
          str(i) = ' '
          write (str(i), '(a3,57x,a)') RCF%prn(isat), rhd(isat, iepo) (1:len_trim(rhd(isat, iepo)))
        end if
      end do
      if (i .eq. 0) cycle
      call timinc(RCF%jd0, RCF%sod0, (iepo - 1)*RCF%dintv, jdb, tb)
      call mjd2date(jdb, tb, iy, imon, id, ih, im, sec)
      write (RCF%lfnrhd, '(a3,i5,4i3,f11.7)') 'TIM', iy, imon, id, ih, im, sec
      do j = 1, i
        write (RCF%lfnrhd, '(a)') str(j) (1:len_trim(str(j)))
      end do
    end do
    close (RCF%lfnrhd)
  end if
!
!! total statistics
  write (*, '(a)') '$$$MESSAGE: TOTAL STATISTICS '
  write (*, '(a,i10,f8.2,a1,/,a,i10,f8.2,a1)') 'NEWLY REMOVED: ', newrem, newrem*1.d2/(avaobs + newrem), '%', &
                                               'NEWLY AMBIGUT: ', newam,   newam*1.d2/(totamb - newam),  '%'
!
!! release memory
50 deallocate (resi)
  if (RCF%pcedit) deallocate (resp)
  deallocate (flag)
  deallocate (rhd)
  deallocate (xamb)
  deallocate (trsi)
end program
