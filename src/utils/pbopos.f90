!
!! pbopos.f90
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
!! Contributor: Jihang Lin, Yuanxin Pan
!!
!!
program pbopos
  implicit none
!
!! variable list
  !! constant
  real*8, parameter :: RAD2DEG = 45.d0/atan(1.d0)
  !! global
  integer*4     narg, nfil, nrec
  integer*4     ierr
  integer*4     ityp
  !! file
  character*4   site, lsit, usit
  integer*4     lfntmp, lfnpos, lfnpbo
  !! record
  integer*4     mjd, year, doy, imon, iday, ih, im
  real*8        fjd, sod, sec
  real*8        refxyz(3), avgxyz(3), tmpxyz(3)
  real*8        stdxyz(3), covxyz(3), difxyz(3)
  real*8        refblh(3), tmpblh(3)
  real*8        stdenu(3), covenu(3)
  real*8        tmpenu(3)
  real*8        sigma0
  real*8        rotmat(3,3), tsqrot(3,3)
  real*8        mqqxyz(3,3), mqqenu(3,3)
  character*15  first_epoch, last_epoch
  character*15  current_epoch
  character*5   solution_type
  character*3   analcent_type
  !! tmporary
  integer*4     i, j
  character*3   filtyp
  character*256 filnam
  character*512 line
!
!! function list
  integer*4 get_valid_unit
!
!! read argument
  narg = iargc()
  if (narg .lt. 2) then
    write (*, '(a)') 'usage: pbopos site path [x_ref y_ref z_ref]'
    write (*, '(a)') ''
    write (*, '(a)') '  convert PRIDE-PPPAR pos files to PBO position series'
    write (*, '(a)') '  created on JAN-11, 2022'
    write (*, '(a)') ''
    write (*, '(a)') 'example:'
    write (*, '(a)') '  1. pbopos alic pos_2020001_alic'
    write (*, '(a)') '  2. pbopos alic 2021/001/pos_2020001_alic'
    write (*, '(a)') ''
    write (*, '(a)') 'notice:'
    write (*, '(a)') '    all position files with standard naming in directory will be recognized '
    write (*, '(a)') '  and found automatically, depends on which kind of path you input: '
    write (*, '(a)') ''
    write (*, '(a)') '  1. ''./''             stop untill no successive pos file exists'
    write (*, '(a)') '  2. ''./yyyy/ddd/      stop untill no successive year folder exists'
    write (*, '(a)') ''
    write (*, '(a)') '    after that, all results will be output into one single PBO file.'
    write (*, '(a)') ''
    write (*, '(a)') '  naming convetion of input files: pos_yyyyddd_site'
    write (*, '(a)') '  naming convetion of output file: SITE.aaa.ttttt_igs14.pos'
    write (*, '(a)') ''
    write (*, '(a)') '  yyyy       4-digit year'
    write (*, '(a)') '  ddd        3-digit day of year'
    write (*, '(a)') '  site/SITE  4-character lower/upper-case site ID'
    write (*, '(a)') '  aaa        3-character analysis center ID'
    write (*, '(a)') '  ttttt      5-character solution type ID'
    call exit(1)
  endif
  call getarg(1, site)
  call getarg(2, filnam)

  refxyz = 0.d0
  refblh = 0.d0
  if (narg .ge. 3) then
    do i = 3, 5
      call getarg(i, line)
      read (line, *, iostat=ierr) refxyz(i-2)
      if (ierr .ne. 0) then
        write (*, '(a,i1,1x,a)') '***ERROR: read argument ', trim(line)
        call exit(1)
      endif
    enddo
    call xyzblh(refxyz, 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, refblh)
  endif
!
!! case conversion
  do i = 1, len(site)
    j = ichar(site(i:i))
    if (j >= 65 .and. j <= 90) then
       lsit(i:i) = char(j+32)
       usit(i:i) = char(j)
    elseif (j >=97 .and. j <= 122) then
       lsit(i:i) = char(j)
       usit(i:i) = char(j-32)
    else
       lsit(i:i) = char(j)
       usit(i:i) = char(j)
    endif
  enddo

  lfnpos = get_valid_unit(10)
  open(lfnpos, file=filnam, status='old', iostat=ierr)
  if (ierr .ne. 0) goto 100

  lfntmp = get_valid_unit(10)
  open(lfntmp, file='pbopos.tmp', form='unformatted')
!
!! read pos files
  nrec = 0
  first_epoch = ''
  last_epoch  = ''
  analcent_type = 'wum'
  solution_type = 'rapid'
  avgxyz = 0.d0
  do while (.true.)
    inquire(unit=lfnpos, name=filnam)
    do while (.true.)
      read (lfnpos, '(a)', end=10, err=200) line
      if (trim(line(61:)) .eq. 'STATION') then
        if (line(1:4) .eq. usit .or. line(1:4) .eq. lsit) cycle
        write (*, '(a,2(1x,a))') '###WARNING: conflicting site name:', line(1:4), filnam
        goto 10
      endif
      if (trim(line(61:)) .eq. 'SAT ORBIT') then
        do i = 1, 3
            j = ichar(line(i:i))
            if (j >= 65 .and. j <= 90) then
              analcent_type(i:i) = char(j+32)
            else
              analcent_type(i:i) = line(i:i)
            endif
        enddo
        select case (line(8:10))
          case ('FIN')
            solution_type = 'final'
          case ('RAP')
            solution_type = 'rapid'
          case ('RTS')
            solution_type = 'suppl'
          case ('ULA')
            solution_type = 'suppf'
        endselect
        select case (line(5:7))
          case ('R03')
            solution_type = 'repro'
          case ('R3T')
            solution_type = 'repro'
        endselect
      endif
      if (trim(line(61:)) .eq. 'POS MODE/PRIORI (meter)') then
        if (line(1:6) .eq. 'Static')    filtyp = 'pos'
      endif
      if (trim(line(61:)) .eq. 'END OF HEADER') then
        do while (.true.)
          read (lfnpos, '(a)', end=10, err=200) line
          if (line(1:1) .eq. '*') cycle
          if (filtyp .eq. 'pos') then
            read (line, '(1x,a4,1x,f11.4,3f15.4,7e25.14)', err=250) &
              site, fjd, tmpxyz(1:3), stdxyz(1:3), covxyz(1:3), sigma0
            if (site .ne. usit .and. site .ne. lsit) cycle
            mjd = int(fjd)
            sod = (fjd - mjd) * 86400.d0
          endif
          write (lfntmp) mjd, sod, tmpxyz, stdxyz, covxyz, sigma0
          if (len_trim(first_epoch) .eq. 0) then
            call mjd2date(mjd, sod, year, imon, iday, ih, im, sec)
            write (first_epoch, '(i4,i0.2,i0.2,1x,i0.2,i0.2,i0.2)') &
              year, imon, iday, ih, im, nint(sec)
            if (all(refxyz .eq. 0.d0)) refxyz = tmpxyz
          endif
          avgxyz = avgxyz + tmpxyz - refxyz
          nrec = nrec + 1
        enddo
      endif
    enddo
10 continue
    call next_posfil(lfnpos, j)
    close(lfnpos)
    if (j .eq. 0) exit
    lfnpos = j
  enddo

  avgxyz = avgxyz/nrec
  avgxyz = refxyz + avgxyz
  if (all(refblh .eq. 0.d0)) then
    refxyz = avgxyz
    call xyzblh(refxyz, 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, refblh)
  endif

  call mjd2date(mjd, sod, year, imon, iday, ih, im, sec)
  write (last_epoch, '(i4,i0.2,i0.2,1x,i0.2,i0.2,i0.2)') &
    year, imon, iday, ih, im, nint(sec)

  call date_and_time(current_epoch, line)
  current_epoch(10:15) = line(1:6)
!
!! write PBO file
  filnam = usit//'.'//analcent_type//'.'//solution_type//'_igs14.pos'
  lfnpbo = get_valid_unit(10)
  open(lfnpbo, file=filnam, status='replace', err=100)
  !! header
  write (lfnpbo, '(a)') 'PBO Station Position Time Series. Reference Frame : ITRF2014'
  write (lfnpbo, '(a)') 'Format Version: 1.1.0'
  write (lfnpbo, '(2a)') '4-character ID: ', usit
  write (lfnpbo, '(2a)') 'Station name  : ', usit
  write (lfnpbo, '(2a)') 'First Epoch   : ', first_epoch
  write (lfnpbo, '(2a)') 'Last Epoch    : ', last_epoch
  write (lfnpbo, '(2a)') 'Release Date  : ', current_epoch
  write (lfnpbo, '(a,3f15.5, a)') 'XYZ Reference position : ', refxyz(1:3), ' (ITRF2014)' 
  write (lfnpbo, '(a,2f16.10,f11.5,a)') 'NEU Reference position : ', refblh(1:2)*RAD2DEG, refblh(3), ' (ITRF2014/WGS84)' 
  write (lfnpbo, '(a)') 'Start Field Description'
  write (lfnpbo, '(a)') 'YYYYMMDD      Year, month, day for the given position epoch'
  write (lfnpbo, '(a)') 'HHMMSS        Hour, minute, second for the given position epoch'
  write (lfnpbo, '(a)') 'JJJJJ.JJJJJ   Modified Julian day for the given position epoch'
  write (lfnpbo, '(a)') 'X             X coordinate, Specified Reference Frame, meters'
  write (lfnpbo, '(a)') 'Y             Y coordinate, Specified Reference Frame, meters'
  write (lfnpbo, '(a)') 'Z             Z coordinate, Specified Reference Frame, meters'
  write (lfnpbo, '(a)') 'Sx            Standard deviation of the X position, meters'
  write (lfnpbo, '(a)') 'Sy            Standard deviation of the Y position, meters'
  write (lfnpbo, '(a)') 'Sz            Standard deviation of the Z position, meters'
  write (lfnpbo, '(a)') 'Rxy           Correlation of the X and Y position'
  write (lfnpbo, '(a)') 'Rxz           Correlation of the X and Z position'
  write (lfnpbo, '(a)') 'Ryz           Correlation of the Y and Z position'
  write (lfnpbo, '(a)') 'Nlat          North latitude, WGS-84 ellipsoid, decimal degrees'
  write (lfnpbo, '(a)') 'Elong         East longitude, WGS-84 ellipsoid, decimal degrees'
  write (lfnpbo, '(a)') 'Height (Up)   Height relative to WGS-84 ellipsoid, m'
  write (lfnpbo, '(a)') 'dN            Difference in North component from NEU reference position, meters'
  write (lfnpbo, '(a)') 'dE            Difference in East component from NEU reference position, meters'
  write (lfnpbo, '(a)') 'du            Difference in vertical component from NEU reference position, meters'
  write (lfnpbo, '(a)') 'Sn            Standard deviation of dN, meters'
  write (lfnpbo, '(a)') 'Se            Standard deviation of dE, meters'
  write (lfnpbo, '(a)') 'Su            Standard deviation of dU, meters'
  write (lfnpbo, '(a)') 'Rne           Correlation of dN and dE'
  write (lfnpbo, '(a)') 'Rnu           Correlation of dN and dU'
  write (lfnpbo, '(a)') 'Reu           Correlation of dE and dU'
  write (lfnpbo, '(a)') 'Soln          "rapid", "final", "suppl/suppf", "campd", or "repro" corresponding to products '// &
                        'generated with rapid or final orbit products, in supplemental processing, campaign data processing'
  write (lfnpbo, '(a)') 'End Field Description'
  write (lfnpbo, '(a)') '*YYYYMMDD HHMMSS JJJJJ.JJJJ         X             Y             Z            Sx        Sy    '// &
                        '   Sz     Rxy   Rxz    Ryz            NLat         Elong         Height         dN        dE '// &
                        '       dU         Sn       Se       Su      Rne    Rnu    Reu  Soln'
  !! record
  rewind lfntmp
  do while (.true.)
    read (lfntmp, end=20) mjd, sod, tmpxyz, stdxyz, covxyz, sigma0
    call mjd2date(mjd, sod, year, imon, iday, ih, im, sec)
    write (current_epoch, '(i4,i0.2,i0.2,1x,i0.2,i0.2,i0.2)') &
      year, imon, iday, ih, im, nint(sec)
    fjd = mjd + sod/86400.d0
    !! xyz to enu
    call xyzblh(tmpxyz, 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, tmpblh)
    call rot_enu2xyz(tmpblh(1), tmpblh(2), tsqrot)
    rotmat = transpose(tsqrot)
    difxyz = tmpxyz - refxyz
    call matmpy(rotmat, difxyz, tmpenu, 3,3,1)
    do i = 1, 3
      mqqxyz(i,i) = stdxyz(i)
    enddo
    mqqxyz(1,2) = covxyz(1)
    mqqxyz(1,3) = covxyz(2)
    mqqxyz(2,3) = covxyz(3)
    mqqxyz(2,1) = mqqxyz(1,2)
    mqqxyz(3,1) = mqqxyz(1,3)
    mqqxyz(3,2) = mqqxyz(2,3)
    call matmpy(rotmat, mqqxyz, mqqenu, 3, 3, 3)
    call matmpy(mqqenu, tsqrot, mqqenu, 3, 3, 3)
    !! cov to corr
    stdxyz = sqrt(stdxyz) * sigma0
    covxyz(1) = covxyz(1)/stdxyz(1)/stdxyz(2)
    covxyz(2) = covxyz(2)/stdxyz(1)/stdxyz(3)
    covxyz(3) = covxyz(3)/stdxyz(2)/stdxyz(3)
    do i = 1, 3
      stdenu(i) = sqrt(mqqenu(i,i)) * sigma0
    enddo
    covenu(1) = mqqenu(1,2)/stdenu(1)/stdenu(2)
    covenu(2) = mqqenu(1,3)/stdenu(1)/stdenu(3)
    covenu(3) = mqqenu(2,3)/stdenu(2)/stdenu(3)
    write (lfnpbo, '(1x,a15,1x,f10.4,3(1x,f14.5),3(1x,f8.5),3(1x,f6.3),3x,'// &
      '2(1x,f15.10),1x,f10.5,2x,3(1x,f9.5),2x,3(1x,f8.5),3(1x,f6.3),1x,a)')   &
      current_epoch, fjd, tmpxyz(1:3), stdxyz(1:3), covxyz(1:3),              &
      tmpblh(1:2)*RAD2DEG, tmpblh(3), tmpenu(2), tmpenu(1), tmpenu(3),        &
      stdenu(2), stdenu(1), stdenu(3), covenu(1), covenu(3), covenu(2),       &
      solution_type
  enddo

20 continue
  close(lfntmp, status='delete')
  close(lfnpbo)

  return

100 continue
  write (*, '(2a)') '***ERROR: open file ', trim(filnam)
  call exit(1)

200 continue
  write (*, '(2a)') '***ERROR: read file ', trim(filnam)
  write (*, '(a)') trim(line)
  call exit(1)

250 continue
  write (*, '(a)') '***ERROR: read line '
  write (*, '(a)') trim(line)
  call exit(1)

contains

  subroutine next_posfil(lfn_this, lfn_next)
    integer*4 lfn_this, lfn_next
  !
  !! local
    character*256 root_dir
    character*256 mid_dir
    character*256 name_this
    character*256 name_next
    character*256 tmp_name
    integer*4 nlen, len0, len1, len2, len3
    integer*4 doy0, year0, mjd0
    integer*4 doy1, year1, mjd1
    integer*4 ierr
    logical*4 is_alive
  !
  !! function used
    integer*4  get_valid_unit
    integer*4  modified_julday
  
    lfn_next = 0
    inquire(unit=lfn_this, name=name_this, exist=is_alive)
400 continue
    if (name_this .eq. '') then
      write (*, *) '###WARNING(next_posfil): open ', trim(name_this)
      return
    endif
  !
  !! split directory and base
    len0 = len_trim(name_this)
    len1 = index(name_this, '/',  BACK = .TRUE.)
    len2 = index(name_this, '\\', BACK = .TRUE.)
    if (len1 .ne. 0) then
      len2 = index(name_this(1:len1-1), '/',  BACK = .TRUE.)
      len3 = index(name_this(1:len2-1), '/',  BACK = .TRUE.)
      root_dir = trim(name_this(1:len1))
      tmp_name = trim(name_this(len1+1:)) 
      mid_dir  = ''
      if (len2-len3 .eq. 5 .and. (len1-len2 .eq. 4 .or. len1-len2 .eq. 8)) then
        mid_dir = root_dir(len3+1:len1)
        if (mid_dir(1:4) .eq. tmp_name(5:8) .and. &
            mid_dir(6:8) .eq. tmp_name(9:11)) then
          root_dir = trim(name_this(1:len3))
        else
          mid_dir = ''
        endif
      endif
      name_this = tmp_name
    elseif (len2 .ne. 0) then
      len1 = len2
      len2 = index(name_this(1:len1-2), '\\', BACK = .TRUE.)
      len3 = index(name_this(1:len2-2), '\\', BACK = .TRUE.)
      root_dir = trim(name_this(1:len1+1))
      tmp_name = trim(name_this(len1+1:)) 
      mid_dir  = ''
      if (len2-len3 .eq. 6 .and. (len1-len2 .eq. 5 .or. len1-len2 .eq. 9)) then
        mid_dir = root_dir(len3+2:len1+1)
        if (mid_dir(1:4) .eq. tmp_name(5:8) .and. &
            mid_dir(7:9) .eq. tmp_name(9:11)) then
          root_dir = trim(name_this(1:len3+1))
        else
          mid_dir = ''
        endif
      endif
      name_this = tmp_name
    else
      root_dir = ''
    endif
  !
  !! try writing successor file name
    nlen = len_trim(name_this)
    if (nlen .ne. 11 .and. nlen .ne. 16) then
      write (*,*) '###WARNING(next_posfil): unrecognized naming convention of pos files: ', trim(name_this)
      return
    endif

    read (name_this(1:nlen), '(4x,i4,i3)', iostat=ierr) year0, doy0
    if (ierr .ne. 0) then
      write (*,*) '###WARNING(next_posfil): illegal naming of pos files: ', trim(name_this)
      return
    endif

    mjd0 = modified_julday(doy0, 0, year0)
    call mjd2doy(mjd0+1, year0, doy0)

    if (len_trim(mid_dir) .gt. 0) then
      select case (len_trim(mid_dir))
        case (9)
          write (mid_dir, '(i4,a1,i0.3,a1)') year0, '/',  doy0, '/' 
        case (11)
          write (mid_dir, '(i4,a2,i0.3,a2)') year0, '\\', doy0, '\\' 
        case (13)
          read (mid_dir(10:12), '(i3)', iostat=ierr) doy1
          if (ierr .eq. 0) then
            year1 = year0
            if (doy1 .lt. doy0) year1 = year0 + 1
            mjd1 = modified_julday(doy1, 0, year1)
            call mjd2doy(mjd1+1, year1, doy1)
            write (mid_dir, '(i4,a1,i0.3,a1,i0.3,a1)') year0, '/',  doy0, '-', doy1, '/' 
          endif
        case (15)
          read (mid_dir(11:13), '(i3)', iostat=ierr) doy1
          if (ierr .eq. 0) then
            year1 = year0
            if (doy1 .lt. doy0) year1 = year0 + 1
            mjd1 = modified_julday(doy1, 0, year1)
            call mjd2doy(mjd1+1, year1, doy1)
            write (mid_dir, '(i4,a2,i0.3,a1,i0.3,a2)') year0, '\\', doy0, '-', doy1, '\\' 
          endif
      endselect
    endif

    write (name_next, '(a4,i4,i0.3,a)') name_this(1:4), year0, doy0, name_this(12:nlen)
    name_next = trim(root_dir)//trim(mid_dir)//trim(name_next)

    lfn_next = get_valid_unit(10)
    open(lfn_next, file=name_next, status='old', iostat=ierr)
    if (ierr .ne. 0) then
      lfn_next = 0
      tmp_name = ''
      len1 = index(name_next, '/',  BACK = .TRUE.)
      len2 = index(name_next, '\\', BACK = .TRUE.)
      if (len1 .gt. 0 ) tmp_name = trim(root_dir)//mid_dir(1:4)//'/.'
      if (len2 .gt. 0 ) tmp_name = trim(root_dir)//mid_dir(1:4)//'\\.'
      inquire(file=trim(tmp_name), exist=is_alive)
      if (is_alive) then
        name_this = name_next
        goto 400
      endif
    endif

    return
  end subroutine
    
end program
