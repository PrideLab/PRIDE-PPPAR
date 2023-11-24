!
!! file_name.f90
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
!! Contributor: Maorong Ge, Jianghui Geng, Songfeng Yang
!! 
!!
!!
!!! purpose  : get the real file name based on its format definition input as variable or
!!            from 'file_name' file. In the file definition variables are indicated
!!            by a pair of `-`, for example -YYYY- for four digital year.
!!
!! parameter: ldefined -- format is defined by `name` instead of that from `file_name`
!!            keyword  -- short indentify for the file.
!!            param_list -- values of variables for file name, except time.
!!               variable_name=variable_value:variable_name=variable_value[:  ]
!!               SNAM=pots:cen=gfz:SATNAM=champ
!!            iyear,imonth,iday ihour -- time
!!            name -- output file name, input file format if ldefined is true
!!
!! warning  : short indentify for the file and variable_name must be defined fixed for
!!            programming.
!!
!
subroutine file_name(ldefined, keyword, param_list, iyear, imonth, iday, ihour, name)
  implicit none

  logical*1 ldefined
  integer*4 iyear, imonth, iday, ihour
  character*(*) keyword, param_list, name
!
!! local
  logical*1 first, ltable, found
  character*256 form
  integer*4 i, j, k, l, mjd, week, wd, idoy, ierr, lun
!
!! default filename definition
  integer*4 nfile, npar, nword
  integer*4, parameter :: MAXFIL = 100
  character*10 fname(MAXFIL)
  character*30 param(40)
  character*256 word(40)
  character*256 fform(MAXFIL)
!
!! function called
  integer*4 modified_julday, get_valid_unit, day_of_year

  data first, ltable/.true., .true./
  save first, ltable, nfile, fname, fform
!
!! read in file name table
  if (first) then
    first = .false.
    lun = get_valid_unit(10)
    open (lun, file='file_name', status='OLD', iostat=ierr)
    nfile = 0
    if (ierr .ne. 0) ltable = .false.
    if (ltable) then
      do while (.true.)
        read (lun, '(a10,1x,a)', end=100) fname(nfile + 1), fform(nfile + 1)
        nfile = nfile + 1
        call left_justify_string(fform(nfile + 1))
        if (nfile .gt. MAXFIL) then
          write (*, '(a,2i4)') '***ERROR(file_name): too manay files &
            defined in `file_name`', nfile, MAXFIL
          call exit(1)
        endif
      enddo
    endif
100 close (lun)
  endif
!
!! first find out the format
  if (ldefined) then
    form = name
  else
    form = ' '
    if (ltable) then
      do i = 1, nfile
        if (fname(i) (1:len_trim(fname(i))) .eq. &
            trim(keyword)) then
          form = fform(i)
          exit
        endif
      enddo
    endif
    if (form(1:1) .eq. ' ') then
      write (*, '(2a)') '***ERROR(file_name): not defined in &
        `file_name`, ', trim(keyword)
      call exit(1)
    endif
  endif
!
!! find out the variables in the format
  call split_string(.false., form, ' ', ' ', '-', nword, word)
  call split_string(.true., param_list, ' ', ' ', ':', npar, param)
!
!! name definition
  mjd = modified_julday(iday, imonth, iyear)
  call gpsweek(iday, imonth, iyear, week, wd)
  idoy = day_of_year(iday, imonth, iyear)
  write (param(npar + 1), '(a,i4.4)') 'YYYY=', iyear
  write (param(npar + 2), '(a,i3.3)') 'DDD=', idoy
  write (param(npar + 3), '(a,i2.2)') 'HH=', ihour
  write (param(npar + 4), '(a,i4.4)') 'WWWW=', week
  write (param(npar + 5), '(a,i5.5)') 'WWWWD=', week*10 + wd
  write (param(npar + 6), '(a,i2.2)') 'YY=', iyear - int(iyear/100)*100
  write (param(npar + 7), '(a,i1.1)') 'Y=', iyear - int(iyear/10)*10
  write (param(npar + 8), '(a,i5.5)') 'MJD=', mjd
  npar = npar + 8

  name = ' '
  do i = 1, nword
    if (len_trim(word(i)) .eq. 0) cycle
    if (int(i/2)*2 .eq. i) then
      l = len_trim(word(i))
      found = .false.
      do j = 1, npar
        if (word(i) (1:l)//'=' .eq. param(j) (1:l + 1)) then
          k = len_trim(name)
          name(k + 1:) = param(j) (l + 2:)
          found = .true.
          exit
        endif
      enddo
      if (.not. found) then
        write (*, '(2a)') '***ERROR(file_name): variable not defined,', word(i) (1:l)
        call exit(1)
      endif
    else
      l = len_trim(name)
      name(l + 1:) = word(i)
    endif
  enddo

  return
end
