!
!! findkey.f90
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
!! Contributor: Maorong Ge
!! 
!!
!!
!! purpose  : get the content of the line start with `keyword` within the
!!            sinex_bracket or in the whole file if the bracket is empty.
!!             Line start with * or # is ignored as comment. Keyword and its
!!             content is seperated by `=`. Comment on the line is seperated
!!             by `!`.
!!
!! parameter: lfn  -- file unit
!!             keyword -- keyword
!!             findkey -- content of keyword. `EMPTY` is returned if the
!!                        keyword is not defined in the file
!!
!
character(*) function findkey(lfn, keyword, sinex_bracket)
  implicit none
  character(*) :: keyword, sinex_bracket
  integer*4 lfn
!
!! local
  logical*1 continous_line
  integer*4 i, j, k, l
  character*65536 line

  rewind lfn
  l = len_trim(sinex_bracket)
  line = ' '
  if (l .ne. 0) then
    do while (index(line, '+'//sinex_bracket(1:l)) .eq. 0)
      read (lfn, '(a)', end=100) line
    enddo
  endif
!
!! `++` at the end of the value, before `!`, means that continous line follows
  findkey = '++'
  continous_Line = .false.
  do while (findkey(len_trim(findkey) - 1:len_trim(findkey)) .eq. '++')
    read (lfn, '(a)', end=100) line
!
!! end of bracket
    if (l .ne. 0 .and. index(line, '-'//sinex_bracket(1:l)) .ne. 0) goto 100
!
!! comment lines
    if (line(1:1) .eq. '#' .or. line(1:1) .eq. '*') cycle
!
!! if keyword there. if not must be continous line
    j = len_trim(keyword)
    i = index(line, keyword(1:j))
!
!! if seperator '=' between keyword and its value there. if not from the beginning
    j = index(line, '=')
    if (j .eq. 0) then
      call left_justify_string(line)
      line = ' '//line
    endif
!
!! comment at the end, starts with '!', ignored
    k = index(line, '!') - 1
    if (i .eq. 0 .and. .not. continous_line) cycle
    i = len(findkey)
    if (k .le. j) k = i
    findkey(len_trim(findkey) - 1:) = line(j + 2:k)
    continous_line = .true.
  enddo
  return

100 continue
  findkey = 'EMPTY'
  write (*, '(a)') '***WARNING(findkey): `'//trim(keyword)//'` not found.'
  write (*, '(a)') '   bracket: ', trim(sinex_bracket)
  return
end
