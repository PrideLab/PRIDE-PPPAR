!
!! split_string.f90
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
!! Contributor: Maorong Ge
!! 
!!
!!
!!! purpose  : split string into words
!!
!! parameter: lnoempty --  ignore blank or empty words  (len_trim(word)==0)
!!            string   --  the string to be splitted
!!            c_start c_end --  starting and ending character, if len_trim(c)==0, the whole string.
!!            seperator  -- seperator
!!            nword -- total number of words
!!            word  -- words
!!
!
subroutine split_string(lnoempty, string, c_start, c_end, seperator, nword, word)
  implicit none
  include '../header/const.h'

  logical*1 lnoempty
  character*(*) string, word(1:*)
  character*1 c_start, c_end, seperator
  integer*4 nword
!
!! local
  integer*4 i0, i1, ilast, i
  character*256 line

  line = string
  if (len(string) .gt. 256) then
    write (*, *) '***ERROR(split_line): input string length > 256'
    call exit(1)
  endif

  i0 = 0
  i1 = len_trim(line) + 1
  if (len_trim(c_start) .ne. 0) i0 = index(line, c_start)
  if (len_trim(c_end) .ne. 0) i1 = index(line, c_end)

  nword = 0
  if (i1 .eq. 0 .or. i0 .gt. i1) return

  ilast = i0 + 1
  line(i1:i1) = seperator
  do i = i0 + 1, i1
    if (line(i:i) .eq. seperator) then
      if (len_trim(line(ilast:i - 1)) .ne. 0) then
        nword = nword + 1
        word(nword) = line(ilast:i - 1)
      else
        if (.not. lnoempty) then
          nword = nword + 1
          word(nword) = ' '
        endif
      endif
      ilast = i + 1
    endif
  enddo

  return
end
