!
!! prn_matbld.f90
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
!! Contributor: Songfeng Yang
!! 
!!
!!
subroutine prn_matbld(prn_mat)
  implicit none
  include '../header/const.h'
  character*3 prn_mat(MAXSAT)
  integer*4 i
  do i=1,MAXSAT
    if(i.gt.0 .and. i.le.MAXSAT_G)then
      write(prn_mat(i),'(a1,i0.2)') 'G',i
    else if(i.gt.MAXSAT_G .and. i.le.MAXSAT_G+MAXSAT_R)then
      write(prn_mat(i),'(a1,i0.2)') 'R',i-MAXSAT_G
    else if(i.gt.MAXSAT_G+MAXSAT_R .and. i.le.MAXSAT_G+MAXSAT_R+MAXSAT_E)then
      write(prn_mat(i),'(a1,i0.2)') 'E',i-MAXSAT_G-MAXSAT_R
    else if(i.gt.MAXSAT_G+MAXSAT_R+MAXSAT_E .and. i.le.MAXSAT_G+MAXSAT_R+MAXSAT_E+MAXSAT_C)then
      write(prn_mat(i),'(a1,i0.2)') 'C',i-MAXSAT_G-MAXSAT_R-MAXSAT_E
    else if(i.gt.MAXSAT_G+MAXSAT_R+MAXSAT_E+MAXSAT_C .and. i.le.MAXSAT_G+MAXSAT_R+MAXSAT_E+MAXSAT_C+MAXSAT_J)then
      write(prn_mat(i),'(a1,i0.2)') 'J',i-MAXSAT_G-MAXSAT_R-MAXSAT_E-MAXSAT_C
    endif
  enddo
return
end
