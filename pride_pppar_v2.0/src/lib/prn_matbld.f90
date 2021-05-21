!
!! prn_matbld.f90
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
!! Contributor: Jianghui Geng, Songfeng Yang
!! 
!!
!!
subroutine prn_matbld(prn_mat)
  implicit none
  include '../header/const.h'
  character*3 prn_mat(MAXSAT)
  character*1 str1
  character*2 str2
  integer*4 i
  do i=1,MAXSAT
    if(i>0 .and. i<=MAXSAT_G)then
      if(i<=9) then
        write(str1,"(I1)")i
        prn_mat(i)='G'//'0'//str1
      else
        write(str2,"(I2)")i
        prn_mat(i)='G'//str2
      endif
    elseif(i>MAXSAT_G .and. i<=MAXSAT_G+MAXSAT_R)then
      if((i-MAXSAT_G)<=9) then
        write(str1,"(I1)")i-MAXSAT_G
        prn_mat(i)='R'//'0'//str1
      else
        write(str2,"(I2)")i-MAXSAT_G
        prn_mat(i)='R'//str2
      endif
    elseif(i>MAXSAT_G+MAXSAT_R .and. i<=MAXSAT_G+MAXSAT_R+MAXSAT_E)then
      if((i-MAXSAT_G-MAXSAT_R)<=9) then
        write(str1,"(I1)")i-MAXSAT_G-MAXSAT_R
        prn_mat(i)='E'//'0'//str1
      else
        write(str2,"(I2)")i-MAXSAT_G-MAXSAT_R
        prn_mat(i)='E'//str2
      endif
    elseif(i>MAXSAT_G+MAXSAT_R+MAXSAT_E .and. i<=MAXSAT_G+MAXSAT_R+MAXSAT_E+MAXSAT_C)then
      if((i-MAXSAT_G-MAXSAT_R-MAXSAT_E)<=9) then
        write(str1,"(I1)")i-MAXSAT_G-MAXSAT_R-MAXSAT_E
        prn_mat(i)='C'//'0'//str1
      else
        write(str2,"(I2)")i-MAXSAT_G-MAXSAT_R-MAXSAT_E
        prn_mat(i)='C'//str2
      endif
    elseif(i>MAXSAT_G+MAXSAT_R+MAXSAT_E+MAXSAT_C .and. i<=MAXSAT_G+MAXSAT_R+MAXSAT_E+MAXSAT_C+MAXSAT_J)then
      if((i-MAXSAT_G-MAXSAT_R-MAXSAT_E-MAXSAT_C)<=9) then
        write(str1,"(I1)")i-MAXSAT_G-MAXSAT_R-MAXSAT_E-MAXSAT_C
        prn_mat(i)='J'//'0'//str1
      else
        write(str2,"(I2)")i-MAXSAT_G-MAXSAT_R-MAXSAT_E-MAXSAT_C
        prn_mat(i)='J'//str2
      endif
    endif
  enddo
return
end
