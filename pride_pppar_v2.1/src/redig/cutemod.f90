!
!! cutemod.f90
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
!! Contributor: Jianghui Geng
!
real*8 function cutemod(d1,d2)
implicit none
include '../header/const.h'

real*8 d1,d2,tmp

cutemod=dmod(d1,d2)
if(dabs(cutemod).lt.MAXWND) then
  cutemod=0.d0
else if(dabs(cutemod-d2).lt.MAXWND) then
  cutemod=0.d0
endif
return
end
