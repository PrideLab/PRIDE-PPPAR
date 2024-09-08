!
!! const.h
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
!! Contributor: Jihang Lin
!!
!!    Absolute Observation Bias Structure for PRIDE PPP-AR v3.0
!

type absbia
  character(3) :: tna    = ''           !! type, number and attribute
  integer(4)   :: period = 0            !! period of each bias record
  integer(4)   :: length = 0            !! length of bias records
  real(8)      :: docb   = 0.d0         !! day boundary discontinuity
  real(8), allocatable :: val(:)        !! bias estimates
  real(8), allocatable :: grd(:)        !! bias gradient/slop estimates 
end type
