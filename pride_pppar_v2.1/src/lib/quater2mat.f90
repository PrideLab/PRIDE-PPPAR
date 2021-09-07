!
!! quater2mat.f90
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
!! Contributor: Jianghui Geng, Songfeng Yang

!!
!!
!! purpose  : calculate rotation matrix from quanternion
!!
!! parameter: q(4)  -- quanternion parameters
!!            q=(q0,q1,q2,q3)     q0: scalar part 
!!                               q1,q2,q3: vectorial part
!!            xmat(3,3) -- rotation matrix for sc-fixed to inertial system
!!
!!

subroutine quater2mat(q,xmat)
implicit none
real*8 q(4),xmat(3,3)
!
real*8 q00,q11,q22,q33,q01,q02,q03,q12,q13,q23

  q00=q(1)*q(1);q11=q(2)*q(2);q22=q(3)*q(3);q33=q(4)*q(4)
  q01=q(1)*q(2);q02=q(1)*q(3);q03=q(1)*q(4)
  q12=q(2)*q(3);q13=q(2)*q(4)
  q23=q(3)*q(4)
  xmat(1,1)=q00+q11-q22-q33
  xmat(1,2)=2.d0*(q12-q03)
  xmat(1,3)=2.d0*(q13+q02)
  xmat(2,1)=2.d0*(q12+q03)
  xmat(2,2)=q00-q11+q22-q33
  xmat(2,3)=2.d0*(q23-q01)
  xmat(3,1)=2.d0*(q13-q02)
  xmat(3,2)=2.d0*(q23+q01)
  xmat(3,3)=q00-q11-q22+q33

return
end
