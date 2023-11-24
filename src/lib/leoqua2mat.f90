!
!! quater2mat.f90
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
!! Contributor: Jianghui Geng, Songfeng Yang, Jing Zeng

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

subroutine leoqua2mat(q, xmat, leoid)
implicit none
real*8 q(4),xmat(3,3)
character*4 leoid
logical*1 lfirst
!
integer*4 i
real*8 q00,q11,q22,q33,q01,q02,q03,q12,q13,q23
  
  data lfirst /.true./
  save lfirst
  xmat = 0.d0
  q00=q(1)*q(1);q11=q(2)*q(2);q22=q(3)*q(3);q33=q(4)*q(4)
  q01=q(1)*q(2);q02=q(1)*q(3);q03=q(1)*q(4)
  q12=q(2)*q(3);q13=q(2)*q(4)
  q23=q(3)*q(4)
  if (index(leoid, 'gra') .ne. 0) then       ! GRACE & GRACE-FO
    xmat(1,1)=q00+q11-q22-q33
    xmat(1,2)=2.d0*(q12+q03)
    xmat(1,3)=2.d0*(q13-q02)
    xmat(2,1)=2.d0*(q12-q03)
    xmat(2,2)=q00-q11+q22-q33
    xmat(2,3)=2.d0*(q23+q01)
    xmat(3,1)=2.d0*(q13+q02)
    xmat(3,2)=2.d0*(q23-q01)
    xmat(3,3)=q00-q11-q22+q33
  else
    xmat(1,1)=1.d0-2.d0*q11-2*q22
    xmat(2,1)=2.d0*(q01+q23)
    xmat(3,1)=2.d0*(q02-q13)
    xmat(1,2)=2.d0*(q01-q23)
    xmat(2,2)=1.d0-2.d0*q00-2.d0*q22
    xmat(3,2)=2.d0*(q12+q03)
    xmat(1,3)=2.d0*(q02+q13)
    xmat(2,3)=2.d0*(q03-q12)
    xmat(3,3)=1.d0-2.d0*q00-2.d0*q11
  endif
  if (lfirst) lfirst=.false.
return
end
