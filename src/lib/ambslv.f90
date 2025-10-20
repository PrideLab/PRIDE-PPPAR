!
!! ambslv.f90
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
!! purpose  : resolve ambiguity
!! parameter:
!!    input : ncad   -- # of candidate ambiguities
!!            q22    -- cofactor matrix
!!    output: bias   -- float / fixed ambiguity estimates
!!            disall -- norm of optimum & suboptimum solutions
!
subroutine ambslv(ncad, q22, bias, disall)
  implicit none

  integer*4 ncad,m,i,j,k,info
  real*8 q22(1:*), bias(1:*), disall(1:*),add(2)
  real*8, dimension(:),allocatable::a,s
  real*8, dimension(:),allocatable::Q,zn

!
!! local
  real*8 dump

  if (ncad .gt. 1) then
    allocate(a(ncad))
    allocate(Q(ncad*ncad))
  
    m=2
    allocate(s(m))
    allocate(zn(ncad*m))
    
    k=0
    do i=1,ncad
      a(i)=bias(i)
      do j=1,ncad
        if(j.ge.i)then
            k=k+1
            Q(j+(i-1)*ncad)=q22(k)
            Q(i+(j-1)*ncad)=q22(k)
        end if
      enddo
    enddo

    call lambda(ncad,m,a,Q,zn,s,add,info)

    do i=1,ncad
      bias(i)=zn(i)
    end do

    do i=1,m
      disall(i)=s(i)
    end do

    deallocate(a)
    deallocate(Q)
    deallocate(s)
    deallocate(zn)
  else
    dump = bias(1)
    bias(1) = nint(bias(1))*1.d0
    dump = bias(1) - dump
    disall(1) = dump/q22(1)*dump
    dump = 1.d0 - dabs(dump)
    disall(2) = dump/q22(1)*dump
  endif

  return
end
