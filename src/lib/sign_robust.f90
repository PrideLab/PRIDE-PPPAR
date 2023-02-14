!
!! get_wgt_mean.f90
!!
!!    Copyright (C) 2022 by Wuhan University
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
!! Contributor: Jianghui Geng, Jihang Lin
!! 
!!
!!
!! purpose  : 1-dimensional sign-constrained least squares robust estimation
!! parameter:
!!    input : nxl,rxl -- element array
!!    output: flg -- flag of each element
!!            ndl -- # of deleted
!!
!
subroutine sign_robust(nxl,rxl,flg,drang,ndl)
implicit none

integer*4 nxl,ndl,flg(1:*)
real*8    drang,rxl(1:*)

!
!! local
integer*4 i,j,k,nmed,nmd
real*8 rmed(2),rmd(2),resi(nxl),pnew(nxl),p(nxl),cost_new,cost

!
!! initialize
p = 1.d9
pnew = 0.d0

!
!! iterate until no elements removed
if(count(flg(1:nxl).lt.2).gt.2) then
  cost=1.d12
!
!! locate median value(s) of original elements
  call med1(nxl,rxl,flg,nmed,rmed)
  do i=1,nmed
!
!! absolute residuals
    do j=1,nxl
      resi(j)=0.d0
      if(flg(j).ge.2) cycle
      resi(j)=dabs(rxl(j)-rmed(i))
    enddo
!
!! locate median value(s) of absolute residuals
    call med1(nxl,resi,flg,nmd,rmd)
    do j=1,nmd
!
!! flag elements possibly contaminated & compute cost values
      cost_new=0.d0
      rmd(j)=rmd(j)/0.6745  !!! chen
      do k=1,nxl
        pnew(k)=0.d0
        if(flg(k).ge.2) cycle
        if(resi(k).lt.drang.or.resi(k).le.3.d0*rmd(j)) then
          pnew(k)=1.d0
          cost_new=cost_new+resi(k)**2
        endif
      enddo
      cost_new=cost_new/count(pnew(1:nxl).eq.1.d0)
!
!! save flags with smallest cost value
      if(cost_new.lt.cost) then
        cost=cost_new
        p(1:nxl)=pnew(1:nxl)
      endif
    enddo
  enddo
!
!! make final decision
  if(count(p(1:nxl).eq.0.d0).gt.count(flg(1:nxl).ge.2)) then
    do i=1,nxl
      if(flg(i).lt.2.and.p(i).eq.0.d0) then
        flg(i)=3
      endif
    enddo
  endif
endif
ndl=count(flg(1:nxl).ge.2)
return
end
!
!! median values
subroutine med1(nxl,rxl,flg,nd,rd)
implicit none
integer*4 i,j,nxl,nd,flg(1:*),ix1
real*8 rxl(1:*),rd(1:*),tmp(nxl),tp
! sort elements from small to large
tmp(1:nxl)=rxl(1:nxl)
do i=1,nxl-1
  if(flg(i).ge.2) cycle
  do j=i+1,nxl
    if(flg(j).ge.2) cycle
    if(tmp(i).gt.tmp(j)) then
      tp    =tmp(i)
      tmp(i)=tmp(j)
      tmp(j)=tp
    endif
  enddo
enddo
! locate median value(s)
j=count(flg(1:nxl).lt.2)
if(mod(j,2).eq.0) then
  nd=2
  rd(1)=tmp(ix1(j/2,nxl,flg))
  rd(2)=tmp(ix1(j/2+1,nxl,flg))
else
  nd=1
  rd(1)=tmp(ix1(j/2+1,nxl,flg))
endif
return
end
!
!! return ith good data
integer*4 function ix1(ind,nxl,flg)
implicit none

integer*4 ind,nxl,flg(1:*),i
i=0
do ix1=1,nxl
  if(flg(ix1).ge.2) cycle
  i=i+1
  if(i.eq.ind) return
enddo
write(*,'(a)') '***ERROR(sign_robust): value not found'
call exit(1)
end
