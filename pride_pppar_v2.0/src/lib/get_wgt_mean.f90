!
!! get_wgt_mean.f90
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
!! Contributor: Maorong Ge, Jianghui Geng, Songfeng Yang
!! 
!!
!!
!! purpose  : compute weighted mean value of an array
!!            & sign-constrained least squares
!! parameter:
!!    input : l_edit -- edit data or not
!!            x,flg  -- data array & flag array
!!            wgt    -- weight array
!!            n      -- # of data
!!    output: ndel -- # of deleted data
!!            mean -- mean value
!!            rms  -- unweighted rms of residuals
!!            sig  -- sigma of mean value
!!
!
subroutine get_wgt_mean(l_edit,x,flg,wgt,n,mrms,ndel,mean,rms,sig)
implicit none

logical*1 l_edit
real*8 x(1:*),sig,mrms,rms,mean,hp,wgt(1:*)
integer*4 flg(1:*)
integer*4 n,ndel

integer*4 i,m,ndel_new
real*8 wgt_sum
!
!! sign-constrained least squares to remove large biases
!..call sign_robust(n,x,flg,ndel)
ndel=0
ndel_new=1
do while(ndel_new.ne.0)
!
!! mean
  m=0; mean=0.d0; wgt_sum=0.d0
  do i=1,n
    if(flg(i).ge.2) cycle
    mean=mean+wgt(i)*x(i); m=m+1; wgt_sum=wgt_sum+wgt(i)
  enddo
  if(m.ne.0) then
    mean=mean/wgt_sum
  else
    rms=999.d0; sig=999.d0; ndel=n; return
  endif
!
!! sigma & standard deviation
  sig=0.d0; rms=0.d0
  do i=1,n
    if(flg(i).ge.2) cycle
    rms=rms+(x(i)-mean)**2*wgt(i)
  enddo
 
  if(m.eq.1) then
    ndel=n-1; rms=999.d0; sig=999.d0; return
  else
    rms=dsqrt(rms/(m-1))
    sig=rms/dsqrt(wgt_sum)
  endif
!
!! if achieve goal
  if(.not.l_edit.or.3.d0*rms.lt.mrms) then
    ndel=count(flg(1:n).ge.2); return
  endif
!
!! edit data according to residuals 
  ndel_new=0; ndel=0
  do i=1,n
    if(flg(i).ge.2) then
      ndel=ndel+1
    else if(dabs(x(i)-mean).ge.3.d0*rms/wgt(i)) then
      ndel=ndel+1
      ndel_new=ndel_new+1
      flg(i)=2
    endif
  enddo
enddo

return
end
