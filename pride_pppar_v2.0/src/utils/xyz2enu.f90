!
!! xyz2enu.f90
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
program xyz2enu
implicit none

integer*4 ierr,k,narg
real*8 xf(3),x(3),blh(3),sod,rot(3,3),x_all(3)
character*50 x_ref(3)
real*8 mean(3),sig(3),std(3)
character*1 cin
character*50 kinfil,difkin
character*256 line
integer*4 num_all,num_err

narg = iargc()
if (narg .ne. 2 .and. narg .ne. 5) then
  write(*,*)"Usage: xyz2enu kin_2020001 enu_2020001 [x_ref y_ref z_ref]"
  write(*,*)"if no x_ref y_ref z_ref, default x_avg y_avg z_avg"
  write(*,*)"output in file (cm): sod E N U"
  write(*,*)"output in screen (cm): STD"
  stop
endif
if (narg .eq. 2 .or. narg .eq. 5) then
  call getarg(1,kinfil)
  call getarg(2,difkin)
endif
if (narg .eq. 5) then
  call getarg(3,x_ref(1))
  call getarg(4,x_ref(2))
  call getarg(5,x_ref(3))
  read(x_ref(1),'(f20.4)') xf(1)
  read(x_ref(2),'(f20.4)') xf(2)
  read(x_ref(3),'(f20.4)') xf(3)
endif

open(15,file= difkin)
open(10,file= kinfil,status='old',iostat=ierr)
if(ierr.ne.0) stop

if (narg .eq. 2) then
  k=0
  x_all=0.d0
  read(10,'(a)',iostat=ierr) line
  read(10,'(a)',iostat=ierr) line
  read(10,'(a)',iostat=ierr) line
  do while(.true.)
    read(10,'(a)',iostat=ierr) line
    if(ierr.ne.0) exit
    read(line,'(5x,f9.2,a1)') sod,cin
    read(line(17:),*) x(1:3)
    if(cin.eq.'*') cycle
    x_all(1:3)=x_all(1:3)+x(1:3)
    k=k+1
  enddo
  xf(1:3)=x_all(1:3)/k !m
endif
rewind(10)

call xyzblh(xf,1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,blh)
call rot_enu2xyz(blh(1),blh(2),rot)
k=0
read(10,'(a)',iostat=ierr) line
read(10,'(a)',iostat=ierr) line
read(10,'(a)',iostat=ierr) line
num_all=0
num_err=0
do while(.true.)
  read(10,'(a)',iostat=ierr) line
  if(ierr.ne.0) exit
  read(line,'(5x,f9.2,a1)') sod,cin
  read(line(17:),*) x(1:3)
  num_all=num_all+1
  if(cin.eq.'*') then
    num_err=num_err+1
    cycle
  endif
  x(1:3)=x(1:3)-xf(1:3)
  call matmpy(x,rot,x,1,3,3)
  k=k+1
  mean(1:3)=mean(1:3)+x(1:3)*1.d2
  sig(1:3)=sig(1:3)+(x(1:3)*1.d2)**2
  write(15,'(f10.3,3f20.2)') sod,x(1:3)*1.d2 !cm
enddo

mean(1:3)=mean(1:3)/k
sig(1:3)=dsqrt(sig(1:3)/k)
std(1:3)=dsqrt(sig(1:3)**2-mean(1:3)*mean(1:3))
write(*,'(a4,1x,3f9.2)') "STD:",std(1:3) !STD:cm

close(15)
close(10)
stop
end
