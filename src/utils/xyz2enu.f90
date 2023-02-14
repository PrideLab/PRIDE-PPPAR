!
!! xyz2enu.f90
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
program xyz2enu
implicit none

integer*4 ierr,k,narg,mjd,satnum,satnum_G,satnum_R,satnum_E,satnum_C,satnum_3,satnum_J
real*8 xf(3),x(3),blh(3),sod,rot(3,3),x_all(3),pdop
character*50 x_ref(3)
real*8 mean(3),sig(3),std(3)
character*1 cin
character*256 kinfil,difkin
character*256 line
integer*4 num_all,num_err
logical*1 tag_header

narg = iargc()
if (narg .ne. 2 .and. narg .ne. 5) then
  write(*,*)"Usage: xyz2enu kin_2020001 enu_2020001 [x_ref y_ref z_ref]"
  write(*,*)"if no x_ref y_ref z_ref, default x_avg y_avg z_avg"
  write(*,*)"output in file (m): sod E N U"
  write(*,*)"output in screen (m): STD"
  call exit(1)
endif
if (narg .eq. 2 .or. narg .eq. 5) then
  call getarg(1,kinfil)
  call getarg(2,difkin)
endif
if (narg .eq. 5) then
  call getarg(3,x_ref(1))
  call getarg(4,x_ref(2))
  call getarg(5,x_ref(3))
  read(x_ref(1),*) xf(1)
  read(x_ref(2),*) xf(2)
  read(x_ref(3),*) xf(3)
endif

open(15,file= difkin)
open(10,file= kinfil,status='old',iostat=ierr)
if(ierr.ne.0) call exit(1)

!
!! read header
tag_header=.false.
do while (.true.)
  read (10, '(a)', end=50) line
  if (index(line, 'END OF HEADER') .ne. 0) then
    tag_header=.true.
    exit
  endif
enddo
read(10,'(a)') line

100 continue
if (narg .eq. 2) then
  k=0
  x_all=0.d0
  do while(.true.)
    read(10,'(a)',iostat=ierr) line
    if(ierr.ne.0) exit
    read(line,'(15x,a1)') cin
    read(line(18:),*) x(1:3)
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
num_all=0
num_err=0
if (tag_header) then
  do while (.true.)
    read (10, '(a)') line
    if (index(line, 'END OF HEADER') .ne. 0) exit
  enddo
  read(10,'(a)') line
endif
do while(.true.)
  read(10,'(a)',iostat=ierr) line
  if(ierr.ne.0) exit
  read(line,'(i5,f10.2,1x,a1,3f13.3,45x,i7,6(1x,i2.2),f7.2)',iostat=ierr) &
       mjd,sod,cin,x(1:3),satnum,satnum_G,satnum_R,satnum_E,satnum_C,satnum_3,satnum_J,pdop
  num_all=num_all+1
  if(cin.eq.'*') then
    num_err=num_err+1
    cycle
  endif
  x(1:3)=x(1:3)-xf(1:3)
  call matmpy(x,rot,x,1,3,3)
  k=k+1
  mean(1:3)=mean(1:3)+x(1:3)
  sig(1:3)=sig(1:3)+x(1:3)**2
  write(15,'(i5,f10.2,3f13.3,i7,6(1x,i2.2),f7.2)') &
        mjd,sod,x(1:3),satnum,satnum_G,satnum_R,satnum_E,satnum_C,satnum_3,satnum_J,pdop
enddo

mean(1:3)=mean(1:3)/k
sig(1:3)=dsqrt(sig(1:3)/k)
std(1:3)=dsqrt(sig(1:3)**2-mean(1:3)*mean(1:3))
write(*,'(a4,1x,3f14.4)') "STD:",std(1:3) !STD:m

close(15)
close(10)
return

50 rewind(10)
goto 100

end
