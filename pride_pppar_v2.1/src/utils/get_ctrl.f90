!
!! get_ctrl.f90
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
!! Contributor: Maorong Ge
!! 
!!
!
      program get_ctrl
! purpose  : get an item from PANDA control file by calling `findkey`.
!            This is mainly for sh_panda
!
!
      character*256 ctrlfile,keyname,keyvalue,bracket
      integer*4 iargc,narg
      logical*4 lexist
!
! function called
      character*256 findkey

      keyvalue='EMPTY'
      narg=iargc()
      if(narg.lt.2) then
        write(*,*) '***ERROR(get_ctrl): wrong args'
        write(*,*) 'get_ctrl ctrl_file control_name [bracket_name] '
        goto 100
      endif
!
! first one is control file
      call getarg(1,ctrlfile)
      inquire(file=ctrlfile,exist=lexist)
      if(.not.lexist) then
        write(*,*) '***ERROR(get_ctrl): control_file not exist',trim(ctrlfile)
        goto 100
      endif
!
! second keyname
      call getarg(2,keyname)
!
! third bracket, if there is
      bracket=' '
      if(narg.eq.3) call getarg(3,bracket)
!
! get it
      open(10,file=ctrlfile,status='old')

      keyvalue=findkey(10,keyname,bracket)
      close(10)
100   continue
      write(*,'(a)') trim(keyvalue)
      end
