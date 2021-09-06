!
!! read_igserp.f90
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
!! Contributor: Maorong Ge, Jianghui Geng, Songfeng Yang
!! 
!!
!!
!! purpose  : read IGS ERP & return TAI-UT1R
!! parameter:
!!    input : 
!!            erpfil -- ERP file
!!            jd,sod -- requested time
!!    output: tmur   -- TAI-UT1R
!!            pole   -- pole movement
!!
!
subroutine read_igserp(erpfil,jd,sod,tmur,pole)
implicit none
include '../header/const.h'

integer*4 jd
real*8 sod,tmur,pole(2)
character*(*) erpfil
!
!! local
logical*1 lfirst,lpcorr,endtag,rewindtag
integer*4 lfn,i,jd_tdt,ierr
real*8 lp,dut1,dlod,domega,sod_tdt,alpha,fjd,mjdx(2),dat(3,2)
character*256 line
!
!! function called
integer*4 get_valid_unit
real*8    taiutc

data lfirst,endtag,rewindtag/.true.,.false.,.false./
save lfirst,lpcorr,lfn,mjdx,dat,endtag,rewindtag

if(endtag) then
  goto 100
endif
if(lfirst) then
  lfirst=.false.
  lfn=get_valid_unit(10)
  open(lfn,file=erpfil,status='old',iostat=ierr)
  if(ierr.eq.0) then
    write(*,'(2a)') '%%%MESSAGE(read_igserp): ERP read ',&
                     erpfil(1:len_trim(erpfil))
  else
    write(*,'(2a)') '***ERROR(read_igserp): open file ',&
                     erpfil(1:len_trim(erpfil))
    call exit(1)
  endif
  read(lfn,'(a)',end=100) line
  if(line(1:9).ne.'version 2'.and.line(1:9).ne.'VERSION 2') then
    write(*,'(2a)') '***ERROR(read_igserp): unknown version ',line(1:9)
    call exit(1)
  endif
  line=' '
  lpcorr=.true.
  do while(index(line,'MJD').eq.0.or.index(line,'UT1').eq.0.or.&
           index(line,'UTC').eq.0.and.index(line,'TAI').eq.0.or.index(line,'LOD').eq.0)
    read(lfn,'(a)',end=100) line
  enddo
  if(index(line,'TAI').ne.0) lpcorr=.false.
  read(lfn,'(a)',end=100) line
!
!! read two records
  mjdx=0
  dat=0.d0
  read(lfn,'(a)',end=100) line
  read(line,*,err=200) mjdx(1),(dat(i,1),i=1,3)
  read(lfn,'(a)',end=100) line
  if(line(1:9).eq.'version 2'.or.line(1:9).eq.'VERSION 2') then
    do while(index(line,'MJD').eq.0.or.index(line,'UT1').eq.0.or.&
             index(line,'UTC').eq.0.and.index(line,'TAI').eq.0.or.index(line,'LOD').eq.0)
      read(lfn,'(a)',end=100) line
    enddo
    read(lfn,'(a)',end=100) line
    read(lfn,'(a)',end=100) line
  endif
  read(line,*,err=200) mjdx(2),(dat(i,2),i=1,3)
  do i=1,2
    fjd=mjdx(i)+GPSTDT/86400.d0
    call RG_ZONT2((fjd-51544.5d0)/36525.d0,dut1,dlod,domega)
    if(lpcorr) then
      lp=taiutc(int(mjdx(i)))
      dat(3,i)=lp+dut1-dat(3,i)*1.d-7
    else
      dat(3,i)=dut1-dat(3,i)*1.d-7
    endif
  enddo
endif
!
!! compare time tag
10 continue
if(jd+sod/86400.d0.ge.mjdx(1).and.jd+sod/86400.d0.le.mjdx(2)) then
  if(mjdx(1).eq.mjdx(2)) then
    do i=1,2
      pole(i)=dat(i,1)
      pole(i)=pole(i)*1.d-6
    enddo
    tmur=dat(3,1)
  else
    alpha=(jd+sod/86400.d0-mjdx(1))/(mjdx(2)-mjdx(1))
    do i=1,2
      pole(i)=dat(i,1)+alpha*(dat(i,2)-dat(i,1))
      pole(i)=pole(i)*1.d-6
    enddo
    tmur=dat(3,1)+alpha*(dat(3,2)-dat(3,1))
  endif
else if(jd+sod/86400.d0.lt.mjdx(1)) then
  if(rewindtag) then
    rewindtag=.true.
    rewind lfn
  endif
  read(lfn,'(a)',end=100) line
  if(line(1:9).eq.'version 2'.or.line(1:9).eq.'VERSION 2') then
    do while(index(line,'MJD').eq.0.or.index(line,'UT1').eq.0.or.&
             index(line,'UTC').eq.0.and.index(line,'TAI').eq.0.or.index(line,'LOD').eq.0)
      read(lfn,'(a)',end=100) line
    enddo
    read(lfn,'(a)',end=100) line
    read(lfn,'(a)',end=100) line
  endif
  read(line,*,err=200) mjdx(1),(dat(i,1),i=1,3)
  read(lfn,'(a)',end=100) line
  if(line(1:9).eq.'version 2'.or.line(1:9).eq.'VERSION 2') then
    do while(index(line,'MJD').eq.0.or.index(line,'UT1').eq.0.or.&
             index(line,'UTC').eq.0.and.index(line,'TAI').eq.0.or.index(line,'LOD').eq.0)
      read(lfn,'(a)',end=100) line
    enddo
    read(lfn,'(a)',end=100) line
    read(lfn,'(a)',end=100) line
  endif
  read(line,*,err=200) mjdx(2),(dat(i,2),i=1,3)
  do i=1,2
    fjd=mjdx(i)+GPSTDT/86400.d0
    call RG_ZONT2((fjd-51544.5d0)/36525.d0,dut1,dlod,domega)
    if(lpcorr) then
      lp=taiutc(int(mjdx(i)))
      dat(3,i)=lp+dut1-dat(3,i)*1.d-7
    else
      dat(3,i)=dut1-dat(3,i)*1.d-7
    endif
  enddo
  goto 10
else if(jd+sod/86400.d0.gt.mjdx(2)) then
  mjdx(1)=mjdx(2)
  do i=1,3
    dat(i,1)=dat(i,2)
  enddo
  read(lfn,'(a)',end=100) line
  if(line(1:9).eq.'version 2'.or.line(1:9).eq.'VERSION 2') then
    do while(index(line,'MJD').eq.0.or.index(line,'UT1').eq.0.or.&
             index(line,'UTC').eq.0.and.index(line,'TAI').eq.0.or.index(line,'LOD').eq.0)
      read(lfn,'(a)',end=100) line
    enddo
    read(lfn,'(a)',end=100) line
    read(lfn,'(a)',end=100) line
  endif
  read(line,*,err=200) mjdx(2),(dat(i,2),i=1,3)
  fjd=mjdx(2)+GPSTDT/86400.d0
  call RG_ZONT2((fjd-51544.5d0)/36525.d0,dut1,dlod,domega)
  if(lpcorr) then
    lp=taiutc(int(mjdx(2)))
    dat(3,2)=lp+dut1-dat(3,2)*1.d-7
  else
    dat(3,2)=dut1-dat(3,2)*1.d-7
  endif
  goto 10
endif

return
100 endtag=.true. 
if(int(mjdx(1)).gt.45150.and.int(mjdx(2)).le.45150) then
  !write(*,'(a)') '###WARNING(read_igserp): end of file igserp'
  mjdx(2)=mjdx(1)
elseif(int(mjdx(1)).le.45150) then
  write(*,'(a)') '***ERROR(read_igserp): end of file igserp'
  call exit(1)
endif
if(mjdx(1).eq.mjdx(2)) then
  do i=1,2
    pole(i)=dat(i,1)
    pole(i)=pole(i)*1.d-6
  enddo
  tmur=dat(3,1)
else
  alpha=(jd+sod/86400.d0-mjdx(1))/(mjdx(2)-mjdx(1))
  do i=1,2
    pole(i)=dat(i,1)+alpha*(dat(i,2)-dat(i,1))
    pole(i)=pole(i)*1.d-6
  enddo
  tmur=dat(3,1)+alpha*(dat(3,2)-dat(3,1))
endif
return
200 write(*,'(2a)') '***ERROR(read_igserp): read file igserp ',line(1:len_trim(line))
call exit(1)
!
!! reset
Entry igserp_reset()
close(lfn)
lfirst=.true.
return
end
