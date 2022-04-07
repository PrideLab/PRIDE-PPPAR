!
!! read_igserp.f90
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
!! Contributor: Maorong Ge, Jianghui Geng, Songfeng Yang, Jihang Lin
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
integer*4 lfn,i,ierr
real*8 mjdx(2),dat(3,2)
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
!
!! initialize
    mjdx=0
    dat=0.d0
!
!! read two records
    read(lfn,'(a)',end=100) line
    do while(len_trim(line).eq.0)
      read(lfn,'(a)',end=100) line
    enddo
    !! the first record
    if(line(1:9).ne.'version 2'.and.line(1:9).ne.'VERSION 2') then
      write(*,'(2a)') '***ERROR(read_igserp): unknown version ',line(1:9)
      call exit(1)
    endif
    line=' '
    lpcorr=.true.
    do while(.not.isheader(line))
      read(lfn,'(a)',end=100) line
    enddo
    if(index(line,'TAI').ne.0) lpcorr=.false.
    do while(.not.isdata(line))
      read(lfn,'(a)',end=100) line
    enddo
    read(line,*,err=200) mjdx(1),(dat(i,1),i=1,3)

    read(lfn,'(a)',end=100) line
    do while(len_trim(line).eq.0)
      read(lfn,'(a)',end=100) line
    enddo
    !! the second record
    if(line(1:9).eq.'version 2'.or.line(1:9).eq.'VERSION 2') then
      do while(.not.isheader(line))
        read(lfn,'(a)',end=100) line
      enddo
      do while(.not.isdata(line))
        read(lfn,'(a)',end=100) line
      enddo
    endif
    read(line,*,err=200) mjdx(2),(dat(i,2),i=1,3)
    call zonaleffect(mjdx(1),dat(:,1),lpcorr)
    call zonaleffect(mjdx(2),dat(:,2),lpcorr)
  endif
!
!! compare time tag
 10 continue
  if(jd+sod/86400.d0.ge.mjdx(1).and.jd+sod/86400.d0.le.mjdx(2)) then
    if(mjdx(1).eq.mjdx(2)) then
      call assignvalue(dat(:,1),tmur,pole)
    else
      call interpolate(mjdx,dat,tmur,pole)
    endif
  else if(jd+sod/86400.d0.lt.mjdx(1)) then
    !! try to find a former record
    if(.not.rewindtag) then
      rewindtag=.true.
      rewind lfn
    else
      !! rewinded but no former record found 
      if(jd+sod/86400.0.lt.mjdx(1)) then
        call assignvalue(dat(:,1),tmur,pole)
        return
      endif
    endif
    read(lfn,'(a)',end=100) line
    do while(len_trim(line).eq.0)
      read(lfn,'(a)',end=100) line
    enddo
    !! the first record is found
    if(line(1:9).eq.'version 2'.or.line(1:9).eq.'VERSION 2') then
      do while(.not.isheader(line))
        read(lfn,'(a)',end=100) line
      enddo
      do while(.not.isdata(line))
        read(lfn,'(a)',end=100) line
      enddo
    endif
    read(line,*,err=200) mjdx(1),(dat(i,1),i=1,3)

    read(lfn,'(a)',end=100) line
    do while(len_trim(line).eq.0)
      read(lfn,'(a)',end=100) line
    enddo
    !! the second record is found
    if(line(1:9).eq.'version 2'.or.line(1:9).eq.'VERSION 2') then
      do while(.not.isheader(line))
        read(lfn,'(a)',end=100) line
      enddo
      do while(.not.isdata(line))
        read(lfn,'(a)',end=100) line
      enddo
    endif
    read(line,*,err=200) mjdx(2),(dat(i,2),i=1,3)
    call zonaleffect(mjdx(1),dat(:,1),lpcorr)
    call zonaleffect(mjdx(2),dat(:,2),lpcorr)
    !! compare time tag again
    goto 10
  else if(jd+sod/86400.d0.gt.mjdx(2)) then
    !! try to find the next record
    mjdx(1)=mjdx(2)
    do i=1,3
      dat(i,1)=dat(i,2)
    enddo
    read(lfn,'(a)',end=100) line
    do while(len_trim(line).eq.0)
      read(lfn,'(a)',end=100) line
    enddo
    !! the next record is found
    if(line(1:9).eq.'version 2'.or.line(1:9).eq.'VERSION 2') then
      do while(.not.isheader(line))
        read(lfn,'(a)',end=100) line
      enddo
      do while(.not.isdata(line))
        read(lfn,'(a)',end=100) line
      enddo
    endif
    read(line,*,err=200) mjdx(2),(dat(i,2),i=1,3)
    call zonaleffect(mjdx(2),dat(:,2),lpcorr)
    !! compare time tag again
    goto 10
  endif
  return

100 continue
  endtag=.true. 
  if(int(mjdx(1)).gt.45150.and.int(mjdx(2)).le.45150) then
    mjdx(2)=mjdx(1)
  elseif(int(mjdx(1)).le.45150) then
    write(*,'(a)') '***ERROR(read_igserp): end of file igserp'
    call exit(1)
  endif
  if(mjdx(1).eq.mjdx(2)) then
    call assignvalue(dat(:,1),tmur,pole)
  else
    call interpolate(mjdx,dat,tmur,pole)
  endif
  return

200 continue
  write(*,'(2a)') '***ERROR(read_igserp): read file igserp ',line(1:len_trim(line))
  call exit(1)

  Entry igserp_reset()
  close(lfn)
  lfirst=.true.
  return

contains

  logical*4 function isheader(line)
    implicit none
    character*256 line
    if (index(line,'MJD').ne.0 .or. index(line,'UT1').ne.0 .or. index(line,'UTC').ne.0) then
      if (index(line,'TAI').ne.0 .or. index(line,'LOD').ne.0) then
        isheader = .true.
        return
      endif
    endif
    if (index(line,'mjd').ne.0 .or. index(line,'ut1').ne.0 .or. index(line,'utc').ne.0) then
      if (index(line,'tai').ne.0 .or. index(line,'lod').ne.0) then
        isheader = .true.
        return
      endif
    endif
    isheader = .false.
  end function

  logical*4 function isdata(line)
    implicit none
    character*256 line
    if (index(line,'MJD').ne.0 .or. index(line,'mjd').ne.0 &
        .or. index(line,'X').ne.0 .or. index(line,'Y').ne.0) then
      isdata = .false.
      return
    endif
    if (len_trim(line).eq.0 .or. index(line,'--').ne.0 &
        .or. index(line,'"').ne.0 .or. index(line,'/').ne.0) then
      isdata = .false.
      return
    endif
    isdata = .true.
  end function

  subroutine zonaleffect(mjd0, dat0, lpcorr)
    implicit none
    real*8 mjd0, dat0(3)
    logical*1 lpcorr
    real*8 fjd, dut1, dlod, domega, lp
    real*8 taiutc
    fjd = mjd0 + GPSTDT/86400.d0
    call RG_ZONT2((fjd-51544.5d0)/36525.d0, dut1, dlod, domega)
    if (lpcorr) then
      lp = taiutc(int(mjd0))
      dat0(3) = lp + dut1 - dat0(3)*1.d-7
    else
      dat0(3) = dut1 - dat0(3)*1.d-7
    endif
  end subroutine

  subroutine assignvalue(dat0, tmur, pole)
    implicit none
    real*8 dat0(3), tmur, pole(2)
    pole(1:2) = dat0(1:2)*1.d-6
    tmur = dat0(3)
  end subroutine

  subroutine interpolate(mjdx, datx, tmur, pole)
    implicit none
    real*8 tmur, pole(2)
    real*8 mjdx(2), datx(3,2)
    real*8 alpha
    alpha=(jd+sod/86400.d0-mjdx(1))/(mjdx(2)-mjdx(1))
    do i=1,2
      pole(i)=datx(i,1)+alpha*(datx(i,2)-datx(i,1))
      pole(i)=pole(i)*1.d-6
    enddo
    tmur=datx(3,1)+alpha*(datx(3,2)-datx(3,1))
  end subroutine

end subroutine
