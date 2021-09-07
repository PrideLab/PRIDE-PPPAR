!
!! rdionex.f90
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
!!
!! read ionex file
subroutine rdionex(flntec,nprn,prn,IM)
implicit none
include '../header/const.h'
include '../header/ionex.h'

type(ionex)   IM
integer*4     nprn
character*3   prn(1:*)
character*20  flntec
!
!! local
integer*4     i,j,k,lfn,ierr
integer*4     iy,imon,id,ih,imi,dumm(73)
real*8        tmp,xscale
character*3   iprn
character*80  str
character*256 line
character*365 strdata
!
!! function used
integer*4     modified_julday,pointer_string
integer*4     get_valid_unit
!
!! open file
lfn=get_valid_unit(10)
open(lfn,file=flntec,status='old',iostat=ierr)
if(ierr.ne.0) then
  write(*,'(2a)') '***ERROR(rdionex): open ',flntec(1:len_trim(flntec))
  call exit(1)
endif
!
!! initialization
IM%lrot=.true. ! rotate maps before interpolation
IM%ctm ='UTC'  ! time system
IM%sdcb=0.d0
!
!! read header
read(lfn,'(a)') line
read(line,*) tmp
if(tmp.ne.1.d0) then
  write(*,'(a,f3.1)') '***ERROR(rdionex): version ',tmp
  call exit(1)
endif
do while(.true.)
  read(lfn,'(a)') line
  if(line(61:73).eq.'END OF HEADER') exit
  if(line(61:68).eq.'INTERVAL') then
    read(line,*) IM%dintv
  else if(line(61:77).eq.'# OF MAPS IN FILE') then
    read(line,*) IM%bmap
    allocate(IM%mjd(1:IM%bmap))
    allocate(IM%sod(1:IM%bmap))
    allocate(IM%tecmp(1:maxlon,1:maxlat,1:IM%bmap))
    allocate(IM%rmsmp(1:maxlon,1:maxlat,1:IM%bmap))
    IM%mjd(1:IM%bmap)=0
    IM%sod(1:IM%bmap)=0.d0
    IM%tecmp(1:maxlon,1:maxlat,1:IM%bmap)=0.d0
    IM%rmsmp(1:maxlon,1:maxlat,1:IM%bmap)=0.d0
  else if(line(61:71).eq.'BASE RADIUS') then
    read(line,*) IM%bradius
  else if(line(61:73).eq.'MAP DIMENSION') then
    read(line,*) i
    if(i.ne.2) then
      write(*,'(a,i4)') '***ERROR(rdionex): map dimension ',i
      call exit(1)
    endif
  else if(line(61:76).eq.'MAPPING FUNCTION') then
    IM%cmf=line(3:6)
  else if(line(61:78).eq.'HGT1 / HGT2 / DHGT') then
    read(line,*) IM%height
  else if(line(61:78).eq.'LAT1 / LAT2 / DLAT') then
    read(line,*) IM%lat1,IM%lat2,IM%dlat
    IM%blat=nint((IM%lat2-IM%lat1)/IM%dlat)+1
    if(IM%blat.gt.maxlat) then
      write(*,'(a,i8)') '***ERROR(rdionex): # of latitudes ',IM%blat
      call exit(1)
    endif
  else if(line(61:78).eq.'LON1 / LON2 / DLON') then
    read(line,*) IM%lon1,IM%lon2,IM%dlon
    IM%blon=nint((IM%lon2-IM%lon1)/IM%dlon)+1
    if(IM%blon.gt.maxlon) then
      write(*,'(a,i8)') '***ERROR(rdionex): # of longitudes ',IM%blon
      call exit(1)
    endif
  else if(line(61:68).eq.'EXPONENT') then
    read(line,*) xscale
    xscale=10.d0**xscale
  else if(line(61:76).eq.'PRN / BIAS / RMS') then
    if(line(4:4).eq.' ') line(4:4)='G'
    if(line(5:5).eq.' ') line(5:5)='0'
    i=pointer_string(nprn,prn,line(4:6))
    if(i.ne.0) then
      read(line(7:),*) IM%sdcb(i)
      IM%sdcb(i)=IM%sdcb(i)*vlight*1.d-9 ! ns to meter
    endif
  endif
enddo
!
!! read data body
read(lfn,'(a)') line
do while(line(61:76).eq.'START OF TEC MAP')
  read(line,*) i
  if(i.gt.IM%bmap) then
    write(*,'(a)') '***ERROR(rdionex): too many TEC maps '
    call exit(1)
  endif
  read(lfn,*) iy,imon,id,ih,imi,tmp
  IM%mjd(i)=modified_julday(id,imon,iy)
  IM%sod(i)=ih*3600.d0+imi*60.d0+tmp
  j=0
  do while(.true.)
    read(lfn,'(a)') line
    if(line(61:80).ne.'LAT/LON1/LON2/DLON/H') exit
    j=j+1
    if(j.gt.IM%blat) then
      write(*,'(a)') '***ERROR(rdionex): too many latitudes '
      call exit(1)
    endif
    strdata=' '
    read(lfn,'(a80)') strdata
    read(lfn,'(a80)') str
    strdata=strdata(1:len_trim(strdata))//str
    read(lfn,'(a80)') str
    strdata=strdata(1:len_trim(strdata))//str
    read(lfn,'(a80)') str
    strdata=strdata(1:len_trim(strdata))//str
    read(lfn,'(a45)') str
    strdata=strdata(1:len_trim(strdata))//str(1:45)
    read(strdata,'(73i5)',iostat=ierr) dumm(1:IM%blon)
    if (ierr .ne. 0) then
      dumm(1:IM%blon)=9999
      write(*,'(a)') '###WARNING(rdionex): err file, set as 9999'
    endif
    do k=1,IM%blon
      IM%tecmp(k,j,i)=dumm(k)
      if(dumm(k).ne.9999) IM%tecmp(k,j,i)=dumm(k)*xscale ! to TECU
    enddo
  enddo
  read(lfn,'(a)',end=100) line
enddo
!
!! read rms
do while(line(61:76).eq.'START OF RMS MAP')
  read(line,*) i
  if(i.gt.IM%bmap) then
    write(*,'(a)') '***ERROR(rdionex): too many RMS maps '
    call exit(1)
  endif
  read(lfn,*) iy,imon,id,ih,imi,tmp
  if(IM%mjd(i).ne.modified_julday(id,imon,iy).or.&
     IM%sod(i).ne.ih*3600.d0+imi*60.d0+tmp) then
    write(*,'(a,i4)') '***ERROR(rdionex): RMS map time ',i
    call exit(1)
  endif
  j=0
  do while(.true.)
    read(lfn,'(a)') line
    if(line(61:80).ne.'LAT/LON1/LON2/DLON/H') exit
    j=j+1
    if(j.gt.IM%blat) then
      write(*,'(a)') '***ERROR(rdionex): too many latitudes '
      call exit(1)
    endif
    strdata=' '
    read(lfn,'(a80)') strdata
    read(lfn,'(a80)') str
    strdata=strdata(1:len_trim(strdata))//str
    read(lfn,'(a80)') str
    strdata=strdata(1:len_trim(strdata))//str
    read(lfn,'(a80)') str
    strdata=strdata(1:len_trim(strdata))//str
    read(lfn,'(a45)') str
    strdata=strdata(1:len_trim(strdata))//str(1:45)
    read(strdata,'(73i5)') dumm(1:IM%blon)
    do k=1,IM%blon
      IM%rmsmp(k,j,i)=dumm(k)
      if(dumm(k).ne.9999) IM%rmsmp(k,j,i)=dumm(k)*xscale ! to TECU
    enddo
  enddo
  read(lfn,'(a)',end=100) line
enddo
close(lfn)

return
100 write(*,'(a)') '***ERROR(rdionex): end of file'
call exit(1)

Entry clean_rdionex(IM)
  if(associated(IM%mjd))  deallocate(IM%mjd)
  if(associated(IM%sod))  deallocate(IM%sod)
  if(associated(IM%tecmp))deallocate(IM%tecmp)
  if(associated(IM%rmsmp))deallocate(IM%rmsmp)
end
