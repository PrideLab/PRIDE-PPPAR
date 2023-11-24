!
!! read_att.f90
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
!!
!! purpose   : read and interpolate LEO satellite attitude for OBX files
!!
!! parameter :
!!    input  : latfil -- LEO att file name
!!             jd,sod -- time requested
!!    output : xscf,yscf,zscf -- unit vector of sc-fixed x-axes in J2000
!!
!
subroutine read_lat(latfil, jd, sod, mate2j, xscf, yscf, zscf, mat)
  implicit none
  include "../header/const.h"
  !
  !! input
  character(*) latfil
  integer*4 jd
  real*8 sod
  real*8 mate2j(3,3)
  !
  !! output
  !real*8 xscf(1:*), yscf(1:*), zscf(1:*)
  real*8 xscf(1:3), yscf(1:3), zscf(1:3)
  !
  !! local
  logical*1 lfirst, lecef
  integer*4 i, j, jj, k, lfn, ierr
  character*4 leoid
  character*3 ii
  integer*4 jdf(2), jdx, iy, imon, id, ih, im, jd0
  real*8 sec, sodf(2), sodx, dt1, dt2, dt, sod0, dwnd
  real*8 dintv, qq(2, 4), q(4), qintp(4), qs(4), qe(4)
  real*8 mat(3,3), mat_qua(3,3), det
  character*256 line
  !
  !! function used
  integer*4 get_valid_unit, modified_julday, pointer_string
  real*8 timdif

  data lfirst, lecef, leoid/.true., .true., ''/
  save lfirst, lfn, jdf, sodf, qq, dintv, jd0, sod0, leoid
  !
  !! first enter and open LEO att file
  if (lfirst) then
    lfirst = .false.
    lfn = get_valid_unit(10)
    open(lfn, file=latfil, status='old', iostat=ierr)
    if (ierr .ne. 0) then
      write(*,'(3a)') '###WARNING(read_lat): open ',trim(latfil),' error. No PCO correction on LEO end'
      lfn=0
      return
    endif
    dintv = 30.d0
    line = ' '
    do while(line(1:3) .ne. '## ')
      read(lfn, '(a)', end=200, err=100) line
      if (index(line, 'FRAME_TYPE') .ne. 0 .and. index (line, 'ECI') .ne. 0) lecef = .false.
      if (index(line, 'LEO_SATELLITE') .ne. 0) leoid = line(22:25)
    enddo
    
    backspace lfn
    ! read two epoch attitudes
    k = 0
    jdf = 0
    do while (k .le. 2)
      read (lfn, '(a)', end=200, err=100) line
      if (line(1:3) .eq. '## ') then
        read (line(4:), *, iostat=ierr) iy, imon, id, ih, im ,sec, j 
        !
        !! check time tag
        call yr2year(iy)
        jdx = modified_julday(id, imon, iy)
        sodx = ih*3600.0 + im*60.0 +sec
        k = k + 1
        if (k .gt. 2) then 
            backspace lfn
            exit
        endif
        if (jdf(k) .eq. 0) then
          jdf(k) = jdx
          sodf(k) = sodx
        endif
      else if (line(1:4) .eq. ' ATT') then
        read (line(5:), *, iostat=ierr) jj, q(1), q(2), q(3), q(4)
        !
        !! set prn and attitudes
        qq(k, 1) = q(1)
        qq(k, 2) = q(2)
        qq(k, 3) = q(3)
        qq(k, 4) = q(4)
      endif
    enddo
    jd0 = jdf(1)
    sod0 = sodf(1)
    dintv = timdif(jdf(2), sodf(2), jdf(1), sodf(1))
    write(*, '(a, f7.1)') '%%%MESSAGE(read_lat): LEO satellite attitude interval ', dintv
  endif
  do i=1,3
    xscf(i) = 0.d0
    yscf(i) = 0.d0
    zscf(i) = 0.d0
  enddo
  if(lfn.eq.0) return

  !
  !! check time tag
  dwnd=0.2d0 ! MAXWND
10 dt1 = timdif(jd, sod, jdf(1), sodf(1))
  dt2 = timdif(jd, sod, jdf(2), sodf(2))
  if (dt1 .lt. -dwnd) then
    write (*, '(a,i5,f10.3,a)') '###WARNING(read_att): t < trefatt ', jd, sod, '. No PCO correction on LEO at '
    return
  else if (dt2 .gt. dwnd) then
    !
    !! transfer attitudes
    qq(1, 1) = qq(2, 1)
    qq(1, 2) = qq(2, 2)
    qq(1, 3) = qq(2, 3)
    qq(1, 4) = qq(2, 4)
    jdf(1) = jdf(2)
    sodf(1) = sodf(2)
    jdf(2) = 0
    !
    !! read next epoch attitude
    read (lfn, '(a)', end=200, err=100) line
    do while(line(1:3) .ne. '## ')
      read(lfn, '(a)', end=200, err=100) line
    enddo
    read (line(4:), *, iostat=ierr) iy, imon, id, ih, im ,sec, j 
    !
    !! check time tag
    call yr2year(iy)
    jdx = modified_julday(id, imon, iy)
    sodx = ih*3600.d0 + im*60.d0 + sec
    if (jdf(2) .eq. 0) then
      jdf(2) =jdx
      sodf(2) = sodx
    endif
    do while (.true.) 
      read (lfn, '(a)', end=200, err=100) line
      if (line(1:4) .ne. ' ATT') exit
      read (line(5:), *, iostat=ierr) jj, q(1), q(2), q(3), q(4)
      !
      !! set prn and attitudes
      qq(2, 1) = q(1)
      qq(2, 2) = q(2)
      qq(2, 3) = q(3)
      qq(2, 4) = q(4)
    enddo
    backspace lfn
    goto 10
  else
  !
  !! check prn exist
    dt=((jd-jdf(1))*86400.d0+(sod-sodf(1)))/&
        ((jdf(2)-jdf(1))*86400.d0+(sodf(2)-sodf(1)))
    do i= 1,4
      qs(i) = qq(1,i)
      qe(i) = qq(2,i)
    enddo
    call slerp(qs,qe,dt,qintp)
    call leoqua2mat(qintp,mat_qua, leoid, lfirst)
    !call leoqua2mat(qintp,mat_qua, lecef)
    call matinv(mat_qua,3,3,det)
    mat = mat_qua
    if(lecef) then
      call matmpy(mate2j, mat_qua, mat, 3, 3, 3)
    endif
    do i=1,3
      xscf(i)=mat(i,1) ! first rank
      yscf(i)=mat(i,2)
      zscf(i)=mat(i,3)
    enddo
  endif

  return

100 write (*, '(2a)') '***ERROR(read_att) : read file ', trim(latfil)
  call exit(1)
200 write (*, '(3a)') '###WARNING(read_att) : end of file ', trim(latfil),'. Nominal attitudes used for all'
  lfn=0
  return
end
