!
!! read_docb.f90
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
!! Contributor: Jihang Lin
!! 
!!
!!
!! purpose  : read DOCB records in BIAS-SINEX format 
!! parameter:
!!    input : DOCB -- discontinuities of orbit, clock, and bias
!
subroutine read_docb(flnosb, nprn, prn, bias, docbjd)
  implicit none
  include '../header/absbia.h'
  include '../header/const.h'

! common
  integer*4     idxfrq(MAXSYS, 2)
  common        idxfrq
! parameter
  character(*)  flnosb
  integer*4     nprn
  character*3   prn(1:*)
  type(absbia)  bias(MAXSAT, MAXTYP)
  integer*4     docbjd
! local
  integer*4     i0, iprn, ntyp, ityp, jtyp, xtyp, imes
  integer*4     lfn, jd, nlen, dsod
  integer*4     iy, idoy, imon, id, ih, im, isod, iepo, jepo, ierr
  real*8        docbval 
  character*2   styp(0:3)
  character*3   cprn, ctyp
  character*256 line
! function used
  integer*4     get_valid_unit
  integer*4     modified_julday
  integer*4     pointer_string
!
!! initialize
  do iprn = 1, MAXSAT
    do ityp = 1, MAXTYP
      bias(iprn, ityp)%tna = ''
      bias(iprn, ityp)%period = 0.d0
      bias(iprn, ityp)%length = 0
      bias(iprn, ityp)%docb   = 0.d0
    end do
  end do
!
!! open file
  lfn = get_valid_unit(500)
  open (lfn, file=flnosb, status='old', iostat=ierr)
  if (ierr .ne. 0) then
    write (*, '(2a)') '###WARNING(read_docb): open file ', trim(flnosb)
    return
  end if
!
!! read header
  read (lfn, '(a)') line
  do while (.true.)
    if (index(line, '+SOLUTION/DAY_BOUNDARY_DISCONTINUITY') .ne. 0) then
      do while (.true.)
 50     read (lfn, '(a)', iostat=ierr) line
        if (ierr .ne. 0) then
          backspace lfn
          exit
        end if
        if (index(line, 'DOCB') .eq. 0) cycle
        if (line(16:19) .ne. '') cycle
        !
        !! read bias
        read (line, '(2(11x,a3),22x,(i4,1x,i3,1x,i5,1x))', err=200) &
            cprn, ctyp, iy, idoy, isod
        if (cprn .eq. '') cycle
        i0 = index(GNSS_PRIO, cprn(1:1))
        iprn = pointer_string(nprn, prn, cprn)
        if (iprn .eq.  0) cycle
        call yeardoy2monthday(iy, idoy, imon, id)
        jd = modified_julday(id, imon, iy)
        if (jd .ne. docbjd) cycle
        read (line, '(70x,f21.15)') docbval
        if (cprn(1:1) .eq. 'R' .and. ctyp(1:1) .eq. 'L') cycle
        !
        !! get index of bias attribute
        ityp = index(OBS_PRIO_SYS(i0), ctyp(3:3))
        if (ityp .eq. 0) cycle
        !
        !! get index of measurement type
        do imes = 1, 2
          jtyp = idxfrq(i0, imes)
          if (jtyp .le. 0 .or. jtyp .gt. MAXFRQ) goto 250
          write (styp(imes-1), '(a1,i1)') 'L', jtyp
          write (styp(imes+1), '(a1,i1)') 'C', jtyp
        end do
        do imes = 0, 3
          if (index(ctyp, styp(imes)) .ne. 0) exit
          if (imes .eq. 3) goto 50
        end do
        !
        !! assign value
        bias(iprn, ityp+9*imes)%length = 1
        bias(iprn, ityp+9*imes)%docb   = docbval*VLIGHT*1.d-9
      end do
    end if
    read (lfn, '(a)', end=100) line
  end do

100 continue
  if (index(line, '%=ENDBIA') .eq. 0) then 
    write (*, '(2a)') '***ERROR(read_bias): end of file ', trim(flnosb)
    call exit(1)
  else
    close (lfn)
  end if
!
!! fill in all types of phase biase
  do iprn = 1, nprn
    cprn = prn(iprn)
    if (cprn .eq. ' ') cycle
    i0 = index(GNSS_PRIO, cprn(1:1))
    ntyp = len_trim(OBS_PRIO_SYS(i0))
    !
    !! complete vacant X code biases
    do imes = 2, 3
      xtyp = index(OBS_PRIO_SYS(i0), 'X') 
      if (xtyp .eq. 0) cycle
      ityp = 0
      jtyp = 0
      if (cprn(1:1) .eq. 'G' .or. cprn(1:1) .eq. 'J') then
        select case (idxfrq(i0, imes-1))
          case (1)
            ityp =  index(OBS_PRIO_SYS(i0), 'S')
            jtyp =  index(OBS_PRIO_SYS(i0), 'L')
          case (2)
            ityp =  index(OBS_PRIO_SYS(i0), 'S')
            jtyp =  index(OBS_PRIO_SYS(i0), 'L')
          case (5)
            ityp =  index(OBS_PRIO_SYS(i0), 'I')
            jtyp =  index(OBS_PRIO_SYS(i0), 'Q')
        end select
      end if
      if (cprn(1:1) .eq. 'E') then
        select case (idxfrq(i0, imes-1))
          case (6)
            ityp =  index(OBS_PRIO_SYS(i0), 'C')
            jtyp =  index(OBS_PRIO_SYS(i0), 'C')
        end select 
      end if
      if (ityp .eq. 0 .or. jtyp .eq. 0) cycle
      xtyp = xtyp + 9*imes
      ityp = ityp + 9*imes
      jtyp = jtyp + 9*imes
      if (bias(iprn, ityp)%length .le. 0) cycle
      if (bias(iprn, jtyp)%length .le. 0) cycle
      if (bias(iprn, ityp)%length .ne. bias(iprn, jtyp)%length) cycle
      if (bias(iprn, xtyp)%length .le. 0) then
        bias(iprn, xtyp)%tna = bias(iprn, ityp)%tna(1:2)//'X'
        bias(iprn, xtyp)%period = bias(iprn, ityp)%period
        bias(iprn, xtyp)%length = bias(iprn, ityp)%length
        bias(iprn, xtyp)%docb = 1.d9
      endif
      if (bias(iprn, ityp)%length .ne. bias(iprn, xtyp)%length) cycle
      do iepo = 1, bias(iprn, xtyp)%length
        if (abs(bias(iprn, xtyp)%docb - 1.d9) .gt. 1.d-3) cycle
        bias(iprn, xtyp)%docb = (bias(iprn, ityp)%docb + &
                                 bias(iprn, jtyp)%docb)/2.d0
      end do
    end do
    !
    !! complete vacant C/Q code biases
    do imes = 2, 3
      if (cprn(1:1) .ne. 'E') cycle
      xtyp = index(OBS_PRIO_SYS(i0), 'X') 
      ityp = index(OBS_PRIO_SYS(i0), 'C') 
      jtyp = index(OBS_PRIO_SYS(i0), 'Q') 
      if (xtyp * ityp * jtyp .eq. 0) cycle
      xtyp = xtyp + 9*imes
      ityp = ityp + 9*imes
      jtyp = jtyp + 9*imes
      if (bias(iprn, xtyp)%length .le. 0) cycle
      nlen = bias(iprn, xtyp)%length
      if (bias(iprn, ityp)%length .le. 0) then
        bias(iprn, ityp)%tna = bias(iprn, xtyp)%tna(1:2)//'C'
        bias(iprn, ityp)%period = bias(iprn, xtyp)%period
        bias(iprn, ityp)%length = nlen
        bias(iprn, ityp)%docb = 1.d9
      end if
      if (abs(bias(iprn, ityp)%docb - 1.d9) .gt. 1.d-3) cycle
      bias(iprn, ityp)%docb = bias(iprn, xtyp)%docb
      if (bias(iprn, jtyp)%length .le. 0) then
        bias(iprn, jtyp)%tna = bias(iprn, xtyp)%tna(1:2)//'Q'
        bias(iprn, jtyp)%period = bias(iprn, xtyp)%period
        bias(iprn, jtyp)%length = nlen
        bias(iprn, jtyp)%docb = 1.d9
      end if
      if (abs(bias(iprn, jtyp)%docb - 1.d9) .gt. 1.d-3) cycle
      bias(iprn, jtyp)%docb = bias(iprn, xtyp)%docb
    end do
    !
    !! complete vacant phase biases 
    do imes = 0, 1
      !
      !! assumes phase biases on the same frequency are identical then find substitute
      do ityp = 9*imes+1, 9*imes+ntyp
        jtyp = ityp 
        if (abs(bias(iprn, jtyp)%docb - 1.d9) .gt. 1.d-3) exit
        if (ityp .eq. 9*imes+ntyp) jtyp = 0
      end do
      if (jtyp .eq. 0) cycle
      do xtyp = 1, ntyp
        ityp = 9*imes+xtyp
        if (bias(iprn, ityp)%length .le. 0) then
          bias(iprn, ityp)%tna = bias(iprn, jtyp)%tna(1:2)//OBS_PRIO_SYS(i0)(xtyp:xtyp)
          bias(iprn, ityp)%period = bias(iprn, jtyp)%period
          bias(iprn, ityp)%length = bias(iprn, jtyp)%length
          bias(iprn, ityp)%docb   = 1.d9
        end if
        if (abs(bias(iprn, ityp)%docb - 1.d9) .gt. 1.d-3) cycle
        bias(iprn, ityp)%docb = bias(iprn, jtyp)%docb
      end do
    end do
  end do
  return

200 continue
  write (*, '(2a)') '***ERROR(read_docb): read file ', trim(flnosb)
  call exit(1)
250 continue
  write (*, '(a,5(1x,a,2i1))') '***ERROR(read_docb): invalid frequency number:', &
    'G', idxfrq(index(GNSS_PRIO, 'G'), 1:2), &
    'R', idxfrq(index(GNSS_PRIO, 'R'), 1:2), &
    'E', idxfrq(index(GNSS_PRIO, 'E'), 1:2), &
    'C', idxfrq(index(GNSS_PRIO, 'C'), 1:2), &
    'J', idxfrq(index(GNSS_PRIO, 'J'), 1:2)
  call exit(1)
end subroutine
