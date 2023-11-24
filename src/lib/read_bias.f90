!
!! read_bias.f90
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
!! Contributor: Jianghui Geng, Songfeng Yang, Xingyu Chen, Jihang Lin
!! 
!!
!!
!! purpose  : read absolute code/phase biases for ambiguity resolution in PPP
!! parameter:
!!    input : OSB -- absolute code/phase biases
!
subroutine read_bias(flnosb, nprn, prn, bias, osbjd0, osbjd1)
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
  real*8        osbjd0, osbjd1
! local
  integer*4     i0, iprn, ntyp, ityp, jtyp, xtyp, imes
  integer*4     lfn, jd0, jd1, nlen, dsod
  integer*4     iy, imon, id, ih, im, iepo, jepo, ierr
  integer*4     iyear0, idoy0, isod0, imon0, id0, fsod0
  integer*4     iyear1, idoy1, isod1, imon1, id1, fsod1
  real*8        osbval, osbgrd
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
      bias(iprn, ityp)%period = 0
      bias(iprn, ityp)%length = 0
      if (allocated(bias(iprn, ityp)%val)) then
        deallocate (bias(iprn, ityp)%val)
        deallocate (bias(iprn, ityp)%grd)
      end if
    end do
  end do
!
!! open file
  lfn = get_valid_unit(500)
  open (lfn, file=flnosb, status='old', iostat=ierr)
  if (ierr .ne. 0) then
    write (*, '(2a)') '###WARNING(read_bias): open file ', trim(flnosb)
    return
  end if
!
!! read header
  read (lfn, '(a)') line
  do while (.true.)
    if (index(line, '+BIAS/SOLUTION') .ne. 0) then
      do while (.true.)
 50     read (lfn, '(a)', iostat=ierr) line
        if (ierr .ne. 0) then
          backspace lfn
          exit
        end if
        if (index(line, 'OSB') .eq. 0) cycle
        if (line(16:19) .ne. '') cycle
        !
        !! read bias
        read (line, '(2(11x,a3),7x,2(i4,1x,i3,1x,i5,1x))', err=200) &
            cprn, ctyp, iyear0, idoy0, isod0, iyear1, idoy1, isod1
        if (cprn .eq. '') cycle
        i0 = index(GNSS_PRIO, cprn(1:1))
        iprn = pointer_string(nprn, prn, cprn)
        if (iprn .eq.  0) cycle
        call yeardoy2monthday(iyear0, idoy0, imon0, id0)
        jd0 = modified_julday(id0, imon0, iyear0)
        call yeardoy2monthday(iyear1, idoy1, imon1, id1)
        jd1 = modified_julday(id1, imon1, iyear1)
        if (jd1 + isod1/864.d2 .le. osbjd0 .or. jd0 + isod0/864.d2 .ge. osbjd1) cycle
        read (line, '(70x,f21.15)') osbval
        if (len_trim(line) .gt. 124) then
          read (line, '(104x,f21.15)') osbgrd
        end if
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
        !! allocate memory
        nlen = bias(iprn, ityp+9*imes)%length
        if ((nlen .le. 0) .or. (.not. allocated(bias(iprn, ityp+9*imes)%val))) then
          bias(iprn, ityp+9*imes)%tna = ctyp(1:3)
          dsod = isod1 - isod0 + (jd1 - jd0) * 86400
          if (mod(dsod, 10) .eq. 9) dsod = dsod + 1 !! rounded to nearest whole tens
          if (dsod .gt. 86400) dsod = 86400         !! period > 1 day
          if (mod(86400, dsod) .ne. 0) dsod = 30    !! irregular periods are divided to 30 s
          nlen = 86400/dsod
          bias(iprn, ityp+9*imes)%period = dsod
          bias(iprn, ityp+9*imes)%length = nlen
          allocate (bias(iprn, ityp+9*imes)%val(1:nlen))
          allocate (bias(iprn, ityp+9*imes)%grd(1:nlen))
          bias(iprn, ityp+9*imes)%val = 1.d9        !! set NaN as 1.d9
          bias(iprn, ityp+9*imes)%grd = 1.d9        !! set NaN as 1.d9
        end if
        !
        !! assign value
        fsod0 = isod0 * 1.d0
        iepo = int(fsod0/bias(iprn, ityp+9*imes)%period)+1
        if (iepo .gt. nlen) iepo = nlen 
        if (iepo .le. 0) iepo = 1
        fsod1 = isod1 + (jd1-jd0)*864.d2
        jepo = int(fsod1/bias(iprn, ityp+9*imes)%period)+1
        if (jepo .gt. nlen) jepo = nlen 
        if (jepo .le. 0) jepo = 1
        bias(iprn, ityp+9*imes)%val(iepo:jepo) = osbval*VLIGHT*1.d-9
        if (nlen .gt. 1) bias(iprn, ityp+9*imes)%grd(iepo:jepo) = osbgrd*VLIGHT*1.d-9
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
      if (.not. allocated(bias(iprn, ityp)%val)) cycle
      if (.not. allocated(bias(iprn, jtyp)%val)) cycle
      if (bias(iprn, ityp)%length .ne. bias(iprn, jtyp)%length) cycle
      if (.not. allocated(bias(iprn, xtyp)%val)) then
        bias(iprn, xtyp)%tna = bias(iprn, ityp)%tna(1:2)//'X'
        bias(iprn, xtyp)%period = bias(iprn, ityp)%period
        bias(iprn, xtyp)%length = bias(iprn, ityp)%length
        nlen = bias(iprn, ityp)%length
        allocate (bias(iprn, xtyp)%val(1:nlen))
        allocate (bias(iprn, xtyp)%grd(1:nlen))
        bias(iprn, xtyp)%val = 1.d9
        bias(iprn, xtyp)%grd = 1.d9
      endif
      if (bias(iprn, ityp)%length .ne. bias(iprn, xtyp)%length) cycle
      do iepo = 1, bias(iprn, xtyp)%length
        if (abs(bias(iprn, xtyp)%val(iepo) - 1.d9) .gt. 1.d-3) cycle
        bias(iprn, xtyp)%val(iepo) = (bias(iprn, ityp)%val(iepo) + &
                                      bias(iprn, jtyp)%val(iepo))/2.d0
        bias(iprn, xtyp)%grd(iepo) = (bias(iprn, ityp)%grd(iepo) + &
                                      bias(iprn, jtyp)%grd(iepo))/2.d0
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
      if (.not. allocated(bias(iprn, xtyp)%val)) cycle
      nlen = bias(iprn, xtyp)%length
      if (.not. allocated(bias(iprn, ityp)%val)) then
        bias(iprn, ityp)%tna = bias(iprn, xtyp)%tna(1:2)//'C'
        bias(iprn, ityp)%period = bias(iprn, xtyp)%period
        bias(iprn, ityp)%length = nlen
        allocate (bias(iprn, ityp)%val(1:nlen))
        allocate (bias(iprn, ityp)%grd(1:nlen))
        bias(iprn, ityp)%val = 1.d9
        bias(iprn, ityp)%grd = 1.d9
      end if
      if (nlen .eq. bias(iprn, ityp)%length) then
        do iepo = 1, nlen
          if (abs(bias(iprn, ityp)%val(iepo) - 1.d9) .gt. 1.d-3) cycle
          bias(iprn, ityp)%val(iepo) = bias(iprn, xtyp)%val(iepo)
          bias(iprn, ityp)%grd(iepo) = bias(iprn, xtyp)%grd(iepo)
        end do
      end if
      if (.not. allocated(bias(iprn, jtyp)%val)) then
        bias(iprn, jtyp)%tna = bias(iprn, xtyp)%tna(1:2)//'Q'
        bias(iprn, jtyp)%period = bias(iprn, xtyp)%period
        bias(iprn, jtyp)%length = nlen
        allocate (bias(iprn, jtyp)%val(1:nlen))
        allocate (bias(iprn, jtyp)%grd(1:nlen))
        bias(iprn, jtyp)%val = 1.d9
        bias(iprn, jtyp)%grd = 1.d9
      end if
      if (nlen .eq. bias(iprn, jtyp)%length) then
        do iepo = 1, nlen
          if (abs(bias(iprn, jtyp)%val(iepo) - 1.d9) .gt. 1.d-3) cycle
          bias(iprn, jtyp)%val(iepo) = bias(iprn, xtyp)%val(iepo)
          bias(iprn, jtyp)%grd(iepo) = bias(iprn, xtyp)%grd(iepo)
        end do
      end if
    end do
    !
    !! complete vacant phase biases 
    do imes = 0, 1
      !
      !! assumes phase biases on the same frequency are identical then find substitute
      do ityp = 9*imes+1, 9*imes+ntyp
        jtyp = ityp 
        if (allocated(bias(iprn, jtyp)%val)) then
          if (any(abs(bias(iprn, jtyp)%val - 1.d9) .gt. 1.d-3)) exit
        end if
        if (ityp .eq. 9*imes+ntyp) jtyp = 0
      end do
      if (jtyp .eq. 0) cycle
      do xtyp = 1, ntyp
        ityp = 9*imes+xtyp
        if (.not. allocated(bias(iprn, ityp)%val)) then
          bias(iprn, ityp)%tna = bias(iprn, jtyp)%tna(1:2)//OBS_PRIO_SYS(i0)(xtyp:xtyp)
          bias(iprn, ityp)%period = bias(iprn, jtyp)%period
          bias(iprn, ityp)%length = bias(iprn, jtyp)%length
          nlen = bias(iprn, jtyp)%length
          allocate (bias(iprn, ityp)%val(1:nlen))
          allocate (bias(iprn, ityp)%grd(1:nlen))
          bias(iprn, ityp)%val = 1.d9
          bias(iprn, ityp)%grd = 1.d9
        end if
        do iepo = 1, bias(iprn, ityp)%length
          if (abs(bias(iprn, ityp)%val(iepo) - 1.d9) .gt. 1.d-3) cycle
          bias(iprn, ityp)%val(iepo) = bias(iprn, jtyp)%val(iepo)
          bias(iprn, ityp)%grd(iepo) = bias(iprn, jtyp)%grd(iepo)
        end do
      end do
    end do
  end do
  return

200 continue
  write (*, '(2a)') '***ERROR(read_bias): read file ', trim(flnosb)
  call exit(1)
250 continue
  write (*, '(a,5(1x,a,2i1))') '***ERROR(read_bias): invalid frequency number:', &
    'G', idxfrq(index(GNSS_PRIO, 'G'), 1:2), &
    'R', idxfrq(index(GNSS_PRIO, 'R'), 1:2), &
    'E', idxfrq(index(GNSS_PRIO, 'E'), 1:2), &
    'C', idxfrq(index(GNSS_PRIO, 'C'), 1:2), &
    'J', idxfrq(index(GNSS_PRIO, 'J'), 1:2)
  call exit(1)
end subroutine
