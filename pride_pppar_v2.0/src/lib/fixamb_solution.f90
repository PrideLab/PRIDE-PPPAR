!
!! fixamb_solution.f90
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
!! purpose  : fix MD ambiguities & update solutions
!! parameter:
!!    input : MD -- ambguity difference struct
!!    output: PM -- parameter struct
!!            QN%-- inversed normal matrix
!
subroutine fixamb_solution(MD, PM, QN, invx, max_del, min_sav, max_chisq, min_ratio)
  implicit none
  include '../header/const.h'
  include '../header/difamb.h'
  include '../header/invnor.h'

  type(difamb) MD(1:*)
  type(pest) PM(1:*)
  type(invm) QN
  integer*4 max_del, min_sav
  real*8 max_chisq, min_ratio, invx(1:*)
!
!! local
  integer*4 i, j, k, ntot, ierr
  real*8 dump, disall(2), chisq, ratio
  real*8, pointer :: bias(:), q22(:), q22h(:, :), q21h(:, :), q12h(:, :), qhlp(:, :)
!
  QN%nfix = 0
  QN%nfix_G = 0
  QN%nfix_E = 0
  QN%nfix_C = 0
  QN%nfix_3 = 0
  QN%nfix_J = 0
  allocate (bias(QN%indp + QN%nxyz))
  allocate (q22(QN%indp*(QN%indp + 1)/2), stat=ierr)
  if (ierr .ne. 0) then
    write (*, '(a)') '***ERROR(fixamb_solution): memory allocation q22 '
    call exit(1)
  endif
  if (QN%ndam .gt. 0) then
    do i = 1, QN%ndam
      MD(i)%id = 0
    enddo
    QN%ncad = QN%ndam
    QN%ncad_G = QN%ndam_G
    QN%ncad_E = QN%ndam_E
    QN%ncad_C = QN%ndam_C
    QN%ncad_3 = QN%ndam_3
    QN%ncad_J = QN%ndam_J
    call sort_invx(MD, QN, invx, bias, q22)
    call ambslv(QN%ncad, q22, bias, disall)
    chisq = (disall(1) + (QN%vtpv*25.d0))/(QN%frdm + QN%ncad)/(QN%vtpv*25.d0)*QN%frdm
    ratio = disall(2)/disall(1)
    write (*, '(a,i4,2f10.3)') 'Bias/Chi/Ratio(Whole): ', QN%ncad, chisq, ratio
!
!! whether the whole can be resolved
    if (chisq .ge. max_chisq .or. ratio .le. min_ratio) then
      call candid_ambi(MD, QN, invx, max_del, min_sav, max_chisq, min_ratio)
      if (QN%ncad .eq. 0) then
        write (*, '(a)') '$$$MESSAGE(fixamb_solution): no more can be fixed '; goto 100
      endif
      call sort_invx(MD, QN, invx, bias, q22)
      call ambslv(QN%ncad, q22, bias, disall)
      chisq = (disall(1) + (QN%vtpv*25.d0))/(QN%frdm + QN%ncad)/(QN%vtpv*25.d0)*QN%frdm
      ratio = disall(2)/disall(1)
    endif
    write (*, '(a,i4,2f10.3)') 'Bias/Chi/Ratio: ', QN%ncad, chisq, ratio
!
!! Q22 inversed matrix
    ntot = QN%nxyz + QN%ndam
    allocate (q22h(QN%ncad, QN%ncad))
    do j = ntot - QN%ncad + 1, ntot
      do i = ntot - QN%ncad + 1, ntot
        if (i .lt. j) then
          q22h(i - ntot + QN%ncad, j - ntot + QN%ncad) = q22h(j - ntot + QN%ncad, i - ntot + QN%ncad)
        else
          q22h(i - ntot + QN%ncad, j - ntot + QN%ncad) = invx(QN%idq(j) + i)
        endif
      enddo
    enddo
    call matinv(q22h, QN%ncad, QN%ncad, dump)
    if (dump .eq. 0.d0) then
      write (*, '(a)') '***ERROR(fixamb_solution): matrix singularity '
      call exit(1)
    endif
!
!! Q21 & Q12
    allocate (q21h(QN%ncad, ntot - QN%ncad))
    allocate (q12h(ntot - QN%ncad, QN%ncad))
    do j = 1, ntot - QN%ncad
      do i = ntot - QN%ncad + 1, ntot
        q21h(i - ntot + QN%ncad, j) = invx(QN%idq(j) + i)
        q12h(j, i - ntot + QN%ncad) = invx(QN%idq(j) + i)
      enddo
    enddo
    allocate (qhlp(ntot - QN%ncad, ntot - QN%ncad))
!
!! update
    do i = ntot - QN%ncad + 1, ntot
      dump = MD(i - QN%nxyz)%rnl
      MD(i - QN%nxyz)%rnl = bias(i - ntot + QN%ncad)               ! update fixed ambiguities
      bias(i - ntot + QN%ncad) = bias(i - ntot + QN%ncad) - dump       ! increment
    enddo
    call matmpy(q12h, q22h, q12h, ntot - QN%ncad, QN%ncad, QN%ncad)
!
!! update inversed normal matrix
    call matmpy(q12h, q21h, qhlp, ntot - QN%ncad, QN%ncad, ntot - QN%ncad)
    do j = 1, ntot - QN%ncad
      do i = j, ntot - QN%ncad
        invx(QN%idq(j) + i) = invx(QN%idq(j) + i) - qhlp(i, j)
      enddo
    enddo
!
!! update estimates & sigma
    call matmpy(q12h, bias, bias, ntot - QN%ncad, QN%ncad, 1)
    do i = 1, QN%nxyz
      PM(i)%xest = PM(i)%xest + bias(i)
    enddo

    QN%vtpv = QN%vtpv + disall(1)*0.04d0
    QN%frdm = QN%frdm + QN%ncad
    do i = 1, QN%ndam - QN%ncad
      MD(i)%rnl = MD(i)%rnl + bias(i + QN%nxyz)
      MD(i)%snl = dsqrt(invx(QN%idq(i + QN%nxyz) + i + QN%nxyz)*QN%vtpv/QN%frdm)
    enddo
!
!! clean
    QN%nfix = QN%nfix + QN%ncad
    QN%nfix_G = QN%nfix_G + QN%ncad_G
    QN%nfix_E = QN%nfix_E + QN%ncad_E
    QN%nfix_C = QN%nfix_C + QN%ncad_C
    QN%nfix_3 = QN%nfix_3 + QN%ncad_3
    QN%nfix_J = QN%nfix_J + QN%ncad_J
    QN%ndam = QN%ndam - QN%ncad
    deallocate (q22h); deallocate (q21h); deallocate (q12h); deallocate (qhlp)
  endif
100 continue
  deallocate (bias)
  deallocate (q22)

  return
end
