!
!! fixamb_search.f90
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
!! Contributor: Maorong Ge, Jianghui Geng, Songfeng Yang
!! 
!!
!!
!! purpose  : fix SD ambiguities & update solutions
!! parameter:
!!    input : SD -- ambguity difference struct
!!    output: PM -- parameter struct
!!            QN%-- inversed normal matrix
!
subroutine fixamb_search(SD, PM, QN, invx, max_del, min_sav, max_chisq, min_ratio, features, nobs, ACF) 
  use iso_c_binding
  implicit none
  include '../header/const.h'
  include 'ambssd.h'
  include 'invnor.h'
  include 'arscfg.h'

  type(ambssd) SD(1:*)
  type(pest) PM(1:*)
  type(invm) QN
  type(arscfg)  ACF

  integer*4 max_del, min_sav
  real*8 max_chisq, min_ratio, invx(1:*) 
  real*8 features(8), nobs
  
  ! interface to C xgb function
  interface
     subroutine callPyXgb(features, flag, prob) bind(C, name="callPyXgb")
        use iso_c_binding
        real(c_double), intent(inout) :: features(8)
        integer(c_int), intent(out)   :: flag
        real(c_double), intent(out)   :: prob
     end subroutine callPyXgb
  end interface
  
  ! local
  integer*4 i, j, k, ntot, ierr
  real*8 dump, disall(2), chisq, ratio  
  integer(c_int) :: flag
  real(c_double) :: prob
  character(len=3) :: label
  ! pointers
  real*8, pointer :: bias(:), q22(:), q22h(:, :), q21h(:, :), q12h(:, :), qhlp(:, :)
  real*8, pointer :: bias_p(:)

  features = 0.D0
  QN%nfix = 0
  allocate (bias(QN%indp + QN%nxyz))
  allocate (bias_p(QN%indp + QN%nxyz))
  allocate (q22(QN%indp*(QN%indp + 1)/2), stat=ierr)

  if (ierr .ne. 0) then
    write (*, '(a)') '***ERROR(fixamb_search): memory allocation q22 '
    call exit(1)
  endif

  if (QN%ndam .gt. 0) then
    do i = 1, QN%ndam
      SD(i)%id = 0
    enddo
    QN%ncad = QN%ndam
    call sort_invx(SD, QN, invx, bias, q22)
    do i = 1, QN%indp + QN%nxyz
      bias_p(i) = bias(i)
    enddo
    QN%ncad = QN%ndam
    call ambslv(QN%ncad, q22, bias, disall)
    chisq = (disall(1) + (QN%vtpv*25.d0))/(QN%frdm + QN%ncad)/(QN%vtpv*25.d0)*QN%frdm
    ratio = disall(2)/disall(1)
    write (*, '(/,a21,a4,2a10)') 'LAMBDA search: ','#AMB','Chi','Ratio'
    write (*, '(a,i4,2f10.3)') 'Narrow-lane AR(all): ', QN%ncad, chisq, ratio

!
!! partial AR
    if (ACF%lambsvm) then
        call cal_features(QN%ncad, q22, bias_p, chisq, ratio, nobs, features)
        call callPyXgb(features, flag, prob)
        label = merge('SUC','FAL', flag .eq. 1)
        write(*,'(a,a,f10.3)') 'AI validation(flag, prob): ', label, prob
    endif
    if ((.not. ACF%lambsvm .and. (chisq .ge. max_chisq .or. ratio .le. min_ratio)) .or. (ACF%lambsvm .and. flag .eq. 0)) then
      call candid_amb(SD, QN, ACF%lambsvm, nobs, invx, max_del, min_sav, max_chisq, min_ratio)
      if (QN%ncad .eq. 0) then
        write (*, '(a)') '$$$MESSAGE(fixamb_search): no more can be fixed '
        goto 100
      endif
      call sort_invx(SD, QN, invx, bias, q22)
      do i = 1, QN%ncad
        bias_p(i) = bias(i)
      enddo
      call ambslv(QN%ncad, q22, bias, disall)
      chisq = (disall(1) + (QN%vtpv*25.d0))/(QN%frdm + QN%ncad)/(QN%vtpv*25.d0)*QN%frdm
      ratio = disall(2)/disall(1)
      if (ACF%lambsvm) then
          call cal_features(QN%ncad, q22, bias_p, chisq, ratio, nobs, features)
          call callPyXgb(features, flag, prob)
          label = merge('SUC','FAL', flag .eq. 1)
          write(*,'(a,a,f10.3)') 'AI validation(flag, prob): ', label, prob
      endif
    endif
 
    if(ACF%lambsvm) then
        do i = 1, QN%ncad
          bias(i) = bias_p(i)
        enddo
    endif
    
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
      write (*, '(a)') '***ERROR(fixamb_search): matrix singularity '
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
      dump = SD(i - QN%nxyz)%rnl
      SD(i - QN%nxyz)%rnl = bias(i - ntot + QN%ncad)               ! update fixed ambiguities
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
      SD(i)%rnl = SD(i)%rnl + bias(i + QN%nxyz)
      SD(i)%snl = dsqrt(invx(QN%idq(i + QN%nxyz) + i + QN%nxyz)*QN%vtpv/QN%frdm)
    enddo
!
!! clean
    QN%nfix = QN%nfix + QN%ncad
    QN%ndam = QN%ndam - QN%ncad
    deallocate (q22h); deallocate (q21h); deallocate (q12h); deallocate (qhlp)
  endif
100 continue
  deallocate (bias)
  deallocate (bias_p)
  deallocate (q22)

  return
end
