!
!! candid_ambi.f90
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
!! Contributor: Jianghui Geng, Songfeng Yang
!! 
!!
!!
!! purpose  : search ambiguity candidates. stack algorithm
!! parameter:
!!    input : QN -- normal matrix
!!            max_del -- maximum deleted ambiguities
!!            min_sav -- minimum saved ambiguiites
!!            max_chisq -- chi square test statistics
!!            min_ratio -- ratio test critical value
!!    output: MD -- multiple difference struct
!
subroutine candid_ambi(MD, QN, invx, max_del, min_sav, max_chisq, min_ratio)
  implicit none
  include '../header/difamb.h'
  include '../header/invnor.h'

  type(difamb) MD(1:*)
  type(invm) QN
  integer*4 max_del, min_sav
  real*8 max_chisq, min_ratio, invx(1:*)
!
!! local
  integer*4, pointer :: candi(:), savi(:)
  real*8, pointer :: bias(:), q22(:)
  integer*4 ncad, ndel, i, j, k, l, m
  integer*4 ncad_G,ncad_E,ncad_C,ncad_3,ncad_J,ndel_G,ndel_E,ndel_C,ndel_3,ndel_J
  real*8 disall(2), chisq, ratio, cdrto
!
!! function called
  integer*4 pointer_int

  QN%ncad = 0
  QN%ncad_G = 0
  QN%ncad_E = 0
  QN%ncad_C = 0
  QN%ncad_3 = 0
  QN%ncad_J = 0
  if (max_del .le. 0 .or. QN%ndam .le. 1 .or. QN%ndam + QN%nfix .le. min_sav) return
!
!! allocation
  allocate (candi(max_del))
  allocate (savi(max_del))
  allocate (bias(1:QN%ndam))
  allocate (q22(QN%ndam*(QN%ndam + 1)/2))
!
!! search best candidate
  ndel = 1
  ndel_G=0
  ndel_E=0
  ndel_C=0
  ndel_3=0
  ndel_J=0
  if(MD(1)%sys .eq. 'G')then
    ndel_G=1
  elseif(MD(1)%sys .eq. 'E')then
    ndel_E=1
  elseif(MD(1)%sys .eq. 'C')then
    ndel_C=1
  elseif(MD(1)%sys .eq. '3')then
    ndel_3=1
  elseif(MD(1)%sys .eq. 'J')then
    ndel_J=1
  endif
  do while (ndel .le. max_del .and. QN%ndam - ndel .gt. 1 .and. QN%ndam + QN%nfix - ndel .ge. min_sav)
    i = ndel
    ncad = QN%ndam - ndel
    ncad_G = QN%ndam_G - ndel_G
    ncad_E = QN%ndam_E - ndel_E
    ncad_C = QN%ndam_C - ndel_C
    ncad_3 = QN%ndam_3 - ndel_3
    ncad_J = QN%ndam_J - ndel_J
    candi(1:max_del) = 0
    cdrto = 0.d0
    do while (i .gt. 0)
      call sel_candi(QN%ndam, i, candi)
      if (i .gt. 0) then
        l = 0; m = 0
        do j = 1, QN%ndam
          if (pointer_int(i, candi, j) .ne. 0) cycle
          l = l + 1
          bias(l) = MD(j)%rnl
          do k = j, QN%ndam
            if (pointer_int(i, candi, k) .ne. 0) cycle
            m = m + 1
            q22(m) = invx(QN%idq(QN%nxyz + j) + QN%nxyz + k)
          enddo
        enddo
        call ambslv(ncad, q22, bias, disall)
        chisq = (disall(1) + (QN%vtpv*25.d0))/(QN%frdm + QN%ncad)/(QN%vtpv*25.d0)*QN%frdm
        ratio = disall(2)/disall(1)
        if (ratio .gt. cdrto .and. chisq .lt. max_chisq) then
          cdrto = ratio
          savi(1:ndel) = candi(1:ndel)
        endif
      endif
    enddo
    if (cdrto .gt. min_ratio) then
      do i = 1, ndel
        MD(savi(i))%id = 1
      enddo
      QN%ncad = ncad
      QN%ncad_G = ncad_G
      QN%ncad_E = ncad_E
      QN%ncad_C = ncad_C
      QN%ncad_3 = ncad_3
      QN%ncad_J = ncad_J
      ndel_G=0
      ndel_E=0
      ndel_C=0
      ndel_3=0
      ndel_J=0
      do i=1,ndel
        if(MD(i)%sys .eq. 'G')then
          ndel_G = ndel_G + 1
        elseif(MD(i)%sys .eq. 'E')then
          ndel_E = ndel_E + 1
        elseif(MD(i)%sys .eq. 'C')then
          ndel_C = ndel_C + 1
        elseif(MD(i)%sys .eq. '3')then
          ndel_3 = ndel_3 + 1
        elseif(MD(i)%sys .eq. 'J')then
          ndel_J = ndel_J + 1
        endif
      enddo
      exit
    else
      ndel = ndel + 1
      ndel_G=0
      ndel_E=0
      ndel_C=0
      ndel_3=0
      ndel_J=0
      do i=1,ndel
        if(MD(i)%sys .eq. 'G')then
          ndel_G = ndel_G + 1
        elseif(MD(i)%sys .eq. 'E')then
          ndel_E = ndel_E + 1
        elseif(MD(i)%sys .eq. 'C')then
          ndel_C = ndel_C + 1
        elseif(MD(i)%sys .eq. '3')then
          ndel_3 = ndel_3 + 1
        elseif(MD(i)%sys .eq. 'J')then
          ndel_J = ndel_J + 1
        endif
      enddo
    endif
  enddo

!
!! clean
  deallocate (candi)
  deallocate (savi)
  deallocate (bias)
  deallocate (q22)

  return
end
!--------------------------------------------------------
!! purpose  : select candidates
!! parameter:
!!    input : ntot -- highest number
!!    output: ndel -- deleted number
!!            candi -- candidate numbers
!
subroutine sel_candi(ntot, ndel, candi)
  implicit none

  integer*4 ntot, ndel, candi(1:*)
!
!! local
  integer*4 i, ic

  if (candi(1) .eq. 0) then
    do i = 1, ndel
      candi(i) = i
    enddo
    return
  endif

  ic = ndel
  do while (ic .gt. 0)
    candi(ic) = candi(ic) + 1
    if (candi(ic) .gt. ntot) then
      ic = ic - 1
      cycle
    endif
    i = ic + 1
    do while (i .le. ndel)
      candi(i) = candi(i - 1) + 1
      if (candi(i) .gt. ntot) then
        ic = ic - 1
        exit
      endif
      i = i + 1
    enddo
    if (i .gt. ndel) exit
  enddo
  if (ic .le. 0) ndel = 0

  return
end
