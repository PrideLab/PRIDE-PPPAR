!
!! arsig.f90
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
!! Ambiguity Resolution for SIngle Station
!
program arsig
  implicit none
  include '../header/const.h'
  include '../header/difamb.h'
  include '../header/invnor.h'
  include 'abfcb.h'
  include 'ambsit.h'

  type(abfcb) FCB
  type(ambsit) AS
  type(difamb) ASD(MAXSD_SIT), SD(MAXSD_SIT)      ! Single-Difference Ambiguities
!! Q matrix
  type(pest) PM(3 + MAXOW_ST)
  type(invm) QN

  integer*4 i, j, k, tfx, tsd
!
!! instruction
  write (*, '(a)') '++++++++++++++++++++++++++++++++++++'
  write (*, '(a)') 'AMBIGUITY RESOLUTION FOR PPP'
  write (*, '(a)') '++++++++++++++++++++++++++++++++++++'
!
!! read arguments & configure options
  call get_arsig_args(FCB)
!
!! read one-way ambiguities
  if (FCB%lsearch) then
    call read_invnormal(FCB, PM, QN, AS)
  else
    call read_ambiguity(FCB, AS)
  endif
!
!! define all the satellite pairs per station
  call define_sat_pairs(FCB, AS, ASD)
!
!! set narrowlane sigma for each single-difference observation
  do i = 1, AS%nsd
    if (FCB%lsearch) then
      j = max(AS%ipt(ASD(i)%ipt(1)), AS%ipt(ASD(i)%ipt(2)))
      k = min(AS%ipt(ASD(i)%ipt(1)), AS%ipt(ASD(i)%ipt(2)))
      if(ASD(i)%sys .eq. 'G')then
        ASD(i)%snl = (FREQ1_G+FREQ2_G)/VLIGHT*dsqrt((QN%invx(j, j) + QN%invx(k, k) - 2.d0*QN%invx(j, k))*QN%vtpv/QN%frdm)
      elseif(ASD(i)%sys .eq. 'E')then
        ASD(i)%snl = (FREQ1_E+FREQ2_E)/VLIGHT*dsqrt((QN%invx(j, j) + QN%invx(k, k) - 2.d0*QN%invx(j, k))*QN%vtpv/QN%frdm)
      elseif(ASD(i)%sys .eq. 'C' .or. ASD(i)%sys .eq. '3')then
        ASD(i)%snl = (FREQ1_C+FREQ2_C)/VLIGHT*dsqrt((QN%invx(j, j) + QN%invx(k, k) - 2.d0*QN%invx(j, k))*QN%vtpv/QN%frdm)
      elseif(ASD(i)%sys .eq. 'J')then
        ASD(i)%snl = (FREQ1_J+FREQ2_J)/VLIGHT*dsqrt((QN%invx(j, j) + QN%invx(k, k) - 2.d0*QN%invx(j, k))*QN%vtpv/QN%frdm)
      endif
    else
      ASD(i)%snl = 0.05d0
    endif
  enddo
!
!! fix widelane ambiguities
  call fix_ambiguity(FCB, AS, ASD(1))
!
!! find indenpendent satellite pairs for each station
  call find_indep(QN%indp,QN%indp_G,QN%indp_E,QN%indp_C,QN%indp_3,QN%indp_J, FCB, SD, AS, ASD)
!
!! map Q matrix
  QN%nfix = 0
  QN%nfix_G = 0
  QN%nfix_E = 0
  QN%nfix_C = 0
  QN%nfix_3 = 0
  QN%nfix_J = 0
  if (FCB%lsearch .and. QN%indp .gt. 0) then
    QN%ndam = QN%indp
    QN%ndam_G = QN%indp_G
    QN%ndam_E = QN%indp_E
    QN%ndam_C = QN%indp_C
    QN%ndam_3 = QN%indp_3
    QN%ndam_J = QN%indp_J
!! mapping
    call map_invnormal(SD, QN, QN%invx)
!! fix ambiguities sequentially
    call fixamb_solution(SD, PM, QN, QN%invx, FCB%maxdel, FCB%minsav, FCB%chisq, FCB%ratio)
!! fixed SD
    if (QN%nfix .ne. 0) then
      do i = QN%ndam + 1, QN%ndam + QN%nfix
        SD(i)%id = PM(SD(i)%ipt(1))%pcode(1)
        SD(i)%ipt(1) = PM(SD(i)%ipt(1))%pcode(2)
        SD(i)%ipt(2) = PM(SD(i)%ipt(2))%pcode(2)
      enddo
      do i = 1, QN%nxyz
        write (*, '(3f15.4)') PM(i)%xest
      enddo
    endif
    write (*, '(a,2i3,f6.1,a1)') 'Narrow-lane AR(search) ', QN%nfix, QN%indp, QN%nfix*1.d2/QN%indp, '%'
    if(QN%indp_G .gt. 0)then
      write (*, '(a,2i3,f6.1,a1)') 'G_Narrow-lane AR(search) ', QN%nfix_G, QN%indp_G, QN%nfix_G*1.d2/QN%indp_G, '%'
    endif
    if(QN%indp_E .gt. 0)then
      write (*, '(a,2i3,f6.1,a1)') 'E_Narrow-lane AR(search) ', QN%nfix_E, QN%indp_E, QN%nfix_E*1.d2/QN%indp_E, '%'
    endif
    if(QN%indp_C .gt. 0)then
      write (*, '(a,2i3,f6.1,a1)') 'C_Narrow-lane AR(search) ', QN%nfix_C, QN%indp_C, QN%nfix_C*1.d2/QN%indp_C, '%'
    endif
    if(QN%indp_3 .gt. 0)then
      write (*, '(a,2i3,f6.1,a1)') '3_Narrow-lane AR(search) ', QN%nfix_3, QN%indp_3, QN%nfix_3*1.d2/QN%indp_3, '%'
    endif
    if(QN%indp_J .gt. 0)then
      write (*, '(a,2i3,f6.1,a1)') 'J_Narrow-lane AR(search) ', QN%nfix_J, QN%indp_J, QN%nfix_J*1.d2/QN%indp_J, '%'
    endif
  endif
!
!! write constraint file
  if (FCB%lsearch) then
    call write_ambcon(QN%nfix, SD(QN%ndam + 1), FCB, AS)
    deallocate (QN%idq)
    deallocate (QN%invx)
  else
    call write_ambcon(QN%indp, SD, FCB, AS)
  endif

  stop
end
