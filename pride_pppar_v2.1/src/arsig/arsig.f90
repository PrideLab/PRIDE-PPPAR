!
!! arsig.f90
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
!
program arsig
  implicit none
  include '../header/const.h'
  include 'ambssd.h'
  include 'invnor.h'
  include 'arscfg.h'
  include 'ambsit.h'

  type(arscfg) ACF
  type(ambsit) AS
  type(ambssd) ASD(MAXSD_SIT), SD(MAXSD_SIT)      ! Single-Difference Ambiguities
!! Q matrix
  type(pest) PM(3 + MAXOW_ST)
  type(invm) QN

  integer*4 i, j, k, tfx, tsd
!
!! read arguments & configure options
  call get_arsig_args(ACF)
!
!! read one-way ambiguities
  if (ACF%lsearch) then
    call read_invnormal(ACF, PM, QN, AS)
  else
    call read_ambiguity(ACF, AS)
  endif
  if (ACF%fcbnprn .eq. 0) then
    write(*,'(a)') '***ERROR(arsig): no resolvable ambiguities '
    call exit(1)
  endif
!
!! define all the satellite pairs per station
  call define_sat_pairs(ACF, AS, ASD)
!
!! set narrowlane sigma for each single-difference observation
  do i = 1, AS%nsd
    if (ACF%lsearch) then
      j = max(AS%ipt(ASD(i)%ipt(1)), AS%ipt(ASD(i)%ipt(2)))
      k = min(AS%ipt(ASD(i)%ipt(1)), AS%ipt(ASD(i)%ipt(2)))
      ASD(i)%snl=dsqrt((QN%invx(j, j) + QN%invx(k, k) - 2.d0*QN%invx(j, k))*QN%vtpv/QN%frdm)
      if(ASD(i)%sys .eq. 'G')then
        ASD(i)%snl = (FREQ1_G+FREQ2_G)/VLIGHT*ASD(i)%snl
      elseif(ASD(i)%sys .eq. 'E')then
        ASD(i)%snl = (FREQ1_E+FREQ2_E)/VLIGHT*ASD(i)%snl
      elseif(ASD(i)%sys .eq. 'C' .or. ASD(i)%sys .eq. '3')then
        ASD(i)%snl = (FREQ1_C+FREQ2_C)/VLIGHT*ASD(i)%snl
      elseif(ASD(i)%sys .eq. 'J')then
        ASD(i)%snl = (FREQ1_J+FREQ2_J)/VLIGHT*ASD(i)%snl
      endif
    else
      ASD(i)%snl = 0.05d0
    endif
  enddo
!
!! fix widelane ambiguities
  call fixamb_rounding(ACF, AS, ASD(1))
!
!! find indenpendent satellite pairs for each station
  call find_indep(QN%indp, ACF, SD, AS, ASD)
!
!! map Q matrix
  QN%nfix   = 0
  if (ACF%lsearch .and. QN%indp .gt. 0) then
    QN%ndam = QN%indp
!! mapping
    call map_invnormal(SD, QN, QN%invx)
!! fix ambiguities sequentially
    call fixamb_search(SD, PM, QN, QN%invx, ACF%maxdel, ACF%minsav, ACF%chisq, ACF%ratio)
!! fixed SD
    if (QN%nfix .ne. 0) then
      do i = QN%ndam + 1, QN%ndam + QN%nfix
        SD(i)%id = 1
        SD(i)%ipt(1) = PM(SD(i)%ipt(1))%psat
        SD(i)%ipt(2) = PM(SD(i)%ipt(2))%psat
      enddo
      if (QN%nxyz .gt. 0) then
        write (*, '(a,3f15.4)') 'Fixed position(XYZ): ',(PM(i)%xest, i=1,QN%nxyz)
      endif
    endif
    write (*, '(a,2i3,f6.1,a1)') 'Narrow-lane AR(Nfix/Nall/%): ', QN%nfix, QN%indp, QN%nfix*1.d2/QN%indp, '%'
  endif
!
!! write constraint file
  if (ACF%lsearch) then
    call write_ambcon(QN%nfix, SD(QN%ndam + 1), ACF, AS)
    deallocate (QN%idq)
    deallocate (QN%invx)
  else
    call write_ambcon(QN%indp, SD, ACF, AS)
  endif

end program
