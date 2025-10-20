!
!! arsig.f90
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
!! Contributor: Maorong Ge, Jianghui Geng, Songfeng Yang, Jihang Lin
!!
!
program arsig
  implicit none
  include '../header/const.h'
  include 'ambssd.h'
  include 'invnor.h'
  include 'arscfg.h'
  include 'ambsit.h'

! common
  integer*4     idxfrq(MAXSYS, 2)
  common        idxfrq
! ambiguity setting
  type(arscfg)  ACF
  type(ambsit)  AS
!! single-differenced
  type(ambssd)  ASD(MAXSD_SIT), SD(MAXSD_SIT)
!! Q matrix
  type(pest)    PM(3 + MAXOW_ST)
  type(invm)    QN
!
  integer*4     i0, i, j, k, tfx, tsd
  real*8        f1(MAXSYS), f2(MAXSYS)
!
!! read arguments & configure options
  call get_arsig_args(ACF)
!
!! frequency
  do i0 = 1, MAXSYS
    f1(i0) = FREQ_SYS(idxfrq(i0, 1), i0)
    if (f1(i0) .eq. 0.d0) goto 100
    f2(i0) = FREQ_SYS(idxfrq(i0, 2), i0)
    if (f2(i0) .eq. 0.d0) goto 100
  end do
!
!! read one-way ambiguities
  if (ACF%lsearch) then
    call read_invnormal(ACF, PM, QN, AS)
  else
    call read_ambiguity(ACF, AS)
  end if

  ACF%nobs =ACF%nobs / (nint((ACF%jd1 - ACF%jd0)*86400.0 + (ACF%sod1 -  ACF%sod0))/ACF%dintv +1)
  
  if (ACF%fcbnprn .eq. 0) then
    write (*, '(a)') '***ERROR(arsig): no resolvable ambiguities '
    call exit(1)
  end if
!
!! define all the satellite pairs per station
  call define_sat_pairs(ACF, AS, ASD)
!
!! set narrowlane sigma for each single-difference observation
  do i = 1, AS%nsd
    if (ACF%lsearch) then
      if (ASD(i)%sys .ne. '3') then
        i0 = index(GNSS_PRIO, ASD(i)%sys)
      else
        i0 = index(GNSS_PRIO, 'C')
      end if
      j = max(AS%ipt(ASD(i)%ipt(1)), AS%ipt(ASD(i)%ipt(2)))
      k = min(AS%ipt(ASD(i)%ipt(1)), AS%ipt(ASD(i)%ipt(2)))
      ASD(i)%snl = dsqrt((QN%invx(j, j) + QN%invx(k, k) - 2.d0*QN%invx(j, k))*QN%vtpv/QN%frdm)
      ASD(i)%snl = (f1(i0) + f2(i0))/VLIGHT*ASD(i)%snl
    else
      ASD(i)%snl = 0.05d0
    end if
  end do
!
!! fix widelane ambiguities
  call fixamb_rounding(ACF, AS, ASD(1))
!
!! find indenpendent satellite pairs for each station
  call find_indep(QN%indp, ACF, SD, AS, ASD)
!
!! map Q matrix
  QN%nfix = 0
  if (ACF%lsearch .and. QN%indp .gt. 0) then
    ACF%nwnfx_G = 0
    ACF%nwnfx_E = 0
    ACF%nwnfx_C = 0
    ACF%nwnfx_3 = 0
    ACF%nwnfx_J = 0
    QN%ndam = QN%indp
!! mapping
    call map_invnormal(SD, QN, QN%invx)
!! fix ambiguities sequentially
    call fixamb_search(SD, PM, QN, QN%invx, ACF%maxdel, ACF%minsav, ACF%chisq, ACF%ratio, ACF%features, ACF%nobs, ACF)
!! fixed SD
    if (QN%nfix .ne. 0) then
      do i = QN%ndam + 1, QN%ndam + QN%nfix
        SD(i)%id = 1
        SD(i)%ipt(1) = PM(SD(i)%ipt(1))%psat
        SD(i)%ipt(2) = PM(SD(i)%ipt(2))%psat
        if (SD(i)%sys .eq. 'G') ACF%nwnfx_G = ACF%nwnfx_G + 1
        if (SD(i)%sys .eq. 'E') ACF%nwnfx_E = ACF%nwnfx_E + 1
        if (SD(i)%sys .eq. 'C') ACF%nwnfx_C = ACF%nwnfx_C + 1
        if (SD(i)%sys .eq. '3') ACF%nwnfx_3 = ACF%nwnfx_3 + 1
        if (SD(i)%sys .eq. 'J') ACF%nwnfx_J = ACF%nwnfx_J + 1
      end do
      if (QN%nxyz .gt. 0) then
        write (*, '(a,3f15.4)') 'Fixed position(XYZ): ', (PM(i)%xest, i=1, QN%nxyz)
      end if
    end if
    if (ACF%ntot_G .ne. 0 .and. ACF%lverbose) then
      write (*, '(a,3i8,2(f6.1,a1,2x))') 'G_Wide/Narrow-lane FR(ind): ', ACF%nwnfx_G, ACF%nwlfx_G, ACF%ntot_G, &
        ACF%nwlfx_G*1.d2/ACF%ntot_G, '%', ACF%nwnfx_G*1.d2/ACF%nwlfx_G, '%'
    end if
    if (ACF%ntot_E .ne. 0 .and. ACF%lverbose) then
      write (*, '(a,3i8,2(f6.1,a1,2x))') 'E_Wide/Narrow-lane FR(ind): ', ACF%nwnfx_E, ACF%nwlfx_E, ACF%ntot_E, &
        ACF%nwlfx_E*1.d2/ACF%ntot_E, '%', ACF%nwnfx_E*1.d2/ACF%nwlfx_E, '%'
    end if
    if (ACF%ntot_C .ne. 0 .and. ACF%lverbose) then
      write (*, '(a,3i8,2(f6.1,a1,2x))') 'C2Wide/Narrow-lane FR(ind): ', ACF%nwnfx_C, ACF%nwlfx_C, ACF%ntot_C, &
        ACF%nwlfx_C*1.d2/ACF%ntot_C, '%', ACF%nwnfx_C*1.d2/ACF%nwlfx_C, '%'
    end if
    if (ACF%ntot_3 .ne. 0 .and. ACF%lverbose) then
      write (*, '(a,3i8,2(f6.1,a1,2x))') 'C3Wide/Narrow-lane FR(ind): ', ACF%nwnfx_3, ACF%nwlfx_3, ACF%ntot_3, &
        ACF%nwlfx_3*1.d2/ACF%ntot_3, '%', ACF%nwnfx_3*1.d2/ACF%nwlfx_3, '%'
    end if
    if (ACF%ntot_J .ne. 0 .and. ACF%lverbose) then
      write (*, '(a,3i8,2(f6.1,a1,2x))') 'mGWide/Narrow-lane FR(ind): ', ACF%nwnfx_J, ACF%nwlfx_J, ACF%ntot_J, &
        ACF%nwlfx_J*1.d2/ACF%ntot_J, '%', ACF%nwnfx_J*1.d2/ACF%nwlfx_J, '%'
    end if
    write (*, '(a,3i8,2(f6.1,a1,2x))') 'mGWide/Narrow-lane FR(ind): ', QN%nfix, QN%indp, ACF%ntot_ind, &
      QN%indp*1.d2/ACF%ntot_ind, '%', QN%nfix*1.d2/QN%indp, '%'
    ACF%nwl_ind = QN%indp
    ACF%nwn_ind = QN%nfix
  end if
!
!! write constraint file
  if (ACF%lsearch) then
    call write_ambcon(QN%nfix, SD(QN%ndam + 1), ACF, AS)
    deallocate (QN%idq)
    deallocate (QN%invx)
  else
    call write_ambcon(QN%indp, SD, ACF, AS)
  end if
  return

100 continue
  write (*, '(a,5(1x,a,2i1))') '***ERROR(lsq_add_obs): invalid frequency number:', &
    'G', idxfrq(index(GNSS_PRIO, 'G'), 1:2), &
    'R', idxfrq(index(GNSS_PRIO, 'R'), 1:2), &
    'E', idxfrq(index(GNSS_PRIO, 'E'), 1:2), &
    'C', idxfrq(index(GNSS_PRIO, 'C'), 1:2), &
    'J', idxfrq(index(GNSS_PRIO, 'J'), 1:2)
  call exit(1)
end program
