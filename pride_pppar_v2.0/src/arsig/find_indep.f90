!
!! find_indep.f90
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
!! purpose  : find independent satellite pairs for ambiguity fixing
!! parameter:
!!    input : AS   -- ambiguity station
!!            ASD  -- single difference ambiguity
!!            FCB  -- fractional cycle biases
!!    output: indp -- # of independent pairs
!!            SD   -- independent ambiguity pairs
!
subroutine find_indep(indp,indp_G,indp_E,indp_C,indp_3,indp_J, FCB, SD, AS, ASD)
  implicit none
  include '../header/const.h'
  include '../header/difamb.h'
  include 'abfcb.h'
  include 'ambsit.h'

  integer*4 indp,indp_G,indp_E,indp_C,indp_3,indp_J
  type(abfcb) FCB
  type(ambsit) AS
  type(difamb) SD(1:*), ASD(MAXSD_SIT)
!
!! local
  integer*4 isit, isd

  indp = 0
  indp_G=0
  indp_E=0
  indp_C=0
  indp_3=0
  indp_J=0
  do isd = 1, AS%nsd
    if (ASD(isd)%id .eq. 0) then
      indp = indp + 1
      if(ASD(isd)%sys .eq. 'G')then
        indp_G=indp_G+1
      elseif(ASD(isd)%sys .eq. 'E')then
        indp_E=indp_E+1
      elseif(ASD(isd)%sys .eq. 'C')then
        indp_C=indp_C+1
      elseif(ASD(isd)%sys .eq. '3')then
        indp_3=indp_3+1
      elseif(ASD(isd)%sys .eq. 'J')then
        indp_J=indp_J+1
      endif
      if (indp .gt. MAXSIT*MAXSD_SIT) then
        write (*, '(a)') '***ERROR(find_indep): too many independent pairs '
        call exit(1)
      endif
      SD(indp) = ASD(isd)
      if (FCB%lsearch) then
!! pointer to inversed normal equation
        SD(indp)%ipt(1) = AS%ipt(ASD(isd)%ipt(1))
        SD(indp)%ipt(2) = AS%ipt(ASD(isd)%ipt(2))
      else
!! pointer to satellite index
        SD(indp)%id = 1
        SD(indp)%ipt(1) = AS%isat(ASD(isd)%ipt(1))
        SD(indp)%ipt(2) = AS%isat(ASD(isd)%ipt(2))
      endif
    else if (ASD(isd)%id .eq. 1) then
      if (FCB%lsearch) then
        indp = indp + 1
        if(ASD(isd)%sys .eq. 'G')then
          indp_G=indp_G+1
        elseif(ASD(isd)%sys .eq. 'E')then
          indp_E=indp_E+1
        elseif(ASD(isd)%sys .eq. 'C')then
          indp_C=indp_C+1
        elseif(ASD(isd)%sys .eq. '3')then
          indp_3=indp_3+1
        elseif(ASD(isd)%sys .eq. 'J')then
          indp_J=indp_J+1
        endif
        if (indp .gt. MAXSIT*MAXSD_SIT) then
          write (*, '(a)') '***ERROR(find_indep): too many independent pairs '
          call exit(1)
        endif
        SD(indp) = ASD(isd)
!! pointer to inversed normal equation
        SD(indp)%ipt(1) = AS%ipt(ASD(isd)%ipt(1))
        SD(indp)%ipt(2) = AS%ipt(ASD(isd)%ipt(2))
      endif
    endif
  enddo

  return
end
