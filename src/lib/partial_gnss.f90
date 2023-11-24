!
!! partial_gnss_G.f90
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
!!! purpose  : observation equaiton of a GNSS observation
!!
!! parameter: npar     -- number of parameters of the station
!!            isat     -- satellite index in the obs. array
!!            ur2s     -- unit vector from rec. to satellite or partial of range wrt satellite position.
!!            drate    -- range rate
!!            parname  -- name of the local parameters
!!            ltog     -- index of local parameter in global parameter table, set in srif_init.
!!            xrot     -- rotation matrix for earthfixed system to inertial system
!!            trpart  -- partial of range wrt ztd and htg
!!            amat     -- observation equaiton according to local parameter table, with ltog to be
!!                        mapped to global parameter table
!!
!!
subroutine partial_gnss(npar, ur2s, drate, parname, ltog, xrot, trpart, amat, isys)
  implicit none
  include '../header/const.h'

  integer*4 npar, ltog(MAXPAR_STA)
  character*(*) parname(MAXPAR_STA)
  real*8 ur2s(3), drate, xrot(3, 3), trpart(1:*), amat(MAXPAR_STA)
!
!! local variables
  logical*1 lfirst
  integer*4 i, ipar
  character*1 isys
!
!! local parameters list see srif_init
  do ipar = 1, MAXPAR_STA
    amat(ipar) = 0.d0
  enddo

  ipar = 0
  do while (ipar .lt. npar)
    ipar = ipar + 1
    if (ltog(ipar) .eq. 0) cycle
!
!! station coordinates
    if (parname(ipar) (1:5) .eq. 'STAPX') then
      do i = 1, 3
        amat(ipar + i - 1) = -(xrot(1, i)*ur2s(1) + xrot(2, i)*ur2s(2) + xrot(3, i)*ur2s(3))
      enddo
      ipar = ipar + 2
!
!! atmospheric delay
    else if (parname(ipar) (1:3) .eq. 'ZTD') then
      amat(ipar) = trpart(1)
!
!! horizontal troposphere gradients
    else if (parname(ipar) (1:4) .eq. 'HTGC') then
      amat(ipar) = trpart(2)
    else if (parname(ipar) (1:4) .eq. 'HTGS') then
      amat(ipar) = trpart(3)
!
!! receiver clock epoch
    else if(parname(ipar)(1:8).eq.'RECCLK_'//isys) then
      amat(ipar)=1.d0-drate
!
!! ambiguity parameters
    else if (parname(ipar) (1:4) .eq. 'AMBC') then
      amat(ipar) = 1.d0
    else

    endif
  enddo

  return
end
