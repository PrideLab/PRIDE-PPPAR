!
!! ef2int.f90
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
!! Contributor: Jianghui Geng, Shuyin Mao, Songfeng Yang
!! 
!!
!!
!! purpose  : rotation from earth-fixed system to inertial system (IERS 2010)
!! parameter:
!!    input : erpfil -- ERP file
!!            jd,sod -- GPS time
!!    output: mate2j -- rotation matrix
!!            gast   -- GAST angle (radian)
!!            xpole,ypole -- polar motion (arcsec)
!!
!
subroutine ef2int(erpfil,jd,sod,mate2j,rmte2j,gast,xpole,ypole)
implicit none
include '../header/const.h'

integer*4 jd
real*8 sod,gast,xpole,ypole,mate2j(3,3),rmte2j(3,3)
character*(*) erpfil
!
!! local
integer*4 jd_tdt,jd_ut1
real*8    sod_tdt,sod_ut1,dpsi,deps
real*8    xhelp(2),taiut1r,dut1_tide,dut1_libr,dut1,dlod,domega,dels,eop(3),pm_libr(2)
real*8    rbpn(3,3),reot(3,3),rpom(3,3)
!
!! fuction called
real*8    iau_GST06,iau_SP00

!
!! time second to arc radian
call timinc(jd,sod,GPSTDT,jd_tdt,sod_tdt)
!
!! read IGS ERP
!! taiut1r: TAI-UT1R
call read_igserp(erpfil,jd,sod,taiut1r,xhelp)
xpole=xhelp(1)
ypole=xhelp(2)
!
!! effect of tidal deformation (zonal tides) on Earth's rotation (5 d to 18.6 y)
call RG_ZONT2((jd_tdt+sod_tdt/86400.d0+MJD2JD-2451545.d0)/36525.d0,dut1_tide,dlod,domega)
!
!! diurnal and semi-diurnal variations due to ocean tides affecting both polar motion and ut1
call ORTHO_EOP(jd_tdt+sod_tdt/86400.d0,eop)
eop=eop*1.d-6
!
!! libration in polar motion
call PMSDNUT2(jd_tdt+sod_tdt/86400.d0,pm_libr)
pm_libr=pm_libr*1.d-6
!
!! libration in ut1
call UTLIBR(jd_tdt+sod_tdt/86400.d0,dut1_libr,dlod)
dut1_libr=dut1_libr*1.d-6
!
!! corrected polar motion and ut1
xpole=xpole+eop(1)+pm_libr(1)    ! in arcseconds
ypole=ypole+eop(2)+pm_libr(2)    ! in arcseconds
dut1 =dut1_tide+eop(3)+dut1_libr ! in seconds
call timinc(jd,sod,19.d0-taiut1r+dut1,jd_ut1,sod_ut1)
!
!! IAU 2006/2000A nutation components
call iau_NUT06A(MJD2JD,jd_tdt+sod_tdt/86400.d0,dpsi,deps)
!
!! GCRS-to-true matrix, Q matrix
call iau_PN06(MJD2JD,jd_tdt+sod_tdt/86400.d0,dpsi,deps,rbpn)
!
!! GST corresponding to the current UT1 and generate R matrix (reot)
gast=iau_GST06(MJD2JD,jd_ut1+sod_ut1/86400.d0,MJD2JD,jd_tdt+sod_tdt/86400.d0,rbpn)
!
!! polar motion matrix: W matrix 
dels=iau_SP00(jd_tdt+sod_tdt/86400.d0,MJD2JD)
call iau_POM00(xpole/3600.d0/180.d0*pi,ypole/3600.d0/180.d0*pi,dels,rpom)
!
!!! final transformation matrix: GCRS-ITRS
call iau_C2TEQX(rbpn,gast,rpom,domega,mate2j,rmte2j)
call iau_TR(mate2j,mate2j)
call iau_TR(rmte2j,rmte2j)

return
end
