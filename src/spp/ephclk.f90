
! Calculate the satellite clock offset --------------------------------------
! satellite clock with broadcast ephemeris ----------------------------------
subroutine ephclk(time, teph, sat, nav, dts, stat)
implicit none
include 'file_para.h'
type(gtime_t), intent(in) :: time,teph
integer*4, intent(in) :: sat
type(nav_t), intent(in) :: nav
real*8, intent(out) :: dts
integer*4, intent(out) :: stat
type(eph_t), pointer :: eph  !,seleph
type(geph_t), pointer :: geph  !,selgeph
!external :: seleph,selgeph
real*8, external :: eph2clk,geph2clk

interface
  function seleph(time, sat, iode, nav)
  implicit none
  include 'file_para.h'
  type(gtime_t), intent(in) :: time
  integer*4, intent(in) :: sat,iode
  type(nav_t), intent(in) :: nav
  type(eph_t), pointer :: seleph
  end function

  function selgeph(time, sat, iode, nav)
  implicit none
  include 'file_para.h'
  type(gtime_t), intent(in) :: time
  integer*4, intent(in) :: sat,iode
  type(nav_t), intent(in) :: nav
  type(geph_t), pointer :: selgeph
  end function
end interface

integer*4 prn,sys
call satsys(sat,prn,sys)

if (sys==SYS_GPS .or. sys==SYS_GAL .or. sys==SYS_QZS .or. sys==SYS_CMP)then
    !eph=seleph(teph,sat,-1,nav)
    eph=>seleph(teph,sat,-1,nav)
    if (.not.associated(eph))then
        stat=0; return
    endif
    dts=eph2clk(time,eph)
elseif(sys==SYS_GLO)then
    !geph=selgeph(teph,sat,-1,nav)
    geph=>selgeph(teph,sat,-1,nav)
    if (.not.associated(geph))then
        stat=0; return
    endif
    dts=geph2clk(time,geph)
else
    stat=0; return
endif
stat=1
end subroutine

! broadcast ephemeris to satellite clock bias -------------------------------
! compute satellite clock bias with broadcast ephemeris (gps, galileo, qzss)
! args   : gtime_t time     I   time by satellite clock (gpst)
!          eph_t *eph       I   broadcast ephemeris
! return : satellite clock bias (s) without relativeity correction
! notes  : see ref (1),(7),(8)
!          satellite clock does not include relativity correction and tdg
!----------------------------------------------------------------------------
real*8 function eph2clk(time, eph)
implicit none
include 'file_para.h'
type(gtime_t),intent(in) :: time
type(eph_t), intent(in) :: eph
real*8 t
integer*4 i
real*8, external :: timediff
t=timediff(time,eph%toc)
do i=1,2
    t=t-(eph%f0+eph%f1*t+eph%f2*t*t)
enddo
eph2clk=eph%f0+eph%f1*t+eph%f2*t*t
end function

! glonass ephemeris to satellite clock bias ---------------------------------
! compute satellite clock bias with glonass ephemeris
! args   : gtime_t time     I   time by satellite clock (gpst)
!          geph_t *geph     I   glonass ephemeris
! return : satellite clock bias (s)
! notes  : see ref (2)
!----------------------------------------------------------------------------
real*8 function geph2clk(time, geph)
implicit none
include 'file_para.h'
type(gtime_t), intent(in) :: time
type(geph_t), intent(in) :: geph
real*8 t
integer*4 i
real*8, external :: timediff
t=timediff(time,geph%toe)
do i=1,2
    t=t-(-geph%taun+geph%gamn*t)
enddo
geph2clk=-geph%taun+geph%gamn*t
end function
