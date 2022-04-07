
! Date function -------------------------------------------------------------
! ydoy ==> mjd
integer*4 function ydoy2mjd_s(iyear, idoy)
implicit none
include 'file_para.h'
integer*4, intent(in) :: iyear, idoy
integer*4 imonth, iday
integer*4, external :: ymd2mjd_s
call ydoy2ymd_s(iyear, idoy, imonth, iday)
ydoy2mjd_s=ymd2mjd_s(iyear, imonth, iday)
end function

! mjd ==> ydoy
subroutine mjd2ydoy_s(jd, iyear, idoy)
implicit none
include 'file_para.h'
integer*4 jd, iyear, idoy
integer*4, external :: ymd2mjd_s
iyear = (jd + 678940)/365
idoy = jd - ymd2mjd_s(iyear, 1, 1)
do while (idoy .le. 0)
  iyear = iyear - 1
  idoy = jd - ymd2mjd_s(iyear, 1, 1) + 1
enddo
return
end subroutine

! ydoy ==> ymd
subroutine ydoy2ymd_s(iyear, idoy, imonth, iday)
implicit none
include 'file_para.h'
integer*4 iyear, idoy, imonth, iday
!! local
integer*4 days_in_month(12), id
data days_in_month/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
days_in_month(2) = 28
if (mod(iyear, 4) .eq. 0 .and. (mod(iyear, 100) .ne. 0 .or. mod(iyear, 400) .eq. 0)) days_in_month(2) = 29
id = idoy
do imonth = 1, 12
  id = id - days_in_month(imonth)
  if (id .gt. 0) cycle
  iday = id + days_in_month(imonth)
  exit
enddo
return
end subroutine

! ymd ==> mjd
integer*4 function ymd2mjd_s(iyear, imonth, iday)
implicit none
include 'file_para.h'
integer*4 imonth, iday, iyear
logical*1, external :: IsDayRight
!! local
integer*4 iyr, doy_of_month(12), modified_julday
data doy_of_month/0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334/
!! check the input data
ymd2mjd_s=-1
if(.not.IsDayRight(iyear,imonth,iday)) return  ! call exit(1)
iyr = iyear
if (imonth .le. 2) iyr = iyr - 1
modified_julday = 365*iyear - 678941 + iyr/4 - iyr/100 + iyr/400 + iday
if (imonth .ne. 0) modified_julday = modified_julday + doy_of_month(imonth)
ymd2mjd_s=modified_julday
return
end function

! convert ymd to mjd --------------------------------------------------------
integer*4 function ymd2mjd(ep)
implicit none
include 'file_para.h'
integer*4, intent(in) :: ep(3)
integer*4 year, mon, day
logical*1, external :: IsDayRight
real*8 mjd
ymd2mjd=-1; year=ep(1); mon=ep(2); day=ep(3)
if(.not.IsDayRight(year,mon,day)) return  ! call exit(1)

if(year<100) year=year+2000
if(mon<=2)then
    mon=mon+12
    year=year-1
endif
mjd=year*365.25d0-dmod(year*365.25d0,1.d0)-679006
mjd=mjd+int(30.6001d0*(mon+1))+2-int(year/100.d0)+int(year/400.d0)+day
ymd2mjd=int(mjd)
end function

! date check ----------------------------------------------------------------
logical*1 function IsDayRight(iyear, imon, iday)
implicit none
include 'file_para.h'
integer*4, intent(in) :: iyear, imon, iday
if(iyear<=0 .or. imon<=0 .or. iday <=0) goto 110
select case(imon)
case(1,3,5,7,8,10,12)
    if(iday>31) goto 110
case(4,6,9,11)
    if(iday>30) goto 110
case(2)
    if(iday>29) goto 110
case default
    goto 110
end select
if(.not.((mod(iyear,4)==0 .and. mod(iyear,100)/=0) .or. mod(iyear,400)==0))then
    if(imon==2 .and. iday>28) goto 110
endif
100 IsDayRight=.true. ; return
110 IsDayRight=.false.; return  ! write(*,"('Incorrect date format : 'I4.4' 'I2.2' 'I2.2)") iyear,imon,iday;
end function

! time check ----------------------------------------------------------------
logical*1 function IsTimeRight(ihour, imin, isec)
implicit none
include 'file_para.h'
integer*4, intent(in) :: ihour, imin
real*8, intent(in) :: isec
if(.not.(0  <=ihour.and.ihour<=23)) goto 110
if(.not.(0  <=imin .and.imin <=59)) goto 110
if(.not.(0d0<=isec .and.isec < 60d0)) goto 110
100 IsTimeRight=.true. ; return
110 IsTimeRight=.false.; return  ! write(*,"('Incorrect time format : 'I3.2' 'I2.2' 'f5.2)") ihour,imin,isec;
end function
