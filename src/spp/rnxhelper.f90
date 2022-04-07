
! Rinex relevant function ---------------------------------------------------
! confirm the type of rnxfile
! rnxtype ==> 0-none, 1-stdnav(long),  2-stdobs(long)
!                     3-stdnav(short), 4-stdobs(short)
integer*4 function getrnxtyp(flname)
implicit none
include 'file_para.h'
character(*), intent(in) :: flname
integer*4 :: len1
len1=len_trim(adjustl(flname))
getrnxtyp=0
if(len1==38)then
    if(index(flname,"O.rnx")/=0) getrnxtyp=2
elseif(len1==34)then
    if(index(flname,"N.rnx")/=0) getrnxtyp=1
elseif(len1==12)then
    if(flname(8:9)=="0." .and. flname(12:12)=="o") getrnxtyp=4
    if(flname(8:9)=="0." .and. flname(12:12)=="n") getrnxtyp=3
    if(flname(8:9)=="0." .and. flname(12:12)=="p") getrnxtyp=3
endif
end function

! get mjd of rnxname
! mjd ==> 0-error, 1-normal
integer*4 function getrnxmjd(flname)
implicit none
include 'file_para.h'
character(*), intent(in) :: flname
integer*4 ityp, info, iyear, idoy
integer*4, external :: getrnxtyp, ydoy2mjd_s
ityp=getrnxtyp(flname)
getrnxmjd=0  ! if(ityp==0)
if(ityp==1 .or. ityp==2)then
    read(flname(13:16),*,iostat=info) iyear
    read(flname(17:19),*,iostat=info) idoy
    getrnxmjd=ydoy2mjd_s(iyear, idoy)
elseif(ityp==3 .or. ityp==4)then
    read(flname(10:11),*,iostat=info) iyear
    read(flname(5:7)  ,*,iostat=info) idoy
    if(0<=iyear.and.iyear<=70) iyear=iyear+2000  ! 1971~2070
    if(70<iyear.and.iyear<=99) iyear=iyear+1900
    getrnxmjd=ydoy2mjd_s(iyear, idoy)
endif
end function

! get the next rnxname
subroutine getrnx_nname(flname, flname2, offset1)
implicit none
include 'file_para.h'
character(*), intent(in) :: flname
character(*), intent(out) :: flname2
integer*4, intent(in) :: offset1
integer*4 mjd1, mjd2, ityp
integer*4 iyear2, idoy2
integer*4, external :: getrnxtyp, getrnxmjd

flname2=""
mjd1=getrnxmjd(flname)
if(mjd1==0) return

mjd2=mjd1+offset1
if(mjd2<=0) mjd2=0
call mjd2ydoy_s(mjd2, iyear2, idoy2)

flname2=flname
ityp=getrnxtyp(flname2)
if(ityp==1 .or. ityp==2)then
    write(flname2(13:19),"(I4.4I3.3)") iyear2, idoy2
elseif(ityp==3 .or. ityp==4)then
    write(flname2(10:11),"(I2.2)") mod(iyear2,100)
    write(flname2(5:7)  ,"(I3.3)") idoy2
endif
end subroutine
