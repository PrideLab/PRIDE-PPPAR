
! Rinex relevant function ---------------------------------------------------
! confirm the type of rnxfile
! rnxtype ==> 0-none, 1-stdnav(long),  2-stdobs(long)
!                     3-stdnav(short), 4-stdobs(short)
!                     5-stdnav(long,  hourly), 6-stdobs(long,  hourly)
!                     7-stdnav(short, hourly), 8-stdobs(short, hourly)
integer*4 function getrnxtyp(flname)
implicit none
include 'file_para.h'
character(*), intent(in) :: flname
integer*4 :: len1
len1=len_trim(adjustl(flname))
getrnxtyp=0
if(len1==38)then
    if(index(flname,"O.rnx")/=0) getrnxtyp=2
    if(index(flname,"_01H_")/=0) getrnxtyp=6
elseif(len1==34)then
    if(index(flname,"N.rnx")/=0) getrnxtyp=1
    if(index(flname,"_01H_")/=0) getrnxtyp=5
elseif(len1==12)then
    if(flname(8:9)=="0.")then
        if(flname(12:12)=="o" .or. flname(12:12)=="O") getrnxtyp=4
        if(flname(12:12)=="n" .or. flname(12:12)=="N") getrnxtyp=3
        if(flname(12:12)=="p" .or. flname(12:12)=="P") getrnxtyp=3
    elseif((flname(8:8)>='a' .and. flname(8:8)<='x') .or. &
           (flname(8:8)>='A' .and. flname(8:8)<='X') .and. flname(9:9)==".")then
        if(flname(12:12)=="o" .or. flname(12:12)=="O") getrnxtyp=8
        if(flname(12:12)=="n" .or. flname(12:12)=="N") getrnxtyp=7
        if(flname(12:12)=="p" .or. flname(12:12)=="P") getrnxtyp=7
    endif
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
if(ityp==1 .or. ityp==2 .or. ityp==5 .or. ityp==6)then
    read(flname(13:16),*,iostat=info) iyear
    read(flname(17:19),*,iostat=info) idoy
    getrnxmjd=ydoy2mjd_s(iyear, idoy)
elseif(ityp==3 .or. ityp==4 .or. ityp==7 .or. ityp==8)then
    read(flname(10:11),*,iostat=info) iyear
    read(flname(5:7)  ,*,iostat=info) idoy
    if(0<=iyear.and.iyear<=70) iyear=iyear+2000  ! 1971~2070
    if(70<iyear.and.iyear<=99) iyear=iyear+1900
    getrnxmjd=ydoy2mjd_s(iyear, idoy)
endif
end function

! get hour of rnxname
integer*4 function getrnxhour(flname)
implicit none
include 'file_para.h'
character(*), intent(in) :: flname
integer*4 ityp, info
integer*4, external :: getrnxtyp
ityp=getrnxtyp(flname)
getrnxhour=0
if(ityp==1 .or. ityp==2 .or. ityp==5 .or. ityp==6)then
    read(flname(20:21),*,iostat=info) getrnxhour
elseif(ityp==3 .or. ityp==4 .or. ityp==7 .or. ityp==8)then
    getrnxhour=iachar(flname(8:8))-97
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
integer*4 iyear2, idoy2, ihour1, ihour2
integer*4, external :: getrnxtyp, getrnxmjd, getrnxhour

flname2=''
ityp=getrnxtyp(flname)

if(ityp<=4)then
    mjd1=getrnxmjd(flname)
    if(mjd1==0) return
    mjd2=mjd1+offset1
    if(mjd2<=0) mjd2=0
    call mjd2ydoy_s(mjd2, iyear2, idoy2)
    flname2=flname
    if(ityp==1 .or. ityp==2)then
        write(flname2(13:19),"(I4.4,I3.3)") iyear2, idoy2
    elseif(ityp==3 .or. ityp==4)then
        write(flname2(10:11),"(I2.2)") mod(iyear2,100)
        write(flname2(5:7)  ,"(I3.3)") idoy2
    endif
elseif(ityp<=8)then
    ihour1=getrnxhour(flname)
    mjd1=getrnxmjd(flname)
    if(mjd1==0) return
    ihour2=ihour1+offset1
    mjd2=mjd1+ihour2/24
    if(mjd2<=0) mjd2=0
    call mjd2ydoy_s(mjd2, iyear2, idoy2)
    ihour2=ihour2-(mjd2-mjd1)*24
    flname2=flname
    if(ityp==5 .or. ityp==6)then
        write(flname2(13:21),"(I4.4,I3.3,I2.2)") iyear2, idoy2, ihour2
    elseif(ityp==7 .or. ityp==8)then
        write(flname2(10:11),"(I2.2)") mod(iyear2,100)
        write(flname2(5:7)  ,"(I3.3)") idoy2
        flname2(8:8)=achar(ihour2+97)
    endif
endif
end subroutine
