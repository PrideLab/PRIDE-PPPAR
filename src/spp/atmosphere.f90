
! Atmospheric function ------------------------------------------------------
! ionosphere model ----------------------------------------------------------
! compute ionospheric delay by broadcast ionosphere model (klobuchar model)
! args   : gtime_t t        I   time (gpst)
!          real*8 *ion      I   iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3}
!          real*8 *pos      I   receiver position {lat,lon,h} (rad,m)
!          real*8 *azel     I   azimuth/elevation angle {az,el} (rad)
! return : ionospheric delay (L1) (m)
!----------------------------------------------------------------------------
real*8 function ionmodel(t, ion_in, pos, azel)
implicit none
include 'file_para.h'
type(gtime_t), intent(in) :: t
real*8, intent(in) :: ion_in(8),pos(3),azel(2)
real*8 :: ion_default(8)=(/&
    0.1118d-07,-0.7451d-08,-0.5961d-07, 0.1192d-06,&
    0.1167d+06,-0.2294d+06,-0.1311d+06, 0.1049d+07 /)
real*8 tt,f,psi,phi,lam,amp,per,x,ion(8),gpst
real*8, external :: norm,rcond
integer*4 week
if (pos(3)<-1d3 .or. azel(2)<=0.d0)then
    ionmodel=0.d0; return
endif
ion=ion_in
if (norm(ion,8)<=0.d0) ion=ion_default

! earth centered angle (semi-circle) 
psi=0.0137d0/(azel(2)/PI+0.11d0)-0.022d0

! subionospheric latitude/longitude (semi-circle) 
phi=pos(1)/PI+psi*dcos(azel(1))
if (phi> 0.416d0) phi= 0.416d0
if (phi<-0.416d0) phi=-0.416d0
lam=pos(2)/PI+psi*dsin(azel(1))/dcos(phi*PI)

! geomagnetic latitude (semi-circle) 
phi=phi+0.064d0*dcos((lam-1.617d0)*PI)

! local time (s) 
call time2gpst(t,week,gpst)
tt=43200.d0*lam+gpst
tt=tt-floor(tt/86400.d0)*86400.d0;  ! 0<=tt<86400 

! slant factor 
f=1.d0+16.d0*(0.53-azel(2)/PI)**3.d0

! ionospheric delay 
amp=ion(1)+phi*(ion(2)+phi*(ion(3)+phi*ion(4)))
per=ion(5)+phi*(ion(6)+phi*(ion(7)+phi*ion(8)))
amp=rcond(amp<    0.d0,    0.d0,amp)
per=rcond(per<72000.d0,72000.d0,per)
x=2.d0*PI*(tt-50400.d0)/per

ionmodel=CLIGHT*f*rcond(dabs(x)<1.57d0,5d-9+amp*(1.d0+x*x*(-0.5+x*x/24.d0)),5d-9)
end function

! ionospheric correction ----------------------------------------------------
! compute ionospheric correction
! args   : gtime_t time      I   time
!          nav_t  *nav       I   navigation data
!          integer*4 sat     I   satellite number
!          real*8 *pos       I   receiver position {lat,lon,h} (rad|m)
!          real*8 *azel      I   azimuth/elevation angle {az,el} (rad)
!          integer*4 ionoopt I   ionospheric correction option (IONOOPT_???)
!          real*8 *ion       O   ionospheric delay (L1) (m)
!          real*8 *var       O   ionospheric delay (L1) variance (m^2)
! return : status(1:ok,0:error)
!----------------------------------------------------------------------------
subroutine ionocorr(time, nav, sat, pos, azel, ionoopt, ion, var, stat)
implicit none
include 'file_para.h'
type(gtime_t), intent(in) :: time
type(nav_t), intent(in) :: nav
integer*4, intent(in) :: sat,ionoopt
real*8, intent(in) :: pos(3),azel(2)
real*8, intent(out) :: ion,var
integer*4, intent(out) :: stat
real*8, external :: ionmodel,SQR,rcond
! broadcast model 
if(ionoopt==IONOOPT_BRDC)then
    ion=ionmodel(time,nav%ion_gps,pos,azel)
    var=SQR(ion*ERR_BRDCI)
    stat=1; return
endif
ion=0.d0
var=rcond(ionoopt==IONOOPT_OFF,SQR(ERR_ION),0.d0)
stat=1
end subroutine

! troposphere model ---------------------------------------------------------
! compute tropospheric delay by standard atmosphere and saastamoinen model
! args   : gtime_t time     I   time
!          real*8 *pos      I   receiver position {lat,lon,h} (rad,m)
!          real*8 *azel     I   azimuth/elevation angle {az,el} (rad)
!          real*8 humi      I   relative humidity
! return : tropospheric delay (m)
!----------------------------------------------------------------------------
real*8 function tropmodel(time, pos, azel, humi)
implicit none
include 'file_para.h'
type(gtime_t), intent(in) :: time
real*8, intent(in) :: pos(3),azel(2),humi
real*8 :: hgt,pres,temp,e,z,trph,trpw,temp0
real*8, external :: rcond
temp0=15.d0
if (pos(3)<-100.d0 .or. 1d4<pos(3) .or. azel(2)<=0)then
    tropmodel=0.d0; return
endif

! standard atmosphere 
hgt=rcond(pos(3)<0.d0,0.d0,pos(3))
pres=1013.25d0*(1.d0-2.2557d-5*hgt)**5.2568d0
temp=temp0-6.5d-3*hgt+273.16d0
e=6.108d0*humi*dexp((17.15d0*temp-4684.d0)/(temp-38.45d0))

! saastamoninen model 
z=PI/2.d0-azel(2)
trph=0.0022768d0*pres/(1.d0-0.00266d0*dcos(2.d0*pos(1))-0.00028d0*hgt/1d3)/dcos(z)
trpw=0.002277d0*(1255.d0/temp+0.05d0)*e/dcos(z)
tropmodel=trph+trpw
end function

! tropospheric correction ---------------------------------------------------
! compute tropospheric correction
! args   : gtime_t time     I   time
!          nav_t  *nav      I   navigation data
!          real*8 *pos      I   receiver position {lat,lon,h} (rad|m)
!          real*8 *azel     I   azimuth/elevation angle {az,el} (rad)
!          integer*4 tropopt I  tropospheric correction option (TROPOPT_???)
!          real*8 *trp      O   tropospheric delay (m)
!          real*8 *var      O   tropospheric delay variance (m^2)
! return : status(1:ok,0:error)
!----------------------------------------------------------------------------
subroutine tropcorr(time, nav, pos, azel, tropopt, trp, var, stat)
implicit none
include 'file_para.h'
type(gtime_t), intent(in) :: time
type(nav_t), intent(in) :: nav
real*8, intent(in) :: pos(3),azel(2)
integer*4, intent(in) :: tropopt
real*8, intent(out) :: trp, var
integer*4, intent(out) :: stat
real*8, external :: tropmodel,SQR,rcond
! saastamoinen model 
if(tropopt==TROPOPT_SAAS .or. tropopt==TROPOPT_EST .or. tropopt==TROPOPT_ESTG)then
    trp=tropmodel(time,pos,azel,REL_HUMI)
    var=SQR(ERR_SAAS/(dsin(azel(2))+0.1d0))
    stat=1; return
endif
! no correction 
trp=0.d0
var=rcond(tropopt==TROPOPT_OFF,SQR(ERR_TROP),0.d0)
stat=1
end subroutine
