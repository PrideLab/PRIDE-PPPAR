subroutine gpt3_1(mjd,dlat,dlon,hell,it,p,T,dT,Tm,e,la,undu)

! gpt3_1.f90
!
! (c) Department of Geodesy and Geoinformation, Vienna University of
! Technology, 2018
!
!
! This subroutine determines pressure, temperature, temperature lapse rate, 
! mean temperature of the water vapor, water vapour pressure, hydrostatic 
! and wet mapping function coefficients ah and aw, water vapour decrease
! factor, geoid undulation and empirical tropospheric gradients for 
! specific sites near the earth's surface. All output values are valid for
! the specified ellipsoidal height hell.
! GPT3_1 is based on an external grid file ('gpt3_1.grd') with mean 
! values as well as sine and cosine amplitudes for the annual and
! semiannual variation of the coefficients.
! Because the .grd file has to be opened anew every time this function is 
! called, the process is fairly time-consuming for a longer set of 
! mjd's. 
!
! In order to get the respective mapping function, the ah and aw coefficients 
! must be input to vmf3_ht.f90.
! In order to get the respective zenith hydrostatic delay, pressure p must be
! input to saasthyd.f.
! In order to get an approximate value for the respective zenith wet delay,
! water vapor pressure e, mean temperature Tm and water vapor decrease factor
! la must be input to asknewet.f.
!
!
! Reference:
! D. Landskron, J. Boehm (2018), VMF3/GPT3: Refined Discrete and Empirical Troposphere Mapping Functions, 
! J Geod (2018) 92: 349., doi: 10.1007/s00190-017-1066-2. 
! Download at: https://link.springer.com/content/pdf/10.1007%2Fs00190-017-1066-2.pdf
!
!
! Input parameters:
!
! mjd:   modified Julian date (scalar, only one epoch per call is possible)
! dlat:   ellipsoidal latitude in radians [-pi/2:+pi/2] (vector)
! dlon:   longitude in radians [-pi:pi] or [0:2pi] (vector)
! hell: ellipsoidal height in m (vector)
! it:    case 1: no time variation but static quantities
!        case 0: with time variation (annual and semiannual terms)
! 
! Output parameters:
!
! p:    pressure in hPa (vector) 
! T:    temperature in degrees Celsius (vector)
! dT:   temperature lapse rate in degrees per km (vector)
! Tm:   mean temperature weighted with the water vapor in degrees Kelvin (vector) 
! e:    water vapour pressure in hPa (vector)
! la:   water vapour decrease factor (vector)
! undu: geoid undulation in m (vector)
!
!
! File created by Daniel Landskron, 2018/02/22
!
! ==========================================================

! Input arguments

integer, intent(in) :: it
double precision, intent(in) :: mjd
double precision, intent(in) :: dlat
double precision, intent(in) :: dlon
double precision, intent(in) :: hell


! Output arguments

double precision, intent(out) :: p
double precision, intent(out) :: T
double precision, intent(out) :: dT
double precision, intent(out) :: Tm
double precision, intent(out) :: e
double precision, intent(out) :: la
double precision, intent(out) :: undu


! Internal variables

! pi
double precision :: pi

! file open
character(len=1000) :: comment_line
integer :: i_grid, l

! variables for the gpt3_5.grd
double precision, dimension(:), allocatable :: lat_grid
double precision, dimension(:), allocatable :: lon_grid
double precision, dimension(:,:), allocatable :: pgrid
double precision, dimension(:,:), allocatable :: Tgrid
double precision, dimension(:,:), allocatable :: Qgrid
double precision, dimension(:,:), allocatable :: dTgrid
double precision, dimension(:), allocatable :: u
double precision, dimension(:), allocatable :: Hs
double precision, dimension(:,:), allocatable :: lagrid
double precision, dimension(:,:), allocatable :: Tmgrid

! variable for the conversion from mjd to doy
double precision :: sec, jd
integer :: hour, minu, day, month, year, jd_int
double precision :: aa, bb, cc, dd, ee, mm
integer, dimension(12) :: days
integer :: leapYear
double precision :: doy

! constants
double precision :: gm
double precision :: dMtr
double precision :: Rg

! coefficients for amplitudes
double precision :: cosfy
double precision :: coshy
double precision :: sinfy
double precision :: sinhy

! variables for coordinate conversions
double precision :: plon, ppod, difflon, diffpod
integer :: ilon, ipod, ilon1, ipod1
double precision :: dnpod1, dnpod2, dnlon1, dnlon2

! variable to decide whether ibilinear or nearest neighbor shall be applied
integer :: ibilinear

! variable for the indices in the .grd file
integer, dimension(4) :: indx
integer :: ix

! further variables
double precision :: hgt, redh, c, Hs1
double precision :: Tv, T0, p0, Q, e0
double precision, dimension(4) :: undul, Ql, dTl, pl, Tl, lal, Tml, el
double precision :: R1, R2

!------------------------------------

integer*4 get_valid_unit,lfn,ierr
logical*1 lfirst

data lfirst/.true./
save lfirst,gm,dMtr,Rg,pi
save pgrid,Tgrid,Qgrid,dTgrid,u,Hs,lagrid,Tmgrid

if(lfirst) then
  lfirst=.false.
  !% mean gravity in m/s**2
  gm = 9.80665d0;
  !% molar mass of dry air in kg/mol
  dMtr = 28.965d-3
  !% universal gas constant in J/K/mol
  Rg = 8.3143d0
  pi = 3.1415926535d0

  !% read gridfile
  lfn=get_valid_unit(10)
  open(lfn,file='gpt3_1.grd',status='old',iostat=ierr)
  if(ierr.ne.0) then
    write(*,'(a)') '***ERROR(gpt3_l): no gpt3_1.grd'
    call exit(1)
  endif
  !% read first comment line
  allocate(lat_grid(64800))
  allocate(lon_grid(64800))
  allocate(pgrid(64800,5))
  allocate(Tgrid(64800,5))
  allocate(Qgrid(64800,5))
  allocate(dTgrid(64800,5))
  allocate(u(64800))
  allocate(Hs(64800))
  allocate(lagrid(64800,5))
  allocate(Tmgrid(64800,5))
  read (unit=lfn, fmt=*) comment_line
  !% loop over grid points
  do i_grid = 1,64800
    !% read data line
    read(unit=lfn, fmt=*) lat_grid(i_grid), &
        lon_grid(i_grid), &
        pgrid(i_grid,1), pgrid(i_grid,2), pgrid(i_grid,3), pgrid(i_grid,4), pgrid(i_grid,5), &                  ! pressure in Pascal
        Tgrid(i_grid,1), Tgrid(i_grid,2), Tgrid(i_grid,3), Tgrid(i_grid,4), Tgrid(i_grid,5), &                  ! temperature in Kelvin
        Qgrid(i_grid,1), Qgrid(i_grid,2), Qgrid(i_grid,3), Qgrid(i_grid,4), Qgrid(i_grid,5), &                  ! specific humidity in kg/kg
        dTgrid(i_grid,1), dTgrid(i_grid,2), dTgrid(i_grid,3), dTgrid(i_grid,4), dTgrid(i_grid,5), &             ! temperature lapse rate in Kelvin/m
        u(i_grid),   &                                                                                          ! geoid undulation in m
        Hs(i_grid),   &                                                                                         ! orthometric grid height in m
        lagrid(i_grid,1), lagrid(i_grid,2), lagrid(i_grid,3), lagrid(i_grid,4), lagrid(i_grid,5), &             ! water vapor decrease factor, dimensionless
        Tmgrid(i_grid,1), Tmgrid(i_grid,2), Tmgrid(i_grid,3), Tmgrid(i_grid,4), Tmgrid(i_grid,5)             ! mean temperature in Kelvin
  enddo
  close(lfn)
  ! divide by constants in order to get the real values
  Qgrid = Qgrid/1000
  dTgrid = dTgrid/1000
endif
      
! convert mjd to doy

hour = floor((mjd-floor(mjd))*24)   ! get hours
minu = floor((((mjd-floor(mjd))*24)-hour)*60)   ! get minutes
sec = (((((mjd-floor(mjd))*24)-hour)*60)-minu)*60   ! get seconds

! change secs, min hour whose sec==60 and days, whose hour==24
if (sec==60) then
    minu = minu +1
end if
if (minu==60) then
    hour = hour +1
end if

! if hr==24, correct jd and set hour==0
jd = mjd+2400000.5d0
if (hour==24) then
    jd = jd + 1
end if

! integer Julian date
jd_int = floor(jd+0.5d0)

aa = jd_int+32044
bb = floor((4*aa+3)/146097)
cc = aa-floor((bb*146097)/4)
dd = floor((4*cc+3)/1461)
ee = cc-floor((1461*dd)/4)
mm = floor((5*ee+2)/153)

day = ee-floor((153*mm+2)/5)+1
month = mm+3-12*floor(mm/10)
year = bb*100+dd-4800+floor(mm/10)

! first check if the specified year is leap year or not (logical output)
if ( (modulo(year,4) == 0   .AND.   modulo(year,100) /= 0)   .OR.   modulo(year,400) == 0 ) then
    leapYear = 1
else
    leapYear = 0
end if

days = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
doy = sum(days(1:month-1)) + day
if (leapYear == 1   .AND.   month > 2) then
    doy = doy + 1
end if
doy = doy + mjd-floor(mjd)   ! add decimal places

!% factors for amplitudes
if(it.eq.1) then ! constant parameters
  cosfy = 0.d0
  coshy = 0.d0
  sinfy = 0.d0
  sinhy = 0.d0
else 
  cosfy = dcos(doy/365.25d0*2.d0*pi)
  coshy = dcos(doy/365.25d0*4.d0*pi)
  sinfy = dsin(doy/365.25d0*2.d0*pi)
  sinhy = dsin(doy/365.25d0*4.d0*pi)
endif
      
!% only positive longitude in degrees
if(dlon.lt.0.d0) then
  plon = (dlon + 2.d0*pi)*180.d0/pi
else
  plon = dlon*180.d0/pi
endif
!% transform to polar distance in degrees
ppod = (-dlat + pi/2.d0)*180.d0/pi 

!% find the index (line in the grid file) of the nearest point
ipod = floor(ppod+1.d0) 
ilon = floor(plon+1.d0)

!% normalized (to one) differences, can be positive or negative
diffpod = ppod - (ipod - 0.5d0)
difflon = plon - (ilon - 0.5d0)
!% added by HCY
if(ipod.eq.181) then
  ipod = 180
endif
!% added by GP
if(ilon.eq.361) then
  ilon = 1
endif
if(ilon.eq.0) then
  ilon = 360
endif

!% get the number of the corresponding line
indx(1) = (ipod - 1)*360 + ilon

!% near the poles: nearest neighbour interpolation, otherwise: ibilinear
!% with the 1 degree grid the limits are lower and upper (GP)
ibilinear = 0
if((ppod.gt.(0.5d0)).and.(ppod.lt.(179.5d0))) then 
  ibilinear = 1
endif          

!% case of nearest neighborhood
if(ibilinear.eq.0) then

  ix = indx(1);
        
  !% transforming ellipsoidal height to orthometric height
  undu = u(ix)
  hgt = hell-undu
    
  !% pressure, temperature at the height of the grid
  T0 = Tgrid(ix,1) + Tgrid(ix,2)*cosfy + Tgrid(ix,3)*sinfy+ Tgrid(ix,4)*coshy + Tgrid(ix,5)*sinhy
  p0 = pgrid(ix,1) + pgrid(ix,2)*cosfy + pgrid(ix,3)*sinfy+ pgrid(ix,4)*coshy + pgrid(ix,5)*sinhy
       
  !% specific humidity
  Q = Qgrid(ix,1) + Qgrid(ix,2)*cosfy + Qgrid(ix,3)*sinfy+ Qgrid(ix,4)*coshy + Qgrid(ix,5)*sinhy
          
  !% lapse rate of the temperature
  dT = dTgrid(ix,1) + dTgrid(ix,2)*cosfy + dTgrid(ix,3)*sinfy+ dTgrid(ix,4)*coshy + dTgrid(ix,5)*sinhy 

  !% station height - grid height
  redh = hgt - Hs(ix)

  !% temperature at station height in Celsius
  T = T0 + dT*redh - 273.15d0

  !% temperature lapse rate in degrees / km
  dT = dT*1000.d0

  !% virtual temperature in Kelvin
  Tv = T0*(1+0.6077d0*Q)

  c = gm*dMtr/(Rg*Tv)

  !% pressure in hPa
  p = (p0*exp(-c*redh))/100.d0
    
  !% water vapour decrease factor la - added by GP
  la = lagrid(ix,1) + lagrid(ix,2)*cosfy + lagrid(ix,3)*sinfy + lagrid(ix,4)*coshy + lagrid(ix,5)*sinhy

  !% mean temperature of the water vapor Tm - added by GP
  Tm = Tmgrid(ix,1) + Tmgrid(ix,2)*cosfy + Tmgrid(ix,3)*sinfy + Tmgrid(ix,4)*coshy + Tmgrid(ix,5)*sinhy

  !% water vapor pressure in hPa - changed by GP
  e0 = Q*p0/(0.622d0+0.378d0*Q)/100.d0  !% on the grid
  e = e0*(100.d0*p/p0)**(la+1) ! % on the station height - (14) Askne and Nordius, 1987
                
else !% ibilinear interpolation
      
  ipod1 = ipod + int(sign(1.d0,diffpod))
  ilon1 = ilon + int(sign(1.d0,difflon))
  !% changed for the 1 degree grid (GP)
  if(ilon1.eq.361) then
     ilon1 = 1
  endif
  if(ilon1.eq.0) then
    ilon1 = 360
  endif

  !% get the number of the line
  indx(2) = (ipod1 - 1)*360 + ilon   !% along same longitude
  indx(3) = (ipod  - 1)*360 + ilon1  !% along same polar distance
  indx(4) = (ipod1 - 1)*360 + ilon1  !% diagonal

  do l = 1,4
      
    !% transforming ellipsoidal height to orthometric height:
    !% Hortho = -N + Hell
    undul(l) = u(indx(l))
    hgt = hell-undul(l)

    !% pressure, temperature at the height of the grid
    T0 = Tgrid(indx(l),1) + Tgrid(indx(l),2)*cosfy + Tgrid(indx(l),3)*sinfy &
       + Tgrid(indx(l),4)*coshy + Tgrid(indx(l),5)*sinhy
    p0 = pgrid(indx(l),1) + pgrid(indx(l),2)*cosfy + pgrid(indx(l),3)*sinfy &
       + pgrid(indx(l),4)*coshy + pgrid(indx(l),5)*sinhy

    !% humidity 
    Ql(l) = Qgrid(indx(l),1) + Qgrid(indx(l),2)*cosfy + Qgrid(indx(l),3)*sinfy &
          + Qgrid(indx(l),4)*coshy + Qgrid(indx(l),5)*sinhy

    !% reduction = stationheight - gridheight
    Hs1 = Hs(indx(l))
    redh = hgt - Hs1

    !% lapse rate of the temperature in degree / m
    dTl(l) = dTgrid(indx(l),1) + dTgrid(indx(l),2)*cosfy + dTgrid(indx(l),3)*sinfy &
           + dTgrid(indx(l),4)*coshy + dTgrid(indx(l),5)*sinhy 

    !% temperature reduction to station height
    Tl(l) = T0 + dTl(l)*redh - 273.15d0

    !% virtual temperature
    Tv = T0*(1+0.6077d0*Ql(l))  
    c = gm*dMtr/(Rg*Tv)
    
    !% pressure in hPa
    pl(l) = (p0*exp(-c*redh))/100.d0
    
    !% water vapor decrease factor la - added by GP
    lal(l) = lagrid(indx(l),1) + lagrid(indx(l),2)*cosfy + lagrid(indx(l),3)*sinfy &
           + lagrid(indx(l),4)*coshy + lagrid(indx(l),5)*sinhy

    !% mean temperature of the water vapor Tm - added by GP
    Tml(l) = Tmgrid(indx(l),1) + Tmgrid(indx(l),2)*cosfy + Tmgrid(indx(l),3)*sinfy &
           + Tmgrid(indx(l),4)*coshy + Tmgrid(indx(l),5)*sinhy

    !% water vapor pressure in hPa - changed by GP
    e0 = Ql(l)*p0/(0.622d0+0.378d0*Ql(l))/100.d0 !% on the grid
    el(l) = e0*(100.d0*pl(l)/p0)**(lal(l)+1.d0)  !% on the station height  (14) Askne and Nordius, 1987
    
  enddo
          
  dnpod1 = abs(diffpod) !% distance nearer point
  dnpod2 = 1.d0 - dnpod1   !% distance to distant point
  dnlon1 = abs(difflon)
  dnlon2 = 1.d0 - dnlon1

  !% pressure
  R1 = dnpod2*pl(1)+dnpod1*pl(2)
  R2 = dnpod2*pl(3)+dnpod1*pl(4)
  p = dnlon2*R1+dnlon1*R2
    
  !% temperature
  R1 = dnpod2*Tl(1)+dnpod1*Tl(2)
  R2 = dnpod2*Tl(3)+dnpod1*Tl(4)
  T = dnlon2*R1+dnlon1*R2

  !% temperature in degree per km
  R1 = dnpod2*dTl(1)+dnpod1*dTl(2)
  R2 = dnpod2*dTl(3)+dnpod1*dTl(4)
  dT = (dnlon2*R1+dnlon1*R2)*1000.d0
    
  !% water vapor pressure in hPa - changed by GP
  R1 = dnpod2*el(1)+dnpod1*el(2)
  R2 = dnpod2*el(3)+dnpod1*el(4)
  e = dnlon2*R1+dnlon1*R2
    
  !% undulation
  R1 = dnpod2*undul(1)+dnpod1*undul(2)
  R2 = dnpod2*undul(3)+dnpod1*undul(4)
  undu = dnlon2*R1+dnlon1*R2

  !% water vapor decrease factor la - added by GP
  R1 = dnpod2*lal(1)+dnpod1*lal(2)
  R2 = dnpod2*lal(3)+dnpod1*lal(4)
  la = dnlon2*R1+dnlon1*R2

  !% mean temperature of the water vapor Tm - added by GP
  R1 = dnpod2*Tml(1)+dnpod1*Tml(2)
  R2 = dnpod2*Tml(3)+dnpod1*Tml(4)
  Tm = dnlon2*R1+dnlon1*R2
endif

return      

Entry clean_gpt3_1()
  if (allocated(lat_grid)) deallocate(lat_grid)
  if (allocated(lon_grid)) deallocate(lon_grid)
  if (allocated(pgrid)) deallocate(pgrid)
  if (allocated(Tgrid)) deallocate(Tgrid)
  if (allocated(Qgrid)) deallocate(Qgrid)
  if (allocated(dTgrid)) deallocate(dTgrid)
  if (allocated(u)) deallocate(u)
  if (allocated(Hs)) deallocate(Hs)
  if (allocated(lagrid)) deallocate(lagrid)
  if (allocated(Tmgrid)) deallocate(Tmgrid)

end subroutine
