!
!! Constant parameter

!! frequency
real*8, parameter    :: freq1_G=1.575420d9, freq2_G=1.22760d9
real*8, parameter    :: freq1_E=1.575420d9, freq2_E=1.17645d9
real*8, parameter    :: freq1_C=1.561098d9, freq2_C=1.26852d9                         ! BDS-2 B1I&B3I
real*8, parameter    :: freq1_J=1.575420d9, freq2_J=1.22760d9

!! physical constant
real*8, parameter    :: vlight=299792458.d0, pi=3.141592654d0
real*8, parameter    :: gm=398600.4418d0, gms=1.32712442076d11, gmm=4.90280070054d3   ! km3/s-2
real*8, parameter    :: erad=6378.1366d0                                              ! km
real*8, parameter    :: omega=7.292115d-5                                             ! rad/s

!! time system
real*8, parameter    :: gpstai=19.d0, taitdt=32.184d0
real*8, parameter    :: gpstdt=gpstai+taitdt, mjd2jd=2400000.5d0
real*8, parameter    :: maxwnd=0.01d0

!! signal type
character*9, parameter :: obs_prio_G='XSLCPW'
character*9, parameter :: obs_prio_R='CP'
character*9, parameter :: obs_prio_E='CQX'
character*9, parameter :: obs_prio_C='I'
character*9, parameter :: obs_prio_J='XSLC'

!! allocation limit
integer*4, parameter :: maxtyp=36
integer*4, parameter :: maxday=100
integer*4, parameter :: maxsat_G=40, maxsat_R=30, maxsat_E=40, maxsat_C=70, maxsat_J=10
integer*4, parameter :: maxsat=maxsat_G+maxsat_R+maxsat_E+maxsat_C+maxsat_J
integer*4, parameter :: maxeph=(maxsat*24+maxsat_E*120)*maxday

!! parameter
integer*4, parameter :: maxpar_sta=maxday*5
integer*4, parameter :: maxow_st=maxsat*maxday*5
integer*4, parameter :: maxpar=maxpar_sta+maxow_st
integer*4, parameter :: maxsd_sit=maxsat*maxday*40

!! ionospheric grid file
integer*4, parameter ::  maxlat=71
integer*4, parameter ::  maxlon=73
