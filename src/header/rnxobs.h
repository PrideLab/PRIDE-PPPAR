!
!! rinex header information  
type rnxhdr 
  logical*1    lc1p1,lc2p2
  integer*4    ver,fact1,fact2,nobstyp
  integer*4    nobstyp3_sys(MAXSYS)
  integer*4    t0(5),t1(5)
  real*8       intv,x,y,z,h,e,n,t0s,t1s
  character*1  sys
  character*3  obstyp(MAXTYP)
  character*3  cobstyp3_sys(MAXTYP, MAXSYS)
  character*4  mark
  character*6  recnum,antnum
  character*20 rectyp,anttyp
end type 
!
!! rinex observation records
type rnxobr 
  integer*4    jd
  integer*4    nprn
  integer*4    itypuse(maxsat,4)
  character*3  typuse(maxsat,4),prn(maxsat)
  real*8       tsec,dtrcv
  real*8       obs(maxsat,6)
  integer*1    lli(maxsat,2) ! Lost of lock indicator (LLI) for carrier-phase observations
!
!! editing
  integer*4    lfnrhd,amb_epo,amb_tot,ava_obs,rem_obs
  character*20 rhdfil
!
!! model
  integer*2    flag(maxsat)
  integer*4    npar,ltog(maxpar_sta,maxsat)
  real*8       ztdc,zage,zdd,zwd,nhtg,ehtg
  real*8       lifamb(maxsat,2),azim(maxsat),elev(maxsat),dmap(maxsat),wmap(maxsat)
  real*8       delay(maxsat),omc(maxsat,4),var(maxsat,4),amat(maxpar_sta,maxsat)
  character*15 pname(maxpar_sta)
!
!! ion model
  real*8       ilat(maxsat)  ! IPP coordicantes in radian
  real*8       ilon(maxsat)  ! IPP coordicantes in radian
  real*8       zeni(maxsat)  ! zenith distance in radians at IPP 
  real*8       ib0(maxsat),itheta(maxsat)

  integer*4    glschn(MAXSAT_R)
end type
