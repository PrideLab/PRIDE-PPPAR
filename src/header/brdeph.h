!!
!! purpose : GNSS broadcast format
type brdeph
  character*3 svn
  real*8 a0,a1,a2,aode,crs,dn,m0,cuc,e,cus,roota
  real*8 toe,cic,omega0,cis,i0,crc,omega,omegadot
  real*8 idot,resvd0,week,resvd1,accu,hlth,tgd,aodc
  integer*4 jd
  real*8 sod
  ! R
  real*8    pos(3),vel(3),acc(3) ! m,m/s,m/s^2
  integer*4 bn,fn,ae
end type
