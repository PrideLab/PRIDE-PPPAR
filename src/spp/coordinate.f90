
! Coordinate function -------------------------------------------------------
! ecef to local coordinate transfromation matrix ----------------------------
! * * compute ecef to local coordinate transfromation matrix
! * * args   : real*8 *pos      I   geodetic position {lat,lon} (rad)
! * *          real*8 *E        O   ecef to local coord transformation matrix (3x3)
! * * return : none
! * * notes  : matirix stored by column-major order (fortran convention)
! * *------------------------------------------------------------------------
subroutine xyz2enu(pos, E)
implicit none
include 'file_para.h'
real*8, intent(in) :: pos(3)
real*8, intent(out) :: E(3,3)
real*8 sinp,cosp, sinl,cosl
sinp=dsin(pos(1)); cosp=dcos(pos(1))
sinl=dsin(pos(2)); cosl=dcos(pos(2))
E(1,1)=-sinl;      E(1,2)=cosl;       E(1,3)=0.d0
E(2,1)=-sinp*cosl; E(2,2)=-sinp*sinl; E(2,3)=cosp
E(3,1)=cosp*cosl;  E(3,2)=cosp*sinl;  E(3,3)=sinp
end subroutine

! transform ecef vector to local tangental coordinate -----------------------
! * * transform ecef vector to local tangental coordinate
! * * args   : real*8 *pos      I   geodetic position {lat,lon} (rad)
! * *          real*8 *r        I   vector in ecef coordinate {x,y,z}
! * *          real*8 *e        O   vector in local tangental coordinate {e,n,u}
! * * return : none
! * *------------------------------------------------------------------------
subroutine ecef2enu(pos, r, e)
implicit none
include 'file_para.h'
real*8, intent(in) :: pos(3)
real*8, intent(in) :: r(3,1)
real*8, intent(out) :: e(3,1)
real*8 E1(3,3)
call xyz2enu(pos,E1)
call matmul2('NN',3,1,3,E1,r,e)
end subroutine

! transform ecef to geodetic postion ----------------------------------------
! * * transform ecef position to geodetic position
! * * args   : real*8 *r        I   ecef position {x,y,z} (m)
! * *          real*8 *pos      O   geodetic position {lat,lon,h} (rad,m)
! * * return : none
! * * notes  : WGS84, ellipsoidal height
! * *------------------------------------------------------------------------
subroutine ecef2pos(r, pos)
implicit none
include 'file_para.h'
real*8, intent(in) :: r(3)
real*8, intent(out) :: pos(3)
real*8 e2,r2,z,zk,v,sinp
real*8, external :: dot,rcond
e2=FE_WGS84*(2.d0-FE_WGS84)
r2=dot(r,r,2)
v=RE_WGS84
z=r(3); zk=0.d0
do while(dabs(z-zk)>=1d-4)
    zk=z
    sinp=z/dsqrt(r2+z*z)
    v=RE_WGS84/dsqrt(1.d0-e2*sinp*sinp)
    z=r(3)+v*e2*sinp
enddo
pos(1)=rcond(r2>1d-12,datan(z/dsqrt(r2)),rcond(r(3)>0.d0,PI/2.d0,-PI/2.d0))
pos(2)=rcond(r2>1d-12,datan2(r(2),r(1)),0.d0)
pos(3)=dsqrt(r2+z*z)-v
end subroutine

! transform geodetic to ecef position ---------------------------------------
! transform geodetic position to ecef position
! args   : real*8 *pos      I   geodetic position {lat,lon,h} (rad,m)
!          real*8 *r        O   ecef position {x,y,z} (m)
! return : none
! notes  : WGS84, ellipsoidal height
!----------------------------------------------------------------------------
subroutine pos2ecef(pos, r)
implicit none
include 'file_para.h'
real*8, intent(in) :: pos(3)
real*8, intent(out) :: r(3)
real*8 sinp, cosp, sinl, cosl
real*8 e2, v
sinp=dsin(pos(1)); cosp=dcos(pos(1))
sinl=dsin(pos(2)); cosl=dcos(pos(2))
e2=FE_WGS84*(2.d0-FE_WGS84)
v =RE_WGS84/dsqrt(1.d0-e2*sinp*sinp)
r(1)=(v+pos(3))*cosp*cosl
r(2)=(v+pos(3))*cosp*sinl
r(3)=(v*(1.d0-e2)+pos(3))*sinp
end subroutine

! geometric distance --------------------------------------------------------
! compute geometric distance and receiver-to-satellite unit vector
! args   : real*8 *rs       I   satellilte position (ecef at transmission) (m)
!          real*8 *rr       I   receiver position (ecef at reception) (m)
!          real*8 *e        O   line-of-sight vector (ecef)
! return : geometric distance (m) (0>:error/no satellite position)
! notes  : distance includes sagnac effect correction
!----------------------------------------------------------------------------
subroutine geodist(rs, rr, e, dist)
implicit none
include 'file_para.h'
real*8, intent(in) :: rs(3),rr(3)
real*8, intent(out) :: e(3),dist
real*8 r
real*8, external :: norm
integer*4 i
if(norm(rs,3)<RE_WGS84)then
    dist=-1.d0; return
endif
e=rs-rr
r=norm(e,3)
e=e/r
dist=r+OMGE*(rs(1)*rr(2)-rs(2)*rr(1))/CLIGHT
end subroutine
