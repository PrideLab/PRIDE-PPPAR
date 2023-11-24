!!
!! parameter definition
!!
type prmt
  character*15 pname                 !! parameter name  
  character*1 ptype                  !! parameter type
  integer*4 ipt                      !! pointer to normal equation  
  integer*4 iepo                     !! number of epochs 
  integer*4 iobs,iobs_G,iobs_R,iobs_E,iobs_C,iobs_3,iobs_J  !! number of observations 
  real*8 xsig,xrwl,xrms,xswl,mele    !! sigma of process & ambiguities
  integer*4 psat                     !! corresponding satellite id
  real*8 ptime(2),ptbeg,ptend        !! living time span
  real*8 xini,xcor                   !! a priori value & correction
  real*8 map,rw,zw                   !! process noise
end type

!!
!! normal matrix
!!
type norm
  integer*4 nc,np,ns,imtx,nmtx,ipm,npm,nuk,nobs
  integer*4,pointer :: iptp(:),iptx(:)
  real*8 ltpl,sig0
  real*8,pointer :: norx(:,:)
end type
