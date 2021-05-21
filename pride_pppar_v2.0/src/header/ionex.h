! ionex
! Single-shell model
type ionex
  logical*1   lrot
  character*4 cmf,ctm
  real*8      bradius,height,dintv
  real*8      lat1,lat2,dlat
  real*8      lon1,lon2,dlon
  real*8      sdcb(maxsat)
!
!! ionosphere map data
  integer*4 blon,blat,bmap
  integer*4,pointer :: mjd(:)
  real*8,pointer    :: sod(:)
  real*8,pointer :: tecmp(:,:,:)
  real*8,pointer :: rmsmp(:,:,:)
end type
