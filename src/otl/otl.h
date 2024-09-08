!!
!! Ocean Thermal Loading Setting Definition
!!
type otls
  !! Calculate Settings
  integer*4  mjd0,mjd1    !! Modified Julian day start/end
  real*8  sod0,sod1       !! second of day start/end
  real*8  plon,plat,phgt  !! Position longitude/latitude/height
  !! Output Settings
  real*8  interval        !! interval of calculate ocean thermal loading
  character*256 otlname   !! Ocean Thermal Loading file name
end type