
! Struct of station ---------------------------------------------------------
type sta_t                   ! station parameter type 
    sequence
    character(MAXANT) name   ! marker name 
    character(MAXANT) marker ! marker number 
    character(MAXANT) antdes ! antenna descriptor 
    character(MAXANT) antsno ! antenna serial number 
    character(MAXANT) rectype ! receiver type descriptor 
    character(MAXANT) recver ! receiver firmware version 
    character(MAXANT) recsno ! receiver serial number 
    integer*4 antsetup       ! antenna setup id 
    integer*4 itrf           ! ITRF realization year 
    integer*4 deltype        ! antenna delta type (0:enu,1:xyz) 
    real*8 pos(3)            ! station position (ecef) (m) 
    real*8 del(3)            ! antenna position delta (e/n/u or x/y/z) (m) 
    real*8 hgt               ! antenna height (m) 
end type
