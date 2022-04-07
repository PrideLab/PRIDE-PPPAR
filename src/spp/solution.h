
! Struct of solution --------------------------------------------------------
type sol_t                  ! solution type 
    sequence
    type(gtime_t) time      ! time (GPST) 
    type(gtime_t) time0     ! time (GPST) of observation
    real*8 rr(6)            ! position/velocity (m|m/s) 
                            ! {x,y,z,vx,vy,vz} or {e,n,u,ve,vn,vu} 
    real*8 qr(6)            ! position variance/covariance (m^2) 
                            ! {c_xx,c_yy,c_zz,c_xy,c_yz,c_zx} or 
                            ! {c_ee,c_nn,c_uu,c_en,c_nu,c_ue} 
    real*8 dtr(6)           ! receiver clock bias to time systems (s) 
    integer*4 mytype        ! type (0:xyz-ecef,1:enu-baseline)
    integer*4 stat          ! solution status (SOLQ_???) 
    integer*4 ns            ! number of valid satellites 
    real*8 dop(4)           ! dilution of precision (gdop, pdop, hdop, vdop)
end type
