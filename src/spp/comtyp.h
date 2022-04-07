
! Common type ---------------------------------------------------------------
type gtime_t                 ! time struct 
    sequence
    integer*4 time           ! time (s) expressed by standard time_t
    real*8 sec               ! fraction of second under 1 s 
end type

type ssat_t                  ! satellite status type 
    sequence
    integer*4 sys            ! navigation system 
    integer*4 vs             ! valid satellite flag single 
    real*8 azel(2)           ! azimuth/elevation angles {az,el} (rad) 
    real*8 resp(NFREQ)       ! residuals of pseudorange (m) 
    real*8 resc(NFREQ)       ! residuals of carrier-phase (m) 
    real*8 icbias(NFREQ)     ! glonass IC bias (cycles) 
    integer*4 vsat(NFREQ)    ! valid satellite flag 
    integer*4 snr (NFREQ)    ! signal strength (0.25 dBHz) 
    integer*4 lock (NFREQ)   ! lock counter of phase 
end type
