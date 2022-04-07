
! Struct of observation -----------------------------------------------------
type obsd_t                  ! observation data record 
    sequence
    type(gtime_t) time       ! receiver sampling time (GPST) 
    integer*4 sat,rcv        ! satellite/receiver number 
    integer*4 SNR  (NFREQ+NEXOBS) ! signal strength (0.25 dBHz) 
    integer*4 LLI  (NFREQ+NEXOBS) ! loss of lock indicator 
    integer*4 code (NFREQ+NEXOBS) ! code indicator (CODE_???) 
    integer*4 qualL(NFREQ+NEXOBS) ! quality of carrier phase measurement 
    integer*4 qualP(NFREQ+NEXOBS) ! quality of pseudorange measurement 
    real*8 L(NFREQ+NEXOBS)   ! observation data carrier-phase (cycle) 
    real*8 P(NFREQ+NEXOBS)   ! observation data pseudorange (m) 
    real*8 D(NFREQ+NEXOBS)   ! observation data doppler frequency (Hz) 
end type

type obs_t                   ! observation data 
    sequence
    integer*4 n,nmax         ! number of obervation data/allocated 
    real*8 tint              ! time interval
    type(obsd_t), pointer :: mydata(:) ! observation data records 
end type
