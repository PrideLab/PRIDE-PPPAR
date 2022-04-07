
! Signal index type ---------------------------------------------------------
type sigind_t                                             ! signal index type
    sequence
    integer*4 n                                           ! number of index 
    integer*4 frq(MAXOBSTYPE)                             ! signal frequency (1:L1,2:L2,...) 
    integer*4 pos(MAXOBSTYPE)                             ! signal index in obs data (-1:no) 
    integer*4 pri(MAXOBSTYPE)                             ! signal priority (15-0) 
    integer*4 mytype(MAXOBSTYPE)                          ! type (0:C,1:L,2:D,3:S) 
    integer*4 code(MAXOBSTYPE)                            ! obs code (CODE_L??) 
    real*8 shift(MAXOBSTYPE)                              ! phase shift (cycle) 
end type
