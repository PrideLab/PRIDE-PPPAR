
! Struct of rtk control -----------------------------------------------------
type rtk_t                    ! RTK control/result type 
    sequence
    type(sol_t) sol           ! RTK solution 
    real*8 rb(6)              ! base position/velocity (ecef) (m|m/s) 
    integer*4 nx,na           ! number of float states/fixed states 
    real*8 tt                 ! time difference between current and previous (s) 
    real*8, pointer :: x(:), P(:,:)  ! float states and their covariance 
    real*8, pointer :: xa(:),Pa(:,:) ! fixed states and their covariance 
    integer*4 nfix            ! number of continuous fixes of ambiguity 
    integer*4 excsat          ! index of next satellite to be excluded for partial ambiguity resolution 
    real*8 com_bias           ! phase bias common between all sats (used to be distributed to all sats 
    type(ssat_t) ssat(MAXSAT) ! satellite status 
    integer*4 neb             ! bytes in error message buffer 
    type(prcopt_t) opt        ! processing options 
    integer*4 initial_mode    ! initial positioning mode 
end type
