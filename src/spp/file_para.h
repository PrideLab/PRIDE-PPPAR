
! Included files ------------------------------------------------------------
include 'const.h'
include 'comtyp.h'
include 'option.h'
include 'rnxnav.h'
include 'rnxobs.h'
include 'sigind_t.h'
include 'solution.h'
include 'station.h'
include 'rtktyp.h'

! Global variables ----------------------------------------------------------
real*8 :: chisqr(100)               ! chi-sqr(n) table (alpha=0.001)
real*8 :: ura_value(15)             ! ura max values
real*8 :: timeoffset_               ! time offset (s)
type(prcopt_t) :: prcopt_default    ! defaults processing options
type(solopt_t) :: solopt_default    ! defaults solution output options

type(sol_t), pointer :: allsol_(:)  ! all solutions
integer*4 :: solindex_              ! solution index
integer*4 :: headwritten_           ! output the header (0-none, 1-written)
character(1024) :: clkfile_

! postpos.f90 ---------------------------------------------------------------
type(obs_t) :: obss                 ! observation data 
type(nav_t) :: navs                 ! navigation data 
type(sta_t) :: stas(MAXRCV)         ! station infomation 
integer*4 :: nepoch                 ! number of observation epochs 
integer*4 :: iobsu                  ! current rover observation data index 
integer*4 :: iobsr                  ! current reference observation data index 
integer*4 :: revs                   ! analysis direction (0:forward,1:backward) 

common chisqr, ura_value, timeoffset_, prcopt_default, solopt_default, &
       allsol_, solindex_, headwritten_, clkfile_, obss, navs, stas, &
       nepoch, iobsu, iobsr, revs
