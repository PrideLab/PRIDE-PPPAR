!/* ========================================================================= */
! name  : do lagrange interpolation
! aim   : do lagrange interpolation for ocean load corrections
! args  : I  : type(otls) OTL : struct of OTL, control used
!         I  : integer    fjd : find modified Julian Day
!         I  : real      fsod : find socond of mjd 
!         I  : integer   dim  : dimension of lagrange interpolation
!         I  : integer lfnotl : the lfn of 30 minute interval otl
!         I  : real otlfil_INT: the interval of output oceanload 
!         O  : real       enu : the corrections of enu
! return: the specific interpolation result you want get
! Contributor: Yingda Deng
!/* ========================================================================= */
subroutine lag_interp_otl(OTL, fjd, fsod, dim, lfnotl, otlfil_INT, enu)
  implicit none
  include 'otl.h'
  type(otls) :: OTL
  integer*4, intent(in) :: fjd, dim
  real*8, intent(in)    :: fsod
  integer*4, intent(in) :: lfnotl
  real*8, intent(in)    :: otlfil_INT
  real*8, intent(out)   :: enu(3)
  !
  !! local
  integer*4 :: jdt, stat, i, j
  real*8 :: sodt, enut(3), x_intp(1:dim), y_intp(1:dim, 1:3)
  ! function called
  real*8, external :: timdif
  real*8, external :: lagrange

  x_intp=0.d0;y_intp=0.d0;enu=0.d0
  !
  !! read result
  if (timdif(fjd,fsod,OTL%mjd0,OTL%sod0) .lt. 0.d0 .and. timdif(fjd,fsod,OTL%mjd1,OTL%sod1) .gt. 0.d0) then
    write(*,'(a)') "***ERROR(lagrange_interp_oceanload): Outside the session"
  endif
  !
  !! find interpolate epoch
  do
    read(lfnotl,*,end=902,iostat=stat) jdt,sodt,(enut(i), i=1,3)
    if(stat .ne. 0) exit
    if (timdif(fjd,fsod,jdt,sodt) .le. dim/2*otlfil_INT) then
      exit
    endif
  enddo

  x_intp(1)=(jdt-fjd)*864.d2+sodt-fsod;
  y_intp(1, 1:3)=enut(1:3)
  do i=2, dim
    read(lfnotl,*,end=902,iostat=stat) jdt,sodt,(y_intp(i,j), j=1,3)
    if(stat .ne. 0) exit
    x_intp(i)=(jdt-fjd)*864.d2+sodt-fsod;
  enddo
  !
  !! do lagrange interpolation
  do i=1, 3
    enu(i) = lagrange(x_intp(1:dim), y_intp(1:dim, i), 0.d0, dim)
  enddo

902 continue
  return
end subroutine