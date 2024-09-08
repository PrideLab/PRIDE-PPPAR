!/* ========================================================================= */
! name	: do lagrange interpolation of sp3
! aim 	: do lagrange interpolation in a given time, based on blocks data
! args  : I  : integer   fjd  : find modified Julian day
!         I  : real      fsod : find socond of mjd 
!         I  : character fprn : find satellite prn 
!         I  : orbhdr   header: the header information of sp3
!         I  : sp3block blocks: all the block information of sp3
!         I  : integer  nepoch: the number of epoch of sp3
!         I  : integer  dim   : dimension of lagrange interpolation ,should be even
!         O  : logical  iflag : the flag of interpolation
!         O  : real     fdata : the data that have benn interpolated
! return: the specific interpolation result you want get
!/* ========================================================================= */
subroutine lagrange_interp_sp3(fjd, fsod, fprn, header, blocks, nepoch, dim, iflag, fdata)
  implicit none
  include '../header/const.h'
  include '../header/orbit.h'
  integer*4, intent(in)   :: fjd
  real*8, intent(in)      :: fsod
  character*3, intent(in) :: fprn
  type(orbhdr), intent(in)   :: header
  type(sp3block), intent(in) :: blocks(1:nepoch)
  integer*4, intent(in)   :: nepoch, dim
  real*8, intent(out)     :: fdata(1:3)
  logical*1, intent(out)  :: iflag
!
!! local
  integer*4,parameter :: MAXSESSION = 3600 ! 60 minute
  integer*4 :: i, j
  integer*4 :: left, right, pos  ! back epoch number, forward epoch number, postion of fprn
  integer*4 :: ln, ln0, rn, count
  real*8    :: x_intp(1:dim), xi
  real*8    :: y_intp(1:dim, 1:3)
  logical*1 :: has_prn
!
!! function call
  real*8, external :: timdif
  real*8, external :: lagrange
  
  if (mod(dim, 2) .ne. 0) then
    write(*, *) "***WARNING(lagrange_interp_sp3): The dimension of interpolation should be even."
    return
  endif
!
!! assigning variables
  iflag = .false.
  left  = 0
  right = 0
  ln    = dim / 2
  rn    = dim / 2
  fdata(1:3) = 1.d15
!
!! find fprn position in all prns
  has_prn = .false.
  do pos = 1, header%nprn
    if(fprn .eq. header%prn(pos)) then
      has_prn = .true.
      exit
    endif
  enddo
  if(has_prn .eqv. .false.)then
    write(*, *) "***WARNING(lagrange_interp_sp3): don't have prn in prn list", fprn
    return
  endif
!
!! find total back and total forward epoch number during the MAXSESSION time
  do i = 1, nepoch
    if(((abs(timdif(blocks(i)%jd, blocks(i)%sod, fjd, fsod))) .gt. MAXSESSION) .or. (blocks(i)%flag(pos) .eqv. .false.)) then
      cycle
    endif
    if(timdif(blocks(i)%jd, blocks(i)%sod, fjd, fsod) .lt. -1E-5) then
      left = left + 1
    else if(timdif(blocks(i)%jd, blocks(i)%sod, fjd, fsod) .gt. 1E-5) then
      right = right + 1
    endif
  enddo
!
!! adjust left and right interpolation record number
  if ((left + right) < ln + rn .or. left < 1 .or. right < 1) then
    write(*, '(2a)') "***WARNING(lagrange_interp_sp3): don't have enough epoch of ", fprn
    return
  else if (ln .gt. left .and. rn .lt. right) then
    ln = left
    rn = dim - ln
  else if (ln .lt. left .and. rn .gt. right) then
    rn = right
    ln = dim - rn
  endif
!
!! find fsod time ( i=epoch number of fsod )
  do i = 1, nepoch
    if(abs(timdif(blocks(i)%jd, blocks(i)%sod, fjd, fsod)) .lt. 1E-5) exit
  enddo
  xi = real(i)
!
!! counting use assigning
  count = ln
  ln0 = ln
  j = 0
!
!! find closest epoch from left epoch
  do while(j .ne. ln0)
    j = j + 1
    if(blocks(i-j)%flag(pos) .eqv. .false.) then
      ln0 = ln0 + 1
      cycle
    endif
    x_intp(count) = (timdif(blocks(i-j)%jd, blocks(i-j)%sod, header%jd0, header%sod0)) / header%dintv + 1
    y_intp(count, 1:3) = blocks(i-j)%x(1:3, pos)
    count = count - 1
  enddo
!
!! find closest epoch from right epoch
  count = 1
  j = 0
  do while(j .ne. rn)
    j = j + 1
    if(blocks(i + j)%flag(pos) .eqv. .false.) then
      rn = rn + 1
      cycle
    endif
    x_intp(ln + count) = (timdif(blocks(i+j)%jd, blocks(i+j)%sod, header%jd0, header%sod0)) / header%dintv + 1
    y_intp(ln + count, 1:3) = blocks(i+j)%x(1:3, pos)
    count = count + 1
  enddo
!
!! do lagrange interpolation
  do i = 1, 3
    fdata(i) = lagrange(x_intp(1:dim), y_intp(1:dim, i), xi, dim)
  enddo
  iflag = .true.
  return
end subroutine