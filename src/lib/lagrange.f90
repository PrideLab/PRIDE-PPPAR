!/* ========================================================================= */
! name	: do lagrange interpolation
! aim 	: do lagrange interpolation of data. point coordinate (x,y)
! args  : I  : real        xs : x sequence of data
!         I  : real        ys : y sequence of data
!         I  : real        ix : x coordinate that want to lagrange interpolation
!         I  : integer      n : points number of do lagrange interpolation, usually
!                               is size(xs) of size(ys)
! return: the y coordinate of ix
!/* ========================================================================= */
real*8 function lagrange(xs, ys, ix, n)
  implicit none
  real*8, intent(in)    :: xs(1:n), ys(1:n)
  real*8, intent(in)    :: ix 
  integer*4, intent(in) :: n
  !
  !! local
  real*8    :: term, res
  integer*4 :: i, j
  !
  !! initial
  term = 0.0
  res  = 0.0
  
  do i = 1, n
    term = ys(i)
    do j = 1, n
      if (j .ne. i) then
        term = term * (ix - xs(j)) / (xs(i) - xs(j))
      end if
    end do
    res = res + term
  end do

  lagrange = res
end function lagrange