SUBROUTINE ROT_X(PHI, R)
!+
!
!  Rotate an r-matrix about the x-axis.
!
!  Status:  vector/matrix support routine.
!
!  Given:
!     PHI      d         angle (radians)
!
!  Given and returned:
!     R        d(3,3)    r-matrix
!
!  Sign convention:  The matrix can be used to rotate the
!  reference frame of a vector.  Calling this routine with
!  positive PHI incorporates in the matrix an additional
!  rotation, about the x-axis, anticlockwise as seen looking
!  towards the origin from positive x.
!
!  This revision:  2007 November 7
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  REAL*8 PHI, R(3, 3)
!
!! LOCAL
  INTEGER*4 I
  REAL*8 S, C, A(3, 3)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!  Matrix representing new rotation.
  S = DSIN(PHI)
  C = DCOS(PHI)
  A = 0.D0
  DO I = 1, 3
    A(I, I) = 1.D0
  ENDDO
  A(2, 2) = C
  A(3, 2) = -S
  A(2, 3) = S
  A(3, 3) = C

!  Rotate.
  CALL MATMPY(A, R, R, 3, 3, 3)

!  Finished.
  RETURN
END
