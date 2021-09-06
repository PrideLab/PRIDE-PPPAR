SUBROUTINE ROT_Z(PSI, R)
!+
!
!  Rotate an r-matrix about the z-axis.
!
!  Status:  vector/matrix support routine.
!
!  Given:
!     PSI      d         angle (radians)
!
!  Given and returned:
!     R        d(3,3)    r-matrix, rotated
!
!  Sign convention:  The matrix can be used to rotate the
!  reference frame of a vector.  Calling this routine with
!  positive PSI incorporates in the matrix an additional
!  rotation, about the z-axis, anticlockwise as seen looking
!  towards the origin from positive z.
!
!  This revision:  2007 November 7
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  REAL*8 PSI, R(3, 3)
!
!! LOCAL
  INTEGER*4 I
  REAL*8 S, C, A(3, 3)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!  Matrix representing new rotation.
  S = SIN(PSI)
  C = COS(PSI)
  A = 0.D0
  DO I = 1, 3
    A(I, I) = 1.D0
  ENDDO
  A(1, 1) = C
  A(2, 1) = -S
  A(1, 2) = S
  A(2, 2) = C

!  Rotate.
  CALL MATMPY(A, R, R, 3, 3, 3)

!  Finished.
  RETURN
END
