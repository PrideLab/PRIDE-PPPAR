!
!! MATMPY.f90
!!
!!    Copyright (C) 2021 by Wuhan University
!!
!!    This program is an open source software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License (version 3) as
!!    published by the Free Software Foundation.
!!
!!    This program is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License (version 3) for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with this program.  If not, see <https://www.gnu.org/licenses/>.
!!
!! Contributor: Maorong Ge
!! 
!!
!
SUBROUTINE MATMPY(A, B, C, NROW, NCOLA, NCOLB)
!
!     PREMULTIPLY MATRIX B BY MATRIX A WITH RESULTS IN MATRIX C
!
!     PARAMETERS
!          A       I  INPUT MATRIX
!          B       I  INPUT MATRIX
!          C       O  OUTPUT MATRIX
!          NROW    I  NUMBER OF ROWS OF MATRIX A
!          NCOLA   I  NUMBER OF COLUMNS OF MATRIX A
!          NCOLB   I  NUMBER OF COLUMNS OF MATRIX B
!
  implicit none
  INTEGER*4 NROW, NCOLA, NCOLB
  REAL*8 A(NROW, NCOLA), B(NCOLA, NCOLB), C(NROW, NCOLB)
!
!! local
  INTEGER*4 I, J, K
  REAL*8 SUM
  REAL*8, POINTER :: C1(:, :)

  ALLOCATE (C1(1:NROW, 1:NCOLB))
!
  DO I = 1, NROW
    DO J = 1, NCOLB
      SUM = 0.D0
      DO k = 1, NCOLA
        SUM = SUM + A(I, k)*B(k, J)
      enddo
      C1(I, J) = SUM
    enddo
  enddo

  DO I = 1, NROW
    DO J = 1, NCOLB
      C(I, J) = C1(I, J)
    ENDDO
  ENDDO
!
  DEALLOCATE (C1)
  RETURN
END
