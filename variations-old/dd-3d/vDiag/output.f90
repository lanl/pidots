SUBROUTINE output

!-------------------------------------------------------------
!
!    Echo the flux output
!
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER :: i, j, k, g

DO g = 1, ng
   WRITE (8,*)
   WRITE (8,112) "========== Group ", g, " Converged Scalar Flux Zero Moment =========="
   DO k = 1, nz
      DO j = 1, ny
         WRITE (8,*) " Plane(z) : ", k, " Row(j) : ", j
         WRITE (8,113) (phi(i,j,k,g), i = 1, nx)
      END DO
   END DO
END DO

OPEN (UNIT=20, FILE="flux", STATUS="REPLACE")
g = 1
DO k = 1, nz
   DO j = 1, ny
      DO i = 1, nx
         WRITE(20,'(ES14.6)') phi(i,j,k,g)
      END DO
   END DO
END DO

112 FORMAT(1X, A, I4, A)
113 FORMAT(2X, 8ES14.6)

! Report the number of warnings
WRITE (8,'(//,1X,I2,A,/)') warn, " warnings have been issued"
      
! End the input echo
WRITE (8,*)
WRITE (8,*) "-------------------------- END SOLUTION ----------------------------------"
WRITE (8,*)

RETURN
END SUBROUTINE output
