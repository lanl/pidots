SUBROUTINE output

!-------------------------------------------------------------
!
!    Echo the flux output
!
!-------------------------------------------------------------

USE invar
USE totvar
IMPLICIT NONE
INTEGER :: i, j, k, oct, n, ieq

! Start the echo of the output for each group
! Print scalar flux as desired

! No printing if sfp is 0

! Print the full solution if sfp is 1
IF (sfp == 1) THEN
   WRITE (8,*)
   WRITE (8,112) "====================== Converged Scalar Flux Zero Moment ========================"
   DO k = 1, nzt
      DO j = 1, nyt
         WRITE (8,*) " Plane(z) : ", k, " Row(j) : ", j
         WRITE (8,113) (fluxt(i,j,k), i = 1, nxt)
      END DO
   END DO
END IF

! Print the root process' solution if sfp is 2
IF (sfp == 2) THEN
   WRITE (8,*)
   WRITE (8,112) "====================== Converged Scalar Flux Zero Moment ========================"
   DO k = 1, nz
      DO j = 1, ny
         WRITE (8,*) " Plane(z) : ", k, " Row(j) : ", j
         WRITE (8,113) (flux(i,j,k), i = 1, nx)
      END DO
   END DO
END IF

112 FORMAT(1X, A)
113 FORMAT(2X, 8ES14.6)
! 114 FORMAT(1X, A, I3, A, I3, A, I3, A, I3, A)
! 115 FORMAT(1X, A, I4, A)
! 116 FORMAT(1X, A, I3, A, 2I2, A)
! 117 FORMAT(1X, A, I3, A, I3, A)

! Report the number of warnings
WRITE (8,'(//,1X,I2,A,/)') warn, " warnings have been issued"
      
! End the input echo
WRITE (8,*)
WRITE (8,*) "-------------------------- END SOLUTION ----------------------------------"
WRITE (8,*)

RETURN
END SUBROUTINE output
