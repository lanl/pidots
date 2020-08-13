SUBROUTINE cong

!-------------------------------------------------------------
!
!  Construct the gamma matrix from A^-1 * B
! 
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER, DIMENSION((ordcb+3*ordsq)) :: pv2
INTEGER :: ieq, info
REAL*8, DIMENSION((ordcb+3*ordsq)) :: wrk

ieq = ordcb + 3*ordsq

! Factor amat so it can be inverted
CALL dgetrf(ieq, ieq, amat, ieq, pv2, info)
IF (info /= 0) THEN
   WRITE (8,'(//,1X,A)') "ERROR: amat either has illegal value or is singular."
   STOP
END IF

      ! Linpack
      ! CALL dgeco(amat,ieq,ieq,pv2,rcond,wrk)
      ! IF (rcond < 1.0E-10) THEN
      !    WRITE (8,'(//,1X,A)') "WARNING: rcond very small, large condition number"
      !    warn = warn + 1
      ! END IF

! Compute the inverse of amat
CALL dgetri(ieq, amat, ieq, pv2, wrk, ieq, info)
IF (info /= 0) THEN
   WRITE (8,'(//,1X,A)') "ERROR: amat either has illegal value or is singular. Cannot invert."
   STOP
END IF

      ! Linpack
      ! info = 1
      ! CALL dgedi(amat,ieq,ieq,pv2,det,wrk,info)

! Now compute gmat with the matrix matrix multiply
gmat = MATMUL(amat,bmat)

RETURN
END SUBROUTINE cong
