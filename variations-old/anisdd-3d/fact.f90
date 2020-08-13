FUNCTION fact(n)

!-------------------------------------------------------------
!
! Compute the factorial, n!
!
!-------------------------------------------------------------

IMPLICIT NONE
INTEGER, INTENT(IN) :: n
INTEGER :: fact, i

fact = 1
IF (n > 0) THEN
   DO i = 1, n
      fact = fact*i
   END DO
END IF

END FUNCTION fact
