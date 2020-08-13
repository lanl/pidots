SUBROUTINE rotmat(a,b,c,s)

!----------------------------------------------------
!
! Compute the Givens rotation matrix parameters c, s
!  for given a, b
!
!----------------------------------------------------

IMPLICIT NONE
REAL*8, INTENT(IN) :: a, b
REAL*8, INTENT(OUT) :: c, s
REAL*8 :: temp

IF (b == 0.0) THEN
   c = 1.0
   s = 0.0
ELSE IF (abs(b) > abs(a)) THEN
   temp = a/b
   s = 1.0/SQRT(1.0 + temp**2)
   c = temp*s
ELSE
   temp = b/a
   c = 1.0/SQRT(1.0 + temp**2)
   s = temp*c
END IF

RETURN
END SUBROUTINE rotmat
