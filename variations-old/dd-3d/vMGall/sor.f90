SUBROUTINE sor(q,its)

!-------------------------------------------------------------
!
!  Perform a simple SOR solution
! 
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER :: it, i
INTEGER, INTENT(OUT) :: its
REAL*8 :: dfmx
REAL*8, DIMENSION(neq), INTENT(IN) :: q
REAL*8, DIMENSION(neq) :: phio, ones

! Initialize my guess
!phi = 0.0
phi = phiold
ones = 1.0

DO it = 1, itmx
   phio = phi
   DO i = 1, neq
      phi(i) = (q(i) - DOT_PRODUCT(jmat(1:i-1,i),phi(1:i-1)) - DOT_PRODUCT(jmat(i+1:neq,i),phi(i+1:neq)))/jmat(i,i)
      phi(i) = sorw*phi(i) + (1.0 - sorw)*phio(i)
   END DO

   IF (MINVAL(ABS(phio)) >= tolr) THEN
      dfmx = MAXVAL(ABS(phi/phio - ones))
   ELSE
      dfmx = MAXVAL(ABS(phi - phio))
   END IF
   
   IF (dfmx > cgr .AND. it < itmx) THEN
      CYCLE
   ELSE IF (dfmx < cgr) THEN
!      cnvf = 1
      its = it
      EXIT
   ELSE IF (it == itmx) THEN
!      cnvf = 0
      its = it
      EXIT
   END IF
END DO

RETURN
END SUBROUTINE sor
