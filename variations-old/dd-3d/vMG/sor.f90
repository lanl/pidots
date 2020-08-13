SUBROUTINE sor(g,q,its)

!-------------------------------------------------------------
!
!  Perform a simple SOR solution
! 
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER, INTENT(IN) :: g
INTEGER :: it, i
INTEGER, INTENT(OUT) :: its
REAL*8 :: dfmx
REAL*8, DIMENSION(neq), INTENT(IN) :: q
REAL*8, DIMENSION(neq) :: phio, ones

! Initialize my guess
!phi(:,g) = 0.0
phi(:,g) = phiold
ones = 1.0

DO it = 1, itmx
   phio = phi(:,g)
   DO i = 1, neq
      phi(i,g) = (q(i) - DOT_PRODUCT(jmat(1:i-1,i),phi(1:i-1,g)) - DOT_PRODUCT(jmat(i+1:neq,i),phi(i+1:neq,g)))/jmat(i,i)
      phi(i,g) = sorw*phi(i,g) + (1.0 - sorw)*phio(i)
   END DO

   IF (MINVAL(ABS(phio)) >= tolr) THEN
      dfmx = MAXVAL(ABS(phi(:,g)/phio - ones))
   ELSE
      dfmx = MAXVAL(ABS(phi(:,g) - phio))
   END IF
   
   IF (dfmx > cgr .AND. it < itmx) THEN
      CYCLE
   ELSE IF (dfmx < cgr) THEN
      cnvf(g) = 1
      its = it
      EXIT
   ELSE IF (it == itmx) THEN
      cnvf(g) = 0
      its = it
      EXIT
   END IF
END DO

RETURN
END SUBROUTINE sor
