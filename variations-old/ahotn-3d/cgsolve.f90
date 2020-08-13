SUBROUTINE cgsolve(g,q,its)

!-------------------------------------------------------------
!
!  Perform a simple conjugate gradient solution
! 
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER, INTENT(IN) :: g
INTEGER :: it
INTEGER, INTENT(OUT) :: its
REAL*8 :: alf, bet, dfmx, tmp, bnorm
REAL*8, DIMENSION(neq), INTENT(IN) :: q
REAL*8, DIMENSION(neq) :: res, sd, tmp1, phio, ones

! Initialize my guess
!phi(:,g) = 0.0
!res = q

phi(:,g) = phiold
res = q - MATMUL(jmat,phi(:,g))
sd = res
tmp = DOT_PRODUCT(res,res)
ones = 1.0

!bnorm = SQRT(tmp)
!! Want to use just the residual when the bnorm is larger than one, means more iterations
!IF (bnorm > 1.0) THEN 
!   bnorm = 1.0
!ELSE
!   ! Want to use faster multiplications at each iteration, so take inverse and store
!   bnorm = 1.0/bnorm
!END IF

DO it = 1, itmx
   ! Copy old phi vector
   phio = phi(:,g)

   ! The CG steps...should be optimized later or use outside solver
   tmp1 = MATMUL(jmat,sd)
   alf = tmp/(DOT_PRODUCT(sd,tmp1))
   phi(:,g) = phi(:,g) + alf*sd
!   dfmx = SQRT(tmp)*bnorm
   IF (MINVAL(ABS(phio)) >= tolr) THEN
      dfmx = MAXVAL(ABS(phi(:,g)/phio - ones))
   ELSE
      dfmx = MAXVAL(ABS(phi(:,g) - phio))
   END IF
   
   ! Cycle if not converged, less than itmx, and it not equal to the number of
   ! equations -- CG has zero residual when iterations = number of equations
   IF (dfmx > cgr .AND. it < itmx .AND. it /= neq) THEN
      res = res - alf*tmp1
      bet = 1.0/tmp
      tmp = DOT_PRODUCT(res,res)
      bet = tmp*bet
      sd  = res + bet*sd
      CYCLE
   ELSE IF (dfmx < cgr .OR. it == neq) THEN
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
END SUBROUTINE cgsolve
