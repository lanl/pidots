SUBROUTINE cgsolve(q,its)

!-------------------------------------------------------------
!
!  Perform a simple conjugate gradient solution
! 
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER :: it
INTEGER, INTENT(OUT) :: its
REAL*8 :: alf, bet, dfmx, tmp, bnorm
REAL*8, DIMENSION(neq), INTENT(IN) :: q
REAL*8, DIMENSION(neq) :: res, sd, tmp1, phio, ones

! Initialize my guess
!phi = 0.0
!res = q

phi = phiold
res = q - MATMUL(jmat,phi)
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
   phio = phi

   ! The CG steps...should be optimized later or use outside solver
   tmp1 = MATMUL(jmat,sd)
   alf = tmp/(DOT_PRODUCT(sd,tmp1))
   phi = phi + alf*sd
!   dfmx = SQRT(tmp)*bnorm
   IF (MINVAL(ABS(phio)) >= tolr) THEN
      dfmx = MAXVAL(ABS(phi/phio - ones))
   ELSE
      dfmx = MAXVAL(ABS(phi - phio))
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
END SUBROUTINE cgsolve
