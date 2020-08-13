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
INTEGER :: it, neq, i
INTEGER, INTENT(OUT) :: its
REAL*8 :: alf, bet, df, dfmx, tmp, bnorm
REAL*8, DIMENSION((nmom*nx*ny)), INTENT(IN) :: q
REAL*8, DIMENSION((nmom*nx*ny)) :: res, sd, tmp1, phio

neq = nmom*nx*ny

! Initialize my guess
phi(:,g) = 0.0
res = q
sd = res
tmp = DOT_PRODUCT(res,res)

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
   res = res - alf*tmp1
   bet = 1.0/tmp
   tmp = DOT_PRODUCT(res,res)
   bet = tmp*bet
   sd  = res + bet*sd
!   dfmx = SQRT(tmp)*bnorm
   dfmx = -1.0
   DO i = 1, neq, nmom
      IF (ABS(phio(i)) >= tolr) THEN
         df = ABS(phi(i,g)/phio(i) - 1.0)
      ELSE
         df = ABS(phi(i,g) - phio(i))
      END IF
      IF (df > dfmx) dfmx = df
   END DO

   IF (dfmx > cgr .AND. it < itmx .AND. it /= neq) THEN
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
