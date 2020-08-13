SUBROUTINE cgsolvesor(q,its)

!-------------------------------------------------------------
!
!  Perform a preconditioned conjugate gradient solution
!  M = (1/(2-w))(D/w+E)*(D/w)^-1*(D/w+F) -- Symmetric SOR preconditioner
!
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER :: it, i, j
INTEGER, INTENT(OUT) :: its
REAL*8 :: alf, bet, dfmx, tmp, tmp2, bnorm
REAL*8, DIMENSION(neq), INTENT(IN) :: q
REAL*8, DIMENSION(neq) :: res, rt, x, z, sd, tmp1, phio, ones

! Initialize my guess
!phi = 0.0
!res = q

phi = phiold
res = q - MATMUL(jmat,phi)

! Store modified sorw factor
tmp2 = 2.0 - sorw

! Solve the pre-conditioner problem Mz = r.
rt = res
DO j = 1, neq               ! Solve the lower triangular Lx=r
   x(j) = rt(j)*dv(j)
   DO i = j+1, neq
      rt(i) = rt(i) - jmat(i,j)*x(j)
   END DO
END DO
! reset x with the SOR factor
x = tmp2*x/dv
DO j = neq, 1, -1           ! Solve the upper triangular Uz=x
   z(j) = x(j)*dv(j)
   DO i = 1, j-1
      x(i) = x(i) - jmat(i,j)*z(j)
   END DO
END DO

sd = z
tmp = DOT_PRODUCT(res,z)
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

      ! Solve the pre-conditioned system Mz = r.
      rt = res
      DO j = 1, neq
         x(j) = rt(j)*dv(j)
         DO i = j+1, neq
            rt(i) = rt(i) - jmat(i,j)*x(j)
         END DO
      END DO
      x = tmp2*x/dv
      DO j = neq, 1, -1
         z(j) = x(j)*dv(j)
         DO i = 1, j-1
            x(i) = x(i) - jmat(i,j)*z(j)
         END DO
      END DO

      bet = 1.0/tmp
      tmp = DOT_PRODUCT(res,z)
      bet = tmp*bet
      sd  = z + bet*sd
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
END SUBROUTINE cgsolvesor
