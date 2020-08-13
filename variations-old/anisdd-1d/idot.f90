SUBROUTINE idot(g,nsm,its)

!-------------------------------------------------------------
!
!  Integral discrete ordinates transport
!
!  Add sv to kmat*psii to get the new phi
!  Compute psio with jpsi*(phi+src) + kpsi*psii
!
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER, INTENT(IN) :: g, nsm
INTEGER :: i, j
INTEGER, INTENT(OUT) :: its
REAL*8, DIMENSION(nsm) :: q

! Compute the solution directly or with CG iterations
IF (idos == 0) THEN
   ! Solve for the new phi with the sum of sv with the product of kmat and psii
   ! Do according to 'tpose' flag
   IF (tpose == 1) THEN
      f(:,g) = sv
      DO j = 1, 2
         DO i = 1, nsm
            f(i,g) = f(i,g) + DOT_PRODUCT(kmat(:,i,j),psii(:,j))
         END DO
      END DO
   ELSE
      f(:,g) = sv + MATMUL(kmat(:,:,1),psii(:,1)) + MATMUL(kmat(:,:,2),psii(:,2))
   END IF
   cnvf(g) = 1
   its = 0
ELSE
   ! Solve for the new phi with CG iterations from the cgsolve routine
   IF (tpose == 1) THEN
      q = sv
      DO j = 1, 2
         DO i = 1, nsm
            q(i) = q(i) + DOT_PRODUCT(kmat(:,i,j),psii(:,j))
         END DO
      END DO
   ELSE
      q = sv + MATMUL(kmat(:,:,1),psii(:,1)) + MATMUL(kmat(:,:,2),psii(:,2))
   END IF
   ! CG only works if theres is a non-zero RHS
   IF (MAXVAL(q) <= 0.0) THEN
      f(:,g) = 0.0
   ELSE
      CALL cgsolve(g,nsm,q,its)
   END IF
END IF

! Compute the new angular flux values at the boundaries from jpsi and kpsi
q = f(:,g) + src

! Find psio depending on 'tpose' flag
IF (tpose == 1) THEN
   DO j = 1, 2
      DO i = 1, apo
         psio(i,j,g) = DOT_PRODUCT(jpsi(:,i,j),q)
         psio(i,j,g) = psio(i,j,g) + kpsi(i,j)*psii(i,j)
      END DO
   END DO
ELSE
   DO j = 1, 2
      psio(:,j,g) = MATMUL(jpsi(:,:,j),q)
      DO i = 1, apo
         psio(i,j,g) = psio(i,j,g) + kpsi(i,j)*psii(i,j)
      END DO
   END DO
END IF

RETURN
END SUBROUTINE idot
