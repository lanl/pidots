SUBROUTINE idot(bit)

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
INTEGER, INTENT(IN) :: bit
INTEGER :: i, j, k, info
REAL*8, DIMENSION(neq) :: q

! Compute the solution directly
! Always use sv and product of kmat and psii to either solve for phi, or form RHS
phi = sv
DO k = 1, 8
   DO j = 1, apo
      DO i = 1, neq
         phi(i) = phi(i) + DOT_PRODUCT(kmatz(:,i,j,k),psiiz(:,j,k))
      END DO
      DO i = 1, neq
         phi(i) = phi(i) + DOT_PRODUCT(kmaty(:,i,j,k),psiiy(:,j,k))
      END DO
      DO i = 1, neq
         phi(i) = phi(i) + DOT_PRODUCT(kmatx(:,i,j,k),psiix(:,j,k))
      END DO
   END DO
END DO

! Use LAPACK solvers if matrix was factored only
CALL dgetrs('T',neq,1,jmat,neq,piv,phi,neq,info)
IF (info /= 0) THEN
   WRITE (8,'(//,1X,A)') "ERROR: dgetrs cannot solve."
   STOP
END IF

! Compute the new angular flux values at the boundaries from jpsi and kpsi
q = phi + src

! Find psio
DO k = 1, 8
   DO j = 1, apo
      DO i = 1, xys
         psioz(i,j,k) = DOT_PRODUCT(jpsiz(:,i,j,k),q)
      END DO
   END DO
END DO
DO k = 1, 8
   DO j = 1, apo
      DO i = 1, xzs
         psioy(i,j,k) = DOT_PRODUCT(jpsiy(:,i,j,k),q)
      END DO
   END DO
END DO
DO k = 1, 8
   DO j = 1, apo
      DO i = 1, yzs
         psiox(i,j,k) = DOT_PRODUCT(jpsix(:,i,j,k),q)
      END DO
   END DO
END DO

DO k = 1, 8
   DO j = 1, apo
      DO i = 1, xys
         psioz(i,j,k) = psioz(i,j,k) + DOT_PRODUCT(kpsizz(:,i,j,k),psiiz(:,j,k))
      END DO
      DO i = 1, xys
         psioz(i,j,k) = psioz(i,j,k) + DOT_PRODUCT(kpsizy(:,i,j,k),psiiy(:,j,k))
      END DO
      DO i = 1, xys
         psioz(i,j,k) = psioz(i,j,k) + DOT_PRODUCT(kpsizx(:,i,j,k),psiix(:,j,k))
      END DO
   END DO
END DO

DO k = 1, 8
   DO j = 1, apo
      DO i = 1, xzs
         psioy(i,j,k) = psioy(i,j,k) + DOT_PRODUCT(kpsiyz(:,i,j,k),psiiz(:,j,k))
      END DO
      DO i = 1, xzs
         psioy(i,j,k) = psioy(i,j,k) + DOT_PRODUCT(kpsiyy(:,i,j,k),psiiy(:,j,k))
      END DO
      DO i = 1, xzs
         psioy(i,j,k) = psioy(i,j,k) + DOT_PRODUCT(kpsiyx(:,i,j,k),psiix(:,j,k))
      END DO
   END DO
END DO

DO k = 1, 8
   DO j = 1, apo
      DO i = 1, yzs
         psiox(i,j,k) = psiox(i,j,k) + DOT_PRODUCT(kpsixz(:,i,j,k),psiiz(:,j,k))
      END DO
      DO i = 1, xzs
         psiox(i,j,k) = psiox(i,j,k) + DOT_PRODUCT(kpsixy(:,i,j,k),psiiy(:,j,k))
      END DO
      DO i = 1, xzs
         psiox(i,j,k) = psiox(i,j,k) + DOT_PRODUCT(kpsixx(:,i,j,k),psiix(:,j,k))
      END DO
   END DO
END DO

RETURN
END SUBROUTINE idot
