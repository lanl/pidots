SUBROUTINE idot(g,its)

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
INTEGER, INTENT(IN) :: g
INTEGER :: i, j, k, it
INTEGER, INTENT(OUT) :: its
REAL*8 :: df, dfmx, offdiag
REAL*8, DIMENSION(nx,ny,nz) :: phio

phi(:,:,:,g) = 0.0

DO it = 1, itmx
   phio = phi(:,:,:,g)
   DO k = 1, nz
      DO j = 1, ny
         DO i = 1, nx
            offdiag = 0.0
            IF (i /= 1) offdiag = offdiag + jmat(1,i,j,k)*phi(i-1,j,k,g)
            IF (j /= 1) offdiag = offdiag + jmat(2,i,j,k)*phi(i,j-1,k,g)
            IF (k /= 1) offdiag = offdiag + jmat(3,i,j,k)*phi(i,j,k-1,g)
            IF (i /= nx) offdiag = offdiag + jmat(5,i,j,k)*phi(i+1,j,k,g)
            IF (j /= ny) offdiag = offdiag + jmat(6,i,j,k)*phi(i,j+1,k,g)
            IF (k /= nz) offdiag = offdiag + jmat(7,i,j,k)*phi(i,j,k+1,g)
            phi(i,j,k,g) = (sv(i,j,k) - offdiag)/jmat(4,i,j,k)
         END DO
      END DO
   END DO

   dfmx = -1.0
   DO k = 1, nz
      DO j = 1, ny
         DO i = 1, nx
            IF (ABS(phio(i,j,k)) >= tolr) THEN
               df = ABS((phi(i,j,k,g) - phio(i,j,k))/phio(i,j,k))
            ELSE
               df = ABS(phi(i,j,k,g) - phio(i,j,k))
            END IF
            ! Find the largest value
            IF (df > dfmx) THEN
               dfmx = df
            END IF
         END DO
      END DO
   END DO

   IF (dfmx > err .AND. it < itmx) THEN
      CYCLE
   ELSE IF (dfmx < err) THEN
      cnvf(g) = 1
      its = it
      WRITE(8,*)
      WRITE(8,*) "Converged in ", its, " iterations."
      WRITE(8,*)
      EXIT
   ELSE IF (it == itmx) THEN
      cnvf(g) = 0
      its = it
      EXIT
   END IF

END DO

RETURN
END SUBROUTINE idot
