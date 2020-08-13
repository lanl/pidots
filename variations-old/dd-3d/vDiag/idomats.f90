SUBROUTINE idomats(g)

!-------------------------------------------------------------
!
!  Integral discrete ordinates matrices
!
!   Control the construction of jmat, kmat, jpsi, kpsi
!   Matrices nned to be made only once, then idot.f90 will
!    use them to solve the system of equations
!
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER, INTENT(IN) :: g
INTEGER :: i, j, k

! Allocate the iteration Jacobian
ALLOCATE(jmat(7,nx,ny,nz))
jmat = 0.0

!-------------------------------------------------------------------
! Can now start the single sweep for constructing all matrices
CALL matsweep(g)

! Begin constructing the RHS
DO k = 1, nz
   DO j = 1, ny
      DO i = 1, nx
         src(i,j,k) = s(i,j,k,g)/sigs((mat(i,j,k)),g,g)
      END DO
   END DO
END DO

! Multiply the source vector by the Jacobi iteration matrix
DO k = 1, nz
   DO j = 1, ny
      DO i = 1, nx
         sv(i,j,k) = jmat(4,i,j,k)*src(i,j,k)
         IF (i /= 1) sv(i,j,k) = sv(i,j,k) + jmat(1,i,j,k)*src(i-1,j,k)
         IF (j /= 1) sv(i,j,k) = sv(i,j,k) + jmat(2,i,j,k)*src(i,j-1,k)
         IF (k /= 1) sv(i,j,k) = sv(i,j,k) + jmat(3,i,j,k)*src(i,j,k-1)
         IF (i /= nx) sv(i,j,k) = sv(i,j,k) + jmat(5,i,j,k)*src(i+1,j,k)
         IF (j /= ny) sv(i,j,k) = sv(i,j,k) + jmat(6,i,j,k)*src(i,j+1,k)
         IF (k /= nz) sv(i,j,k) = sv(i,j,k) + jmat(7,i,j,k)*src(i,j,k+1)
      END DO
   END DO
END DO

! Form the LHS matrix
jmat = -jmat
DO k = 1, nz
   DO j = 1, ny
      DO i = 1, nx
         jmat(4,i,j,k) = 1.0 + jmat(4,i,j,k)
      END DO
   END DO
END DO

RETURN
END SUBROUTINE idomats
