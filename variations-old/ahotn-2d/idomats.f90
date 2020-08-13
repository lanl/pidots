SUBROUTINE idomats(g,neq,bcs)

!-------------------------------------------------------------
!
!  Integral discrete ordinates matrices
!
!   Control the construction of jmat, kmat, jpsi, kpsi
!   Matrices need to be made only once, then idot.f90 will
!    use them to solve the system of equations
!
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER, INTENT(IN) :: g, neq, bcs
INTEGER :: ieq, i, j, k, l, info
INTEGER, DIMENSION(:), ALLOCATABLE :: piv
REAL*8, DIMENSION(:), ALLOCATABLE :: q, wrk

! Initialize the matrices
kmat = 0.0
jpsi = 0.0
kpsi = 0.0

! Allocate the iteration Jacobian
ALLOCATE(jmat(neq,neq))
jmat = 0.0

! Allocate the temp arrays
ALLOCATE(piv(nx*ny*ordsq), q(nx*ny*ordsq), wrk(nx*ny*ordsq))

! Start the single sweep for constructing all matrices
CALL matsweep(g,bcs)

! Divide source by scattering cross section and save into vector
DO k = 0, lambda
   DO l = 0, lambda
      DO j = 1, ny
         DO i = 1, nx
            ieq = ((i-1) + (j-1)*nx)*ordsq + (l+1) + k*order
            src(ieq) = s(i,j,k,l,g)/sigs((mat(i,j)),g,g)
         END DO
      END DO
   END DO
END DO

! Multiply the source vector by the Jacobi iteration matrix
q = MATMUL(jmat,src)

! Form the LHS matrix
jmat = -jmat
DO k = 1, neq
   jmat(k,k) = 1.0 + jmat(k,k)
END DO

! Invert the new LHS matrix
CALL dgetrf(neq,neq,jmat,neq,piv,info)
IF (info /= 0) THEN
   WRITE (8,'(//,1X,A)') "ERROR: jmat either has illegal value or is singular, cannot factor."
   STOP
END IF
CALL dgetri(neq,jmat,neq,piv,wrk,neq,info)
IF (info /= 0) THEN
   WRITE (8,'(//,1X,A)') "ERROR: jmat either has illegal value or is singular, cannot invert."
   STOP
END IF

! Save the matrix vector product of the inverted matrix and the modified source
sv = MATMUL(jmat,q)

! Save the product of the inverted matrix and kmat
DO i = 1, 4
   kmat(:,:,i) = MATMUL(jmat,kmat(:,:,i))
END DO

! Now have the four needed values: sv, kmat, jpsi, kpsi

DEALLOCATE(jmat)
DEALLOCATE(piv,wrk,q)

RETURN
END SUBROUTINE idomats
