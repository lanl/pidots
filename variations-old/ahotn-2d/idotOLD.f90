SUBROUTINE idot(g)

!-------------------------------------------------------------
!
!  Integral discrete ordinates transport
!
!   Allocate matrices for solution, call for the differential
!   sweep. When matrices are made and return, construct RHS,
!   symmetrize, call for solution
!
!-------------------------------------------------------------

USE invar
USE solvar
USE timevar
IMPLICIT NONE
INTEGER, INTENT(IN) :: g
INTEGER :: neq, bcs, bcs2, ieq, i, j, k, l, m, ndx, jj, info
REAL*8 :: symm, symm2
REAL*8, DIMENSION(ordsq*nx*ny) :: src, q, sol

! Allocate the Iteration Jacobian
neq = ordsq*nx*ny
ALLOCATE(jmat(neq,neq))
jmat = 0.0

! Allocate the BC coefficient matrix
bcs = (nx+ny)*apo*order
ALLOCATE(kmat(neq,bcs,4))
kmat = 0.0

! Allocate the angular flux out matrices
bcs2 = (nx+ny)*order
ALLOCATE(jpsi(bcs2,neq,apo,4))
ALLOCATE(kpsi(bcs2,bcs2,apo,4))
jpsi = 0.0
kpsi = 0.0

! Can now start the single sweep for constructing all matrices
CALL matsweep(g,bcs)

! Begin constructing the RHS
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

! Add the contributions from the four quadrants' angular flux in
q = q + MATMUL(kmat(:,:,1),psi1) + MATMUL(kmat(:,:,2),psi2) + MATMUL(kmat(:,:,3),psi3) + MATMUL(kmat(:,:,4),psi4)

! Form the LHS matrix
jmat = -jmat
DO k = 1, neq
   jmat(k,k) = 1.0 + jmat(k,k)
END DO

! Symmetrize the matrix
DO j = 1, ny
   DO i = 1, nx
      m = mat(i,j)
      symm = sigs(m,g,g)*dx(i)*dy(j)      
      DO k = 0, lambda
         DO l = 0, lambda
            IF (order == 1) symm2 = symm
            IF (order /= 1) symm2 = symm*(2.0*(k+1)-1.0)*(2.0*(l+1)-1.0)
            ndx = ((j-1)*nx + (i-1))*ordsq + (l+1) + k*order
            q(ndx) = symm2*q(ndx)
            DO jj = 1, neq
               jmat(ndx,jj) = symm2*jmat(ndx,jj)
            END DO
         END DO
      END DO
   END DO
END DO

! Determine if the matrix will be put into a file. If so, do it.
READ(7,*) matrix
IF (matrix == 1 .AND. g == 1) THEN
   OPEN (UNIT = 14, FILE = "jmat")
   WRITE(14,*) "! File containing the matrix, RHS, and solution -- GROUP 1 ONLY"
   WRITE(14,*) "! Problem scope: I J lambda, Max #its, Conv. criterion"
   WRITE(14,'(2X,I5,1X,I5,1X,I2,1X,I5,1X,ES9.3)') nx, ny, lambda, itmx, err
   WRITE(14,*)
   WRITE(14,*) "! JMAT-symmetric printed row-wise (i.e., j varies faster than i)"
   DO i = 1, neq
      WRITE(14,116) (jmat(i,j), j = 1, neq)
   END DO
   WRITE(14,*)
   WRITE(14,*) "! RHS vector printed first to last"
   WRITE(14,116) (q(i), i = 1, neq)
   WRITE(14,*)
   WRITE(14,*) "! Solution printed first to last in vector"
END IF
116 FORMAT(2X,8ES14.6)

! Get the time to construct the matrix
! Time will be higher if required to copy matrix to file
CALL CPU_TIME(tjmat)

READ(7,*) itmflag
IF (itmflag == 1) THEN
   ! Try a conjugate gradient solver
   CALL cgsolve(g,q,sol)
   ! Properly place the solution in the f-matrix
   DO k = 0, lambda
      DO l = 0, lambda
         DO j = 1, ny
            DO i = 1, nx
               ieq = ((i-1) + (j-1)*nx)*ordsq + (l+1) + k*order
               f(i,j,k,l,g) = sol(ieq)
            END DO
         END DO
      END DO
   END DO

   ! Write the solution to separate file if requested
   IF (matrix == 1 .AND. g == 1) WRITE(14,116) (sol(i), i = 1, neq)

ELSE
   IF (itmflag /= 0) THEN
      WRITE (8,'(/,1X,A)') "WARNING: Flag for ITM solution type not 0 or 1, direct solver used by default."
      warn = warn + 1
   END IF
   ! Solve the matrix system for the converged scalar flux soluton
   ! Asymmetric solver
     ! CALL dgesv(neq,1,jmat,neq,ipiv,q,neq,info)
   ! Symmetric solver
   CALL dposv('U',neq,1,jmat,neq,q,neq,info)
   IF (info /= 0) THEN
      WRITE (8, '(//,1X,A)') "ERROR: matrix either has illegal value or is singular."
      STOP
   END IF
   cnvf(g) = 1
   ! Properly place the solution in the f-matrix
   DO k = 0, lambda
      DO l = 0, lambda
         DO j = 1, ny
            DO i = 1, nx
               ieq = ((i-1) + (j-1)*nx)*ordsq + (l+1) + k*order
               f(i,j,k,l,g) = q(ieq)
            END DO
         END DO
      END DO
   END DO

   ! Write the solution to separate file if requested
   IF (matrix == 1) WRITE(14,116) (q(i), i = 1, neq)

END IF

! Compute the new angular flux values at the boundaries from jpsi and kpsi
IF (itmflag == 1) THEN
   q = sol + src
ELSE
   q = q + src
END IF
DO j = 1, 4
   DO i = 1, apo
      psio(:,i,j,g) = MATMUL(jpsi(:,:,i,j),q)
      k = (i-1)*(nx+ny)*order
      l = i*(nx+ny)*order
      IF (j == 1) THEN
         psio(:,i,j,g) = psio(:,i,j,g) + MATMUL(kpsi(:,:,i,j),psi1((k+1):l))
      ELSE IF (j == 2) THEN
         psio(:,i,j,g) = psio(:,i,j,g) + MATMUL(kpsi(:,:,i,j),psi2((k+1):l))
      ELSE IF (j == 3) THEN
         psio(:,i,j,g) = psio(:,i,j,g) + MATMUL(kpsi(:,:,i,j),psi3((k+1):l))
      ELSE
         psio(:,i,j,g) = psio(:,i,j,g) + MATMUL(kpsi(:,:,i,j),psi4((k+1):l))
      END IF
   END DO
END DO

DEALLOCATE(jmat,kmat,jpsi,kpsi)

RETURN
END SUBROUTINE idot
