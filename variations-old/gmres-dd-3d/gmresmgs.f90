SUBROUTINE gmresmgs(g)

!-----------------------------------------------------
!
! Solve the parallel system of equations with
!  a GMRES(m) iterative scheme (may have restarts)
!
! Use Modified Gram Schmidt orthogonalization
!  for the Arnoldi method
!
!-----------------------------------------------------

USE invar
USE solvar
USE timevar
IMPLICIT NONE
INTEGER, INTENT(IN) :: g
INTEGER :: ierr, it, i, j, ii, jj, k
REAL*8 :: bnorm, temp, rnorm, df
REAL*8, DIMENSION(rn) :: cs, sn, y
REAL*8, DIMENSION(rn1) :: e1, sh, sht
REAL*8, DIMENSION(neq) :: rs, ws
REAL*8, DIMENSION(bcs,8) :: ra, wa
REAL*8, DIMENSION(rn1,rn) :: h
REAL*8, DIMENSION(neq,rn1) :: vs
REAL*8, DIMENSION(bcs,8) :: vai
REAL*8, DIMENSION(bcs,8,rn1) :: va

INCLUDE 'mpif.h'

! The large GMRES matrix is distributed. Non-zero elements are the operators.
! The solution vector is made up of f and psio, both distributed. r is distributed similarly.
f(:,g) = 0.0
psio = 0.0
rs = bs
ra = ba

! Get rnorm depending on the preconditioner flags
IF (pcf == 0) THEN
   ! bs/ba make up full b. b is distributed.
   bnorm = DOT_PRODUCT(bs,bs)
   DO i = 1, 8
      bnorm = bnorm + DOT_PRODUCT(ba(:,i),ba(:,i))
   END DO
   CALL MPI_ALLREDUCE(bnorm,temp,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
   bnorm = SQRT(temp)
   IF (bnorm == 0.0) THEN
      IF (irank == root) WRITE(8,'(/,2X,A)') "ERROR: RHS vector all zero."
      CALL MPI_FINALIZE(ierr)
      STOP
   END IF
   rnorm = bnorm
ELSE
   CALL precond(pinv,tpose,rs)
   ! Now set rnorm directly
   rnorm = DOT_PRODUCT(rs,rs)
   DO i = 1, 8
      rnorm = rnorm + DOT_PRODUCT(ra(:,i),ra(:,i))
   END DO
   CALL MPI_ALLREDUCE(rnorm,temp,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
   rnorm = SQRT(temp)
END IF

! Other initializations
vs = 0.0
va = 0.0
h = 0.0
cs = 0.0
sn = 0.0
e1 = 0.0
e1(1) = 1.0

!---------------------------------------------------------------------

DO it = 1, itmx
   vs(:,1) = (1.0/rnorm)*rs
   va(:,:,1) = (1.0/rnorm)*ra(:,:)
   sh = rnorm*e1

   ! Start the restart loop
   DO j = 1, rn
      ! Need to communicate elements of current V column to complete matvec op.
      ! Data communicated exactly like psio/psii in PBJ method, so only need va
      CALL gridcomm(va(:,:,j),vai)

      ! Call for matrix vector multiplication. Operations follow those of IDOT from PBJ code.
      CALL avmult(vs(:,j),va(:,:,j),vai,ws,wa)

      ! Apply preconditioner as necessary
      IF (pcf == 1) CALL precond(pinv,tpose,ws)

      ! Orthogonalization phase with Gram-Schmidt
      DO i = 1, j
         h(i,j) = DOT_PRODUCT(ws,vs(:,i))
         DO k = 1,8
            h(i,j) = h(i,j) + DOT_PRODUCT(wa(:,k),va(:,k,i))
         END DO
         CALL MPI_ALLREDUCE(h(i,j),temp,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
         h(i,j) = temp
         ws = ws - h(i,j)*vs(:,i)
         wa = wa - h(i,j)*va(:,:,i)
      END DO
      ! Compute final element in current column of H
      h(j+1,j) = DOT_PRODUCT(ws,ws)
      DO k = 1, 8
         h(j+1,j) = h(j+1,j) + DOT_PRODUCT(wa(:,k),wa(:,k))
      END DO
      CALL MPI_ALLREDUCE(h(j+1,j),temp,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      h(j+1,j) = SQRT(temp)
      ! Form the next vector of V
      vs(:,j+1) = (1.0/h(j+1,j))*ws
      va(:,:,j+1) = (1.0/h(j+1,j))*wa

      ! Perform Givens rotations on H to set up least squares problem
      DO i = 1, j-1
         temp     =  cs(i)*h(i,j) + sn(i)*h(i+1,j)
         h(i+1,j) = -sn(i)*h(i,j) + cs(i)*h(i+1,j)
         h(i,j)   = temp
      END DO
      CALL rotmat(h(j,j),h(j+1,j),cs(j),sn(j))
      temp    =  cs(j)*sh(j)
      sh(j+1) = -sn(j)*sh(j)
      sh(j)   = temp
      h(j,j) = cs(j)*h(j,j) + sn(j)*h(j+1,j)
      h(j+1,j) = 0.0

      ! Set the error/residual norm to the final value in the s vector
      df = ABS(sh(j+1))

      ! Solve the least squares problem on the upper triangular portion of H if
      ! error is low enough
      sht = sh
      IF (df <= err) THEN
         DO jj = j, 1, -1
            IF (h(jj,jj) == 0.0) THEN
               IF (irank == root) WRITE(8,'(/,2X,A)') "ERROR: zero value on H diagonal."
               CALL MPI_FINALIZE(ierr)
               STOP
            END IF
            y(jj) = sht(jj)/h(jj,jj)
            DO ii = 1, jj-1
               sht(ii) = sht(ii) - h(ii,jj)*y(jj)
            END DO
         END DO
         ! And since it is low enough, update the solution
         f(:,g) = f(:,g) + MATMUL(vs(:,1:j),y(1:j))
         DO k = 1, 8
            psio(:,k) = psio(:,k) + MATMUL(va(:,k,1:j),y(1:j))
         END DO

         ! Print messages
         IF (irank == root) THEN
            WRITE(8,'(/,2X,A,I3,A,I5,A,ES11.3)') "Group ", g, " converged in ", (it-1)*rn+j, " iterations."
            WRITE(8,'(2X,A,I5)') "Number of restarts: ", (it-1)
            WRITE(8,'(2X,A,ES11.3)') "Residual norm = ", df
         END IF
         EXIT
      ELSE
         IF (itp == 1 .AND. irank == root) WRITE(8,111) g, (it-1)*rn+j, df
      END IF
   END DO

   ! Break again if converged enough
   IF (df <= err) EXIT

   ! Prepare for next iteration
   DO j = rn, 1, -1
      IF (h(j,j) == 0.0) THEN
         IF (irank == root) WRITE(8,'(/,2X,A)') "ERROR: zero value on H diagonal."
         CALL MPI_FINALIZE(ierr)
         STOP
      END IF
      y(j) = sh(j)/h(j,j)
      DO i = 1, j-1
         sh(i) = sh(i) - h(i,j)*y(j)
      END DO
   END DO

   f(:,g) = f(:,g) + MATMUL(vs(:,1:rn),y(1:rn))
   DO k = 1, 8
      psio(:,k) = psio(:,k) + MATMUL(va(:,k,1:rn),y(1:rn))
   END DO

   ! Communicate psio to psii
   CALL gridcomm(psio,psii)

   ! To compute a new residual, need to do new multiplications of ITMM operators
   ! and solution vectors
   CALL avmult(f(:,g),psio,psii,rs,ra)
   rs = bs - rs
   IF (pcf == 1) CALL precond(pinv,tpose,rs) 
   ra = ba - ra
   rnorm = df
END DO

IF (df > err) THEN
   IF (irank == root) WRITE(8,'(/,2X,A,I2,A,ES11.3)') "ERROR: Group ", g, " failed to converged. Residual norm = ", rnorm
   CALL MPI_FINALIZE(ierr)
   STOP
END IF

111 FORMAT(2X,'Gr',I3,' It ',I5,' Dfmx ',ES11.3)

RETURN
END SUBROUTINE gmresmgs
