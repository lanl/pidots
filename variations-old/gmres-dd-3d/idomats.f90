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
INTEGER :: ieq, i, j, k, info
INTEGER :: ix1, ix2, ix3, ix4, ix5, ix6, kz, lz, ky, ly, kx, lx, jndx
REAL*8, DIMENSION(neq) :: q, wrk

! Initialize the matrices and vectors
jmat = 0.0
kmat = 0.0
jpsi = 0.0
kpsi = 0.0
bs   = 0.0
ba   = 0.0

! Can now start the single sweep for constructing all matrices
CALL matsweep(g)

! Begin constructing the RHS
DO k = 1, nz
   DO j = 1, ny
      DO i = 1, nx
         ieq = i + (j-1)*nx + (k-1)*xys
         q(ieq) = s(i,j,k,g)/sigs((mat(i,j,k)),g,g)
      END DO
   END DO
END DO

! Multiply q by Jphi to form the RHS vector related to phi equations
IF (tpose == 1) THEN
   DO i = 1, neq
      bs(i) = DOT_PRODUCT(jmat(:,i),q)
   END DO
ELSE
   bs = MATMUL(jmat,q)
END IF

! Multiply q by Jpsi to form the RHS vector related to psio equations
IF (tpose == 1) THEN
   DO j = 1, 8
      DO i = 1, bcs
         ba(i,j) = DOT_PRODUCT(jpsi(:,i,j),q)
      END DO
   END DO
ELSE
   DO j = 1, 8
      ba(:,j) = MATMUL(jpsi(:,:,j),q)
   END DO
END IF

! Now add to RHS vectors contribution from non-zero *fixed* boundary conditions
! All P will do for balance, even if not at global domain boundary (where BCs apply)

! Do only if a BC value is 2, otherwise don't do unnecessary operations
IF (MAXVAL(bc) == 2) THEN

   ! bs updates
   IF (tpose == 1) THEN
      DO j = 1, 8
         DO i = 1, neq
            bs(i) = bs(i) + DOT_PRODUCT(kmat(:,i,j),psii(:,j))
         END DO
      END DO
   ELSE
      bs = bs + MATMUL(kmat(:,:,1),psii(:,1)) + MATMUL(kmat(:,:,2),psii(:,2)) + MATMUL(kmat(:,:,3),psii(:,3)) & 
              + MATMUL(kmat(:,:,4),psii(:,4)) + MATMUL(kmat(:,:,5),psii(:,5)) + MATMUL(kmat(:,:,6),psii(:,6)) &
              + MATMUL(kmat(:,:,7),psii(:,7)) + MATMUL(kmat(:,:,8),psii(:,8))
   END IF

   ! ba updates
   ! Set up the loop independent variables
   ix1 = 1
   ix2 = xys
   ix3 = ix2 + 1
   ix4 = ix2 + xzs
   ix5 = ix4 + 1
   ix6 = ix4 + yzs
   IF (tpose == 1) THEN
      DO j = 1, 8
         DO i = 1, apo
            kz = (i-1)*xys + 1
            lz = i*xys
            ky = apo*xys + (i-1)*xzs + 1
            ly = apo*xys + i*xzs
            kx = apo*(xys+xzs) + (i-1)*yzs + 1
            lx = apo*(xys+xzs) + i*yzs
            jndx = 0
            DO k = kz, lz
               jndx = jndx + 1
               ba(k,j) = ba(k,j) + DOT_PRODUCT(kpsi(ix1:ix2,jndx,i,j),psii(kz:lz,j)) &
                                 + DOT_PRODUCT(kpsi(ix3:ix4,jndx,i,j),psii(ky:ly,j)) &
                                 + DOT_PRODUCT(kpsi(ix5:ix6,jndx,i,j),psii(kx:lx,j))
            END DO
            DO k = ky, ly
               jndx = jndx + 1
               ba(k,j) = ba(k,j) + DOT_PRODUCT(kpsi(ix1:ix2,jndx,i,j),psii(kz:lz,j)) &
                                 + DOT_PRODUCT(kpsi(ix3:ix4,jndx,i,j),psii(ky:ly,j)) &
                                 + DOT_PRODUCT(kpsi(ix5:ix6,jndx,i,j),psii(kx:lx,j))
            END DO
            DO k = kx, lx
               jndx = jndx + 1
               ba(k,j) = ba(k,j) + DOT_PRODUCT(kpsi(ix1:ix2,jndx,i,j),psii(kz:lz,j)) &
                                 + DOT_PRODUCT(kpsi(ix3:ix4,jndx,i,j),psii(ky:ly,j)) &
                                 + DOT_PRODUCT(kpsi(ix5:ix6,jndx,i,j),psii(kx:lx,j))
            END DO
         END DO
      END DO
   ELSE
      DO j = 1, 8
         DO i = 1, apo
            kz = (i-1)*xys + 1
            lz = i*xys
            ky = apo*xys + (i-1)*xzs + 1
            ly = apo*xys + i*xzs
            kx = apo*(xys+xzs) + (i-1)*yzs + 1
            lx = apo*(xys+xzs) + i*yzs
            ba(kz:lz,j) = ba(kz:lz,j) + MATMUL(kpsi(ix1:ix2,ix1:ix2,i,j),psii(kz:lz,j)) &
                                      + MATMUL(kpsi(ix1:ix2,ix3:ix4,i,j),psii(ky:ly,j)) &
                                      + MATMUL(kpsi(ix1:ix2,ix5:ix6,i,j),psii(kx:lx,j))

            ba(ky:ly,j) = ba(ky:ly,j) + MATMUL(kpsi(ix3:ix4,ix1:ix2,i,j),psii(kz:lz,j)) &
                                      + MATMUL(kpsi(ix3:ix4,ix3:ix4,i,j),psii(ky:ly,j)) &
                                      + MATMUL(kpsi(ix3:ix4,ix5:ix6,i,j),psii(kx:lx,j))

            ba(kx:lx,j) = ba(kx:lx,j) + MATMUL(kpsi(ix5:ix6,ix1:ix2,i,j),psii(kz:lz,j)) &
                                      + MATMUL(kpsi(ix5:ix6,ix3:ix4,i,j),psii(ky:ly,j)) &
                                      + MATMUL(kpsi(ix5:ix6,ix5:ix6,i,j),psii(kx:lx,j))
         END DO
      END DO
   END IF

END IF

!-----------------------------------------------------------------
! Form the LHS matrix
jmat = -jmat
DO i = 1, neq
   jmat(i,i) = 1.0 + jmat(i,i)
END DO

! Form the diagonal block preconditioner if pcf = 1
IF (pcf == 1) THEN
   ALLOCATE(prec(neq,neq), piv(neq))
   prec = jmat

   ! Factor the preconditioner matrix regardless
   CALL dgetrf(neq,neq,prec,neq,piv,info)
   IF (info /= 0) THEN
      WRITE (8,'(//,1X,A)') "ERROR: jmat either has illegal value or is singular, cannot factor."
      STOP
   END IF

   ! Next check pinv to determine if matrix inversion is necessary
   IF (pinv == 1) THEN
      CALL dgetri(neq,prec,neq,piv,wrk,neq,info)
      IF (info /= 0) THEN
         WRITE (8,'(//,1X,A)') "ERROR: jmat either has illegal value or is singular, cannot invert."
         STOP
      END IF
      DEALLOCATE(piv)
   END IF
END IF

RETURN
END SUBROUTINE idomats
