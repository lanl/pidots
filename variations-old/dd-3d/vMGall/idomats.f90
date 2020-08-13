SUBROUTINE idomats

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
INTEGER :: ieq, i, j, k, info, m, tt, t
REAL*8 :: symm, one, zero
REAL*8, DIMENSION(:), ALLOCATABLE :: q, wrk, symmvec
REAL*8, DIMENSION(:,:), ALLOCATABLE :: tmat

! Initialize the matrices
kmat = 0.0
jpsi = 0.0
kpsi = 0.0

! Set one and zero
one = 1.0
zero = 0.0

! Allocate the iteration Jacobian
ALLOCATE(jmat(neq,neq))
jmat = 0.0

! Allocate the temp source array
ALLOCATE(q(neq))

!-------------------------------------------------------------------
! Can now start the single sweep for constructing all matrices
CALL matsweep

! Begin constructing the RHS
DO k = 1, nz
   DO j = 1, ny
      DO i = 1, nx
         ieq = i + (j-1)*nx + (k-1)*xys
         src(ieq) = s(i,j,k)/sigs(mat(i,j,k))
      END DO
   END DO
END DO

! Multiply the source vector by the Jacobi iteration matrix
IF (tpose == 1) THEN
   DO j = 1, neq
      q(j) = DOT_PRODUCT(jmat(:,j),src)
   END DO
ELSE
   q = MATMUL(jmat,src)
END IF

! Form the LHS matrix
jmat = -jmat
DO t = 1, neq
   jmat(t,t) = 1.0 + jmat(t,t)
END DO

!------------------------------------------------------------------
IF (sym == 1 .OR. idos == 1 .OR. idos == 2) THEN
   ALLOCATE(symmvec(neq))
   IF (tpose == 1) THEN
      ! Symmetrize jmat and put q into sv
      DO k = 1, nz
         DO j = 1, ny
            DO i = 1, nx
               m = mat(i,j,k)
               symm = sigs(m)*dx(i)*dy(j)*dz(k)
               ieq = (k-1)*xys + (j-1)*nx + i
               symmvec(ieq) = symm
               DO tt = 1, neq
                  jmat(tt,ieq) = symmvec(ieq)*jmat(tt,ieq)
               END DO
            END DO
         END DO
      END DO
   ELSE
      DO k = 1, nz
         DO j = 1, ny
            DO i = 1, nx
               m = mat(i,j,k)
               symm = sigs(m)*dx(i)*dy(j)*dz(k)
               ieq = (k-1)*xys + (j-1)*nx + i
               symmvec(ieq) = symm
               DO tt = 1, neq
                  jmat(ieq,tt) = symmvec(ieq)*jmat(ieq,tt)
               END DO
            END DO
         END DO
      END DO
   END IF
END IF

!------------------------------------------------------------------------
! Other operations depend on the IDO Solver value: 0-Direct, 1-CG, 2-SOR
! If idos=0, need to factor jmat
IF (idos == 0) THEN

   ! Call LAPACK routines according to sym flag
   IF (sym == 0) THEN

      ! Allocate the piv vector as it's now needed
      ALLOCATE(piv(neq))

      ! Always factor the jmat matrix
      CALL dgetrf(neq,neq,jmat,neq,piv,info)
      IF (info /= 0) THEN
         WRITE (8,'(//,1X,A)') "ERROR: jmat either has illegal value or is singular, cannot factor."
         STOP
      END IF

      ! Set sv with q. Have needed values src, sv, jmat, kmat, jpsi, kpsi, piv
      sv = q
      DEALLOCATE(q)

   !---------------------------------------------------------------
   ! Symmetric direct
   ELSE

      ! Always factor the new LHS matrix with symmetric routine
      CALL dpotrf('U',neq,jmat,neq,info)
      IF (info /= 0) THEN
         WRITE (8,'(//,1X,A)') "ERROR: jmat either has illegal value or is singular, cannot factor."
         STOP
      END IF

      ! And always apply symmvec to q and to kmat
      q = symmvec*q
      IF (tpose == 1) THEN
         DO t = 1, 8
            DO j = 1, neq
               kmat(:,j,t) = symmvec(j)*kmat(:,j,t)
            END DO
         END DO
      ELSE
         DO t = 1, 8
            DO i = 1, neq
               kmat(i,:,t) = symmvec(i)*kmat(i,:,t)
            END DO
         END DO
      END IF

      ! Set sv as q. Have needed values src, sv, jmat, kmat, jpsi, kpsi
      sv = q
      DEALLOCATE(q,symmvec)

   END IF

!-----------------------------------------------------------------------
! CG or SOR
ELSE
   ! Reset kmat with the symmvec
   IF (tpose == 1) THEN
      DO t = 1, 8
         DO j = 1, neq
            kmat(:,j,t) = symmvec(j)*kmat(:,j,t)
         END DO
      END DO
   ELSE
      DO t = 1, 8
         DO i = 1, neq
            kmat(i,:,t) = symmvec(i)*kmat(i,:,t)
         END DO
      END DO
   END IF

   ! Set the sv vector
   sv = symmvec*q

   ! If a preconditioner is to be used, the inverse diagonal will be needed
   IF (idos == 1 .AND. pcf == 1) THEN
      ALLOCATE(dv(neq))
      DO i = 1, neq
         dv(i) = 1.0/jmat(i,i)
      END DO
      ! Check if the type is SSOR
      IF (pcty == 2) dv = sorw*dv
   END IF

   ! Now have the six needed values: src, sv, jmat, kmat, jpsi, kpsi (and dv)
   DEALLOCATE(q,symmvec)
END IF

RETURN
END SUBROUTINE idomats
