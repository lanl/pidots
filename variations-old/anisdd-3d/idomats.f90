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
INTEGER :: ieq, i, j, k, info, t, l, ll
REAL*8 :: symm
REAL*8, DIMENSION(:), ALLOCATABLE :: q, wrk, symmvec

! Initialize the matrices
kmat = 0.0
jpsi = 0.0
kpsi = 0.0

! Allocate the iteration Jacobian
ALLOCATE(jmat(neq,neq))
jmat = 0.0

! Allocate the temp arrays
ALLOCATE(q(neq))

! Allocate idos-dependent arrays
IF (idos == 0) THEN
   ALLOCATE(piv(neq), wrk(neq))
END IF

IF (sym == 1 .OR. idos == 1) THEN
   ALLOCATE(symmvec(neq))
END IF

! Can now start the single sweep for constructing all matrices
CALL matsweep(g)

! Begin constructing the RHS
DO k = 1, nz
   DO j = 1, ny
      DO i = 1, nx
         t = ((i-1) + (j-1)*nx + (k-1)*nx*ny)*nmom
         DO l = 0, anord
            ieq = t + l**2 + 1
            src(ieq) = sm(ieq,g)/sigs(l,mat(i,j,k),g,g)
            DO ll = 1, l
               ieq = t + l**2 + 2*ll
               src(ieq) = sm(ieq,g)/sigs(l,mat(i,j,k),g,g)
               src(ieq+1) = sm(ieq+1,g)/sigs(l,mat(i,j,k),g,g)
            END DO
         END DO
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

! Symmetrize the matrix if necessary
IF (sym == 1 .OR. idos == 1) CALL symtrz(g,symmvec)

!--------------------------------------------------------------------------------------------------
! Other operations depend on the IDO Solver value: 0-Direct, 1-CG
! If idos=0, need to invert the current jmat and multiply q and kmat
IF (idos == 0) THEN
   IF (sym == 0) THEN
      ! Invert the new LHS matrix
      CALL dgetrf(neq,neq,jmat,neq,piv,info)
      IF (info /= 0) THEN
         WRITE (8,'(//,1X,A)') "ERROR: jmat either has illegal value or is singular, cannot factor."
         STOP
      END IF

      ! Set sv with q. Have needed values src, sv, jmat, kmat, jpsi, kpsi, piv
      sv = q
      DEALLOCATE(q,wrk)

!---------------------------------------------------------------
! Symmetric direct
   ELSE
      ! Invert the new LHS matrix with symmetric routine; not sure about
      ! positive-definiteness-- use dsytrf instead of dpotrf
      CALL dsytrf('U',neq,jmat,neq,piv,wrk,neq,info)
      IF (info /= 0) THEN
         WRITE (8,'(//,1X,A)') "ERROR: jmat either has illegal value or is singular, cannot factor."
         STOP
      END IF

      ! Multiply q by the symmvec and set sv
      ! Then multiply kmat's by symmvec
      sv = symmvec*q
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

      ! Have needed values: src, sv, jmat, kmat, jpsi, kpsi, piv
      DEALLOCATE(q,symmvec,wrk)

   END IF
!------------------------------------------------------
! CG
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

   ! Now have the six needed values: src, sv, jmat, kmat, jpsi, kpsi
   DEALLOCATE(q,symmvec)
END IF

RETURN
END SUBROUTINE idomats
