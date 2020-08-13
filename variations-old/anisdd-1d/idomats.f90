SUBROUTINE idomats(g,nsm)

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
INTEGER, INTENT(IN) :: g, nsm
INTEGER :: i, j, l, info, m, tt, t, ieq
INTEGER, DIMENSION(:), ALLOCATABLE :: piv
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
ALLOCATE(jmat(nsm,nsm))
jmat = 0.0

! Allocate the temp arrays
ALLOCATE(q(nsm))

! Allocate idos-dependent arrays
ALLOCATE(piv(nsm), wrk(nsm))

IF (sym == 1 .OR. idos == 1) THEN
   ALLOCATE(symmvec(nsm))
END IF

! Can now start the single sweep for constructing all matrices
CALL matsweep(g)

! Begin constructing the RHS
DO i = 1, nx
   DO l = 0, anord
      ieq = (i-1)*sord + 1 + l
      src(ieq) = sm(ieq,g)/sigs(mat(i),l,g,g)
   END DO
END DO

! Multiply the source vector by the Jacobi iteration matrix
IF (tpose == 1) THEN
   DO j = 1, nsm
      q(j) = DOT_PRODUCT(jmat(:,j),src)
   END DO
ELSE
   q = MATMUL(jmat,src)
END IF

! Form the LHS matrix
jmat = -jmat
DO t = 1, nsm
   jmat(t,t) = 1.0 + jmat(t,t)
END DO

IF (sym == 1 .OR. idos == 1) THEN
   IF (tpose == 1) THEN
      ! Symmetrize jmat and put q into sv
      DO i = 1, nx
         m = mat(i)
         DO l = 0, anord
            ieq = (i-1)*sord + 1 + l
            symm = sigs(m,l,g,g)*dx(i)*(2.0*l+1.0)
            IF (MOD(l,2) == 1) symm = -1.0*symm        ! Odd order has the negative mult.
            symmvec(ieq) = symm
            DO tt = 1, nsm
               jmat(tt,ieq) = symm*jmat(tt,ieq)
            END DO
         END DO
      END DO
   ELSE
      DO i = 1, nx
         m = mat(i)
         DO l = 0, anord
            ieq = (i-1)*sord + 1 + l
            symm = sigs(m,l,g,g)*dx(i)*(2.0*l+1.0)
            IF (MOD(l,2) == 1) symm = -1.0*symm        ! Odd order has the negative mult.
            symmvec(ieq) = symm
            DO tt = 1, nsm
               jmat(ieq,tt) = symm*jmat(ieq,tt)
            END DO
         END DO
      END DO
   END IF
END IF

! Other operations depend on the IDO Solver value: 0-Direct, 1-CG
! If idos=0, need to invert the current jmat and multiply q and kmat
IF (idos == 0) THEN
   IF (sym == 0) THEN
      ! Invert the new LHS matrix
      CALL dgetrf(nsm,nsm,jmat,nsm,piv,info)
      IF (info /= 0) THEN
         WRITE (8,'(//,1X,A)') "ERROR: jmat either has illegal value or is singular, cannot factor."
         STOP
      END IF
      CALL dgetri(nsm,jmat,nsm,piv,wrk,nsm,info)
      IF (info /= 0) THEN
         WRITE (8,'(//,1X,A)') "ERROR: jmat either has illegal value or is singular, cannot invert."
         STOP
      END IF

      ! Save the matrix vector product of the inverted matrix and the modified source
      ! Then save the product of the inverted matrix and kmat
      IF (tpose == 0) THEN
         sv = MATMUL(jmat,q)
         ALLOCATE(tmat(nsm,apo))
         DO t = 1, 2
            CALL DGEMM('n','n',nsm,apo,nsm,one,jmat,nsm,kmat(:,:,t),nsm,zero,tmat,nsm)
            kmat(:,:,t) = tmat
         END DO
      ELSE
         DO j = 1, nsm
            sv(j) = DOT_PRODUCT(jmat(:,j),q)
         END DO
         ! kmat*jmat
         ALLOCATE(tmat(apo,nsm))
         DO t = 1, 2
            CALL DGEMM('n','n',apo,nsm,nsm,one,kmat(:,:,t),apo,jmat,nsm,zero,tmat,apo)
            kmat(:,:,t) = tmat
         END DO
      END IF

      ! Now have the five needed values: src, sv, kmat, jpsi, kpsi
      DEALLOCATE(jmat,tmat)
      DEALLOCATE(piv,q,wrk)
!---------------------------------------------------------------
! Symmetric direct
   ELSE
      ! Invert the new LHS matrix with symmetric routine.
      ! Symmetric INDEFINITE routines used. Not positive definite system.
      CALL dsytrf('U',nsm,jmat,nsm,piv,wrk,nsm,info)
      IF (info /= 0) THEN
         WRITE (8,'(//,1X,A)') "ERROR: jmat either has illegal value or is singular, cannot factor."
         STOP
      END IF
      CALL dsytri('U',nsm,jmat,nsm,piv,wrk,info)
      IF (info /= 0) THEN
         WRITE (8,'(//,1X,A)') "ERROR: jmat either has illegal value or is singular, cannot invert."
         STOP
      END IF

      ! Symmetric inverter only returns triangular portion of invert, copy rest
      DO j = 1, nsm
         DO i = 1, j
            jmat(j,i) = jmat(i,j)
         END DO
      END DO

      ! Multiply q by the symmvec and set sv
      ! Then multiply kmat's by symmvec
      ! Lastly reset kmat as jmat*kmat
      q = symmvec*q
      IF (tpose == 0) THEN
         sv = MATMUL(jmat,q)
         DO t = 1, 2
            DO i = 1, nsm
               kmat(i,:,t) = symmvec(i)*kmat(i,:,t)
            END DO
         END DO
         ALLOCATE(tmat(nsm,apo))
         DO t = 1, 2
            CALL DGEMM('n','n',nsm,apo,nsm,one,jmat,nsm,kmat(:,:,t),nsm,zero,tmat,nsm)
            kmat(:,:,t) = tmat
         END DO
      ELSE
         DO j = 1, nsm
            sv(j) = DOT_PRODUCT(jmat(:,j),q)
         END DO
         DO t = 1, 2
            DO j = 1, nsm
               kmat(:,j,t) = symmvec(j)*kmat(:,j,t)
            END DO
         END DO
         ! Aided a bit by the symmetry of jmat, use that feature.
         ! tmat comes out as the normal kmat shape instead of transposed
         ALLOCATE(tmat(nsm,apo))
         DO t = 1, 2
            CALL DGEMM('n','t',nsm,apo,nsm,one,jmat,nsm,kmat(:,:,t),apo,zero,tmat,nsm)
            kmat(:,:,t) = TRANSPOSE(tmat)
         END DO
      END IF

      ! Now have the five needed values: src, sv, kmat, jpsi, kpsi
      DEALLOCATE(jmat,tmat,q,symmvec,piv,wrk)
   END IF
!------------------------------------------------------
! CG
ELSE
   ! Reset kmat with the symmvec
   IF (tpose == 1) THEN
      DO t = 1, 2
         DO j = 1, nsm
            kmat(:,j,t) = symmvec(j)*kmat(:,j,t)
         END DO
      END DO
   ELSE
      DO t = 1, 2
         DO i = 1, nsm
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
