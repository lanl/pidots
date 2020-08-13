SUBROUTINE idomats(g,neq,bcs)

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
INTEGER, INTENT(IN) :: g, neq, bcs
INTEGER :: ieq, i, j, info, m, tt, t, l, ll, indx
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
DO j = 1, ny
   DO i = 1, nx
      DO l = 0, anord
         DO ll = 0, l
            ieq = ((i-1) + (j-1)*nx)*nmom + l*(l+1)/2 + ll + 1
            src(ieq) = sm(ieq,g)/sigs(l,mat(i,j),g,g)
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

IF (sym == 1 .OR. idos == 1) THEN
   IF (tpose == 1) THEN
      ! Symmetrize jmat and put q into sv
      DO j = 1, ny
         DO i = 1, nx
            m = mat(i,j)
            indx = ((i-1) + (j-1)*nx)*nmom
            DO l = 0, anord
               symm = sigs(l,m,g,g)*dx(i)*dy(j)*((-1.0)**l)
               DO ll = 0, l
                  ieq = indx + l*(l+1)/2 + ll + 1
                  IF (ll > 0) THEN
                     symmvec(ieq) = symm*2.0
                  ELSE
                     symmvec(ieq) = symm
                  END IF
                  DO tt = 1, neq
                     jmat(tt,ieq) = symmvec(ieq)*jmat(tt,ieq)
                  END DO
               END DO
            END DO
         END DO
      END DO
   ELSE
      DO j = 1, ny
         DO i = 1, nx
            m = mat(i,j)
            indx = ((i-1) + (j-1)*nx)*nmom
            DO l = 0, anord
               symm = sigs(l,m,g,g)*dx(i)*dy(j)*((-1.0)**l)
               DO ll = 0, l
                  ieq = indx + l*(l+1)/2 + ll + 1
                  IF (ll > 0) THEN
                     symmvec(ieq) = 2.0*symm
                  ELSE
                     symmvec(ieq) = symm
                  END IF
                  DO tt = 1, neq
                     jmat(ieq,tt) = symmvec(ieq)*jmat(ieq,tt)
                  END DO
               END DO
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
      ! Then save the product of the inverted matrix and kmat
      IF (tpose == 0) THEN
         sv = MATMUL(jmat,q)
         ALLOCATE(tmat(neq,bcs))
         DO t = 1, 4
            CALL DGEMM('n','n',neq,bcs,neq,one,jmat,neq,kmat(:,:,t),neq,zero,tmat,neq)
            kmat(:,:,t) = tmat
            ! Old loop method abandoned for DGEMM. Keep here just in case.
            !DO k = 1, bcs
            !   c = 0.0
            !   DO j = 1, neq
            !      DO i = 1, neq
            !         c(i) = c(i) + jmat(i,j)*kmat(j,k,t)
            !      END DO
            !   END DO
            !   kmat(:,k,t) = c(:)
            !END DO
         END DO
      ELSE
         DO j = 1, neq
            sv(j) = DOT_PRODUCT(jmat(:,j),q)
         END DO
         ! kmat*jmat
         ALLOCATE(tmat(bcs,neq))
         DO t = 1, 4
            CALL DGEMM('n','n',bcs,neq,neq,one,kmat(:,:,t),bcs,jmat,neq,zero,tmat,bcs)
            kmat(:,:,t) = tmat
            ! Old loop method abandoned for DGEMM. Keep here just in case.
            !DO k = 1, bcs
            !   c = kmat(k,:,t)
            !   DO j = 1, neq
            !      kmat(k,j,t) = DOT_PRODUCT(c,jmat(:,j))
            !   END DO
            !END DO
         END DO
      END IF

      ! Now have the five needed values: src, sv, kmat, jpsi, kpsi
      DEALLOCATE(jmat,tmat)
      DEALLOCATE(piv,q,wrk)
!---------------------------------------------------------------
! Symmetric direct
   ELSE
      ! Invert the new LHS matrix with symmetric routine
      CALL dsytrf('U',neq,jmat,neq,piv,wrk,neq,info)
      IF (info /= 0) THEN
         WRITE (8,'(//,1X,A)') "ERROR: jmat either has illegal value or is singular, cannot factor."
         STOP
      END IF
      CALL dsytri('U',neq,jmat,neq,piv,wrk,info)
      IF (info /= 0) THEN
         WRITE (8,'(//,1X,A)') "ERROR: jmat either has illegal value or is singular, cannot invert."
         STOP
      END IF

      ! Symmetric inverter only returns triangular portion of invert, copy rest
      DO j = 1, neq
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
         DO t = 1, 4
            DO i = 1, neq
               kmat(i,:,t) = symmvec(i)*kmat(i,:,t)
            END DO
         END DO
         ALLOCATE(tmat(neq,bcs))
         DO t = 1, 4
            CALL DGEMM('n','n',neq,bcs,neq,one,jmat,neq,kmat(:,:,t),neq,zero,tmat,neq)
            kmat(:,:,t) = tmat
            ! Old loop method abandoned for DGEMM. Keep here just in case.
            !DO k = 1, bcs
            !   c = 0.0
            !   DO j = 1, neq
            !      DO i = 1, neq
            !         c(i) = c(i) + jmat(i,j)*kmat(j,k,t)
            !      END DO
            !   END DO
            !   kmat(:,k,t) = c(:)
            !END DO
         END DO
      ELSE
         DO j = 1, neq
            sv(j) = DOT_PRODUCT(jmat(:,j),q)
         END DO
         DO t = 1, 4
            DO j = 1, neq
               kmat(:,j,t) = symmvec(j)*kmat(:,j,t)
            END DO
         END DO
         ! Aided a bit by the symmetry of jmat, use that feature.
         ! tmat comes out as the normal kmat shape instead of transposed
         ALLOCATE(tmat(neq,bcs))
         DO t = 1, 4
            CALL DGEMM('n','t',neq,bcs,neq,one,jmat,neq,kmat(:,:,t),bcs,zero,tmat,neq)
            kmat(:,:,t) = TRANSPOSE(tmat)
            ! Old loop method abandoned for DGEMM. Keep here just in case.
            !DO k = 1, bcs
            !   c = kmat(k,:,t)
            !   DO j = 1, neq
            !      kmat(k,j,t) = DOT_PRODUCT(c,jmat(:,j))
            !   END DO
            !END DO
         END DO
      END IF

      DEALLOCATE(piv,wrk)
      ! Now have the five needed values: src, sv, kmat, jpsi, kpsi
      DEALLOCATE(jmat,tmat,q,symmvec)
   END IF
!------------------------------------------------------
! CG
ELSE
   ! Reset kmat with the symmvec
   IF (tpose == 1) THEN
      DO t = 1, 4
         DO j = 1, neq
            kmat(:,j,t) = symmvec(j)*kmat(:,j,t)
         END DO
      END DO
   ELSE
      DO t = 1, 4
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
