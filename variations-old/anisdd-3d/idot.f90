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
INTEGER :: i, j, k, ix1, ix2, ix3, ix4, ix5, ix6, kz, lz, ky, ly, kx, lx, jndx, info
INTEGER, INTENT(OUT) :: its
REAL*8, DIMENSION(neq) :: q

! Compute the solution directly or with CG iterations
IF (idos == 0) THEN
   ! Solve for the new phi with the sum of sv with the product of kmat and psii
   ! Do according to 'tpose' flag
   IF (tpose == 1) THEN
      phi(:,g) = sv
      DO j = 1, 8
         DO i = 1, neq
            phi(i,g) = phi(i,g) + DOT_PRODUCT(kmat(:,i,j),psii(:,j))
         END DO
      END DO
   ELSE
      phi(:,g) = sv + MATMUL(kmat(:,:,1),psii(:,1)) + MATMUL(kmat(:,:,2),psii(:,2)) + MATMUL(kmat(:,:,3),psii(:,3)) &
                    + MATMUL(kmat(:,:,4),psii(:,4)) + MATMUL(kmat(:,:,5),psii(:,5)) + MATMUL(kmat(:,:,6),psii(:,6)) &
                    + MATMUL(kmat(:,:,7),psii(:,7)) + MATMUL(kmat(:,:,8),psii(:,8))
   END IF

   ! Use LAPACK solvers
   IF (tpose == 1) THEN
      IF (sym == 1) THEN
         CALL dsytrs('U',neq,1,jmat,neq,piv,phi(:,g),neq,info)
         IF (info /= 0) THEN
            WRITE (8,'(//,1X,A)') "ERROR: dsytrs cannot solve."
            STOP
         END IF
      ELSE
         CALL dgetrs('T',neq,1,jmat,neq,piv,phi(:,g),neq,info)
         IF (info /= 0) THEN
            WRITE (8,'(//,1X,A)') "ERROR: dgetrs cannot solve."
            STOP
         END IF
      END IF
   ELSE
      IF (sym == 1) THEN
         CALL dsytrs('U',neq,1,jmat,neq,piv,phi(:,g),neq,info)
         IF (info /= 0) THEN
            WRITE (8,'(//,1X,A)') "ERROR: dsytrs cannot solve."
            STOP
         END IF
      ELSE
         CALL dgetrs('N',neq,1,jmat,neq,piv,phi(:,g),neq,info)
         IF (info /= 0) THEN
            WRITE (8,'(//,1X,A)') "ERROR: dgetrs cannot solve."
            STOP
         END IF
      END IF
   END IF
   cnvf(g) = 1
   its = 0
ELSE
   ! Solve for the new phi with CG iterations from the cgsolve routinei
   IF (tpose == 1) THEN
      q = sv
      DO j = 1, 8
         DO i = 1, neq
            q(i) = q(i) + DOT_PRODUCT(kmat(:,i,j),psii(:,j))
         END DO
      END DO
   ELSE
      q = sv + MATMUL(kmat(:,:,1),psii(:,1)) + MATMUL(kmat(:,:,2),psii(:,2)) + MATMUL(kmat(:,:,3),psii(:,3)) &
             + MATMUL(kmat(:,:,4),psii(:,4)) + MATMUL(kmat(:,:,5),psii(:,5)) + MATMUL(kmat(:,:,6),psii(:,6)) &
             + MATMUL(kmat(:,:,7),psii(:,7)) + MATMUL(kmat(:,:,8),psii(:,8))
   END IF
   ! CG only works if theres is a non-zero RHS
   IF (MAXVAL(q) <= 0.0) THEN
      phi(:,g) = 0.0
   ELSE
      CALL cgsolve(g,q,its)
   END IF
END IF

! Compute the new angular flux values at the boundaries from jpsi and kpsi
q = phi(:,g) + src

! Set up the loop independent variables
ix1 = 1
ix2 = nx*ny
ix3 = ix2 + 1
ix4 = ix2 + nx*nz
ix5 = ix4 + 1
ix6 = ix4 + ny*nz
! Find psio depending on 'tpose' flag
IF (tpose == 1) THEN
   DO j = 1, 8
      DO i = 1, bcs
         psio(i,j,g) = DOT_PRODUCT(jpsi(:,i,j),q)
      END DO
      DO i = 1, apo
         kz = (i-1)*(nx*ny) + 1
         lz = i*(nx*ny)
         ky = apo*(nx*ny) + (i-1)*(nx*nz) + 1
         ly = apo*(nx*ny) + i*(nx*nz)
         kx = apo*(nx*ny+nx*nz) + (i-1)*(ny*nz) + 1
         lx = apo*(nx*ny+nx*nz) + i*(ny*nz)
         jndx = 0
         ! zbc's in psio
         DO k = kz, lz
            jndx = jndx + 1
            psio(k,j,g) = psio(k,j,g) + DOT_PRODUCT(kpsi(ix1:ix2,jndx,i,j),psii(kz:lz,j)) &
                                      + DOT_PRODUCT(kpsi(ix3:ix4,jndx,i,j),psii(ky:ly,j)) &
                                      + DOT_PRODUCT(kpsi(ix5:ix6,jndx,i,j),psii(kx:lx,j))
         END DO
         ! ybc's in psio
         DO k = ky, ly
            jndx = jndx + 1
            psio(k,j,g) = psio(k,j,g) + DOT_PRODUCT(kpsi(ix1:ix2,jndx,i,j),psii(kz:lz,j)) &
                                      + DOT_PRODUCT(kpsi(ix3:ix4,jndx,i,j),psii(ky:ly,j)) &
                                      + DOT_PRODUCT(kpsi(ix5:ix6,jndx,i,j),psii(kx:lx,j))
         END DO
         ! xbc's in psio
         DO k = kx, lx
            jndx = jndx + 1
            psio(k,j,g) = psio(k,j,g) + DOT_PRODUCT(kpsi(ix1:ix2,jndx,i,j),psii(kz:lz,j)) &
                                      + DOT_PRODUCT(kpsi(ix3:ix4,jndx,i,j),psii(ky:ly,j)) &
                                      + DOT_PRODUCT(kpsi(ix5:ix6,jndx,i,j),psii(kx:lx,j))
         END DO
      END DO
   END DO
ELSE
   DO j = 1, 8
      psio(:,j,g) = MATMUL(jpsi(:,:,j),q)
      DO i = 1, apo
         kz = (i-1)*(nx*ny) + 1
         lz = i*(nx*ny)
         ky = apo*(nx*ny) + (i-1)*(nx*nz) + 1
         ly = apo*(nx*ny) + i*(nx*nz)
         kx = apo*(nx*ny+nx*nz) + (i-1)*(ny*nz) + 1
         lx = apo*(nx*ny+nx*nz) + i*(ny*nz)
         ! zbc's in psio
         psio(kz:lz,j,g) = psio(kz:lz,j,g) + MATMUL(kpsi(ix1:ix2,ix1:ix2,i,j),psii(kz:lz,j)) &
                                           + MATMUL(kpsi(ix1:ix2,ix3:ix4,i,j),psii(ky:ly,j)) &
                                           + MATMUL(kpsi(ix1:ix2,ix5:ix6,i,j),psii(kx:lx,j))
         ! ybc's in psio
         psio(ky:ly,j,g) = psio(ky:ly,j,g) + MATMUL(kpsi(ix3:ix4,ix1:ix2,i,j),psii(kz:lz,j)) &
                                           + MATMUL(kpsi(ix3:ix4,ix3:ix4,i,j),psii(ky:ly,j)) &
                                           + MATMUL(kpsi(ix3:ix4,ix5:ix6,i,j),psii(kx:lx,j))
         ! xbc's in psio
         psio(kx:lx,j,g) = psio(kx:lx,j,g) + MATMUL(kpsi(ix5:ix6,ix1:ix2,i,j),psii(kz:lz,j)) &
                                           + MATMUL(kpsi(ix5:ix6,ix3:ix4,i,j),psii(ky:ly,j)) &
                                           + MATMUL(kpsi(ix5:ix6,ix5:ix6,i,j),psii(kx:lx,j))
      END DO
   END DO
END IF

RETURN
END SUBROUTINE idot
