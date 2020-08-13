SUBROUTINE idot(g,its,bit)

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
INTEGER, INTENT(IN) :: g, bit
INTEGER :: i, j, k, ix1, ix2, ix3, ix4, ix5, ix6, kz, lz, ky, ly, kx, lx, jndx, info
INTEGER, INTENT(OUT) :: its
REAL*8, DIMENSION(neq) :: q

! Compute the solution directly or with CG iterations
IF (idos == 0) THEN
   ! Always use sv and product of kmat and psii to either solve for phi, or form RHS
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

   ! Use LAPACK solvers if matrix was factored only
   IF (invf == 0) THEN
      IF (tpose == 1) THEN
         IF (sym == 1) THEN
            CALL dpotrs('U',neq,1,jmat,neq,phi(:,g),neq,info)
            IF (info /= 0) THEN
               WRITE (8,'(//,1X,A)') "ERROR: dpotrs cannot solve."
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
            CALL dpotrs('U',neq,1,jmat,neq,phi(:,g),neq,info)
            IF (info /= 0) THEN
               WRITE (8,'(//,1X,A)') "ERROR: dpotrs cannot solve."
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
   END IF
   cnvf(g) = 1
   its = 0
ELSE
   ! Check the number of P
   IF (isize == 1) THEN
      WRITE (8,'(/,1X,A,/)') "WARNING: Depending on BCs, only one global iteration may be performed, &
                              and results will be based on local convergence criterion, not global."
      warn = warn + 1
   END IF
   ! Solve for the new phi with CG iterations from the cgsolve routine
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
   ! SOR and CG only work if there's is a non-zero RHS
   IF (MAXVAL(q) <= 0.0) THEN
      phi(:,g) = 0.0
      cnvf(g) = 0
   ELSE
      IF (idos == 1) THEN
         IF (pcf == 0) THEN
            CALL cgsolve(g,q,its)
         ELSE
            IF (pcty == 0) THEN
               CALL cgsolvejac(g,q,its)
            ELSE IF (pcty == 1) THEN
               CALL cgsolvesgs(g,q,its)
            ELSE
               CALL cgsolvesor(g,q,its)
            END IF
         END IF
      ELSE
         CALL sor(g,q,its)
      END IF
   END IF
END IF

! Compute the new angular flux values at the boundaries from jpsi and kpsi
q = phi(:,g)

! Set up the loop independent variables
ix1 = 1
ix2 = xys
ix3 = ix2 + 1
ix4 = ix2 + xzs
ix5 = ix4 + 1
ix6 = ix4 + yzs
! Find psio depending on 'tpose' flag
IF (tpose == 1) THEN
   DO j = 1, 8
      DO i = 1, bcs
         psio(i,j,g) = DOT_PRODUCT(jpsi(:,i,j),q)
      END DO
      DO i = 1, apo
         kz = (i-1)*xys + 1
         lz = i*xys
         ky = apo*xys + (i-1)*xzs + 1
         ly = apo*xys + i*xzs
         kx = apo*(xys+xzs) + (i-1)*yzs + 1
         lx = apo*(xys+xzs) + i*yzs
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
         kz = (i-1)*xys + 1
         lz = i*xys
         ky = apo*xys + (i-1)*xzs + 1
         ly = apo*xys + i*xzs
         kx = apo*(xys+xzs) + (i-1)*yzs + 1
         lx = apo*(xys+xzs) + i*yzs
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
! Lastly add in the RHS
psio(:,:,g) = psio(:,:,g) + rpsi(:,:)

RETURN
END SUBROUTINE idot
