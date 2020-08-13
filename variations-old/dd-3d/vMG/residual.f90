SUBROUTINE residual(g,rphi,rpsi)

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
INTEGER :: i, j, k, ix1, ix2, ix3, ix4, ix5, ix6, kz, lz, ky, ly, kx, lx, jndx
REAL*8, DIMENSION(neq) :: bphi, lphi
REAL*8, DIMENSION(bcs,8) :: bpsi, lpsi
REAL*8, DIMENSION(neq), INTENT(OUT) :: rphi
REAL*8, DIMENSION(bcs,8), INTENT(OUT) :: rpsi

! Compute the RHS components
! phi
bphi = sv
! psi
IF (tpose == 1) THEN
   DO j = 1, 8
      DO i = 1, bcs
         bpsi(i,j) = DOT_PRODUCT(jpsi(:,i,j),src)
      END DO
   END DO
ELSE
   DO j = 1, 8
      bpsi(:,j) = MATMUL(jpsi(:,:,j),src)
   END DO
END IF


! Compute the LHS components
! phi
IF (tpose == 1) THEN
   DO i = 1, neq
      lphi(i) = DOT_PRODUCT(jmat2(:,i),phi(:,g))
   END DO
   DO j = 1, 8
      DO i = 1, neq
         lphi(i) = lphi(i) - DOT_PRODUCT(kmat(:,i,j),psii(:,j))
      END DO
   END DO
ELSE
   lphi = MATMUL(jmat,phi(:,g))
   DO j = 1, 8
      lphi = lphi - MATMUL(kmat(:,:,j),psii(:,j))
   END DO
END IF

! psi
ix1 = 1
ix2 = xys
ix3 = ix2 + 1
ix4 = ix2 + xzs
ix5 = ix4 + 1
ix6 = ix4 + yzs
IF (tpose == 1) THEN
   DO j = 1, 8
      DO i = 1, bcs
         lpsi(i,j) = DOT_PRODUCT(jpsi(:,i,j),phi(:,g))
      END DO
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
            lpsi(k,j) = lpsi(k,j) + DOT_PRODUCT(kpsi(ix1:ix2,jndx,i,j),psii(kz:lz,j)) &
                                  + DOT_PRODUCT(kpsi(ix3:ix4,jndx,i,j),psii(ky:ly,j)) &
                                  + DOT_PRODUCT(kpsi(ix5:ix6,jndx,i,j),psii(kx:lx,j))
         END DO
         DO k = ky, ly
            jndx = jndx + 1
            lpsi(k,j) = lpsi(k,j) + DOT_PRODUCT(kpsi(ix1:ix2,jndx,i,j),psii(kz:lz,j)) &
                                  + DOT_PRODUCT(kpsi(ix3:ix4,jndx,i,j),psii(ky:ly,j)) &
                                  + DOT_PRODUCT(kpsi(ix5:ix6,jndx,i,j),psii(kx:lx,j))
         END DO
         DO k = kx, lx
            jndx = jndx + 1
            lpsi(k,j) = lpsi(k,j) + DOT_PRODUCT(kpsi(ix1:ix2,jndx,i,j),psii(kz:lz,j)) &
                                  + DOT_PRODUCT(kpsi(ix3:ix4,jndx,i,j),psii(ky:ly,j)) &
                                  + DOT_PRODUCT(kpsi(ix5:ix6,jndx,i,j),psii(kx:lx,j))
         END DO
      END DO
   END DO
ELSE
   DO j = 1, 8
      lpsi(:,j) = MATMUL(jpsi(:,:,j),phi(:,g))
      DO i = 1, apo
         kz = (i-1)*xys + 1
         lz = i*xys
         ky = apo*xys + (i-1)*xzs + 1
         ly = apo*xys + i*xzs
         kx = apo*(xys+xzs) + (i-1)*yzs + 1
         lx = apo*(xys+xzs) + i*yzs
         lpsi(kz:lz,j) = lpsi(kz:lz,j) + MATMUL(kpsi(ix1:ix2,ix1:ix2,i,j),psii(kz:lz,j)) &
                                       + MATMUL(kpsi(ix1:ix2,ix3:ix4,i,j),psii(ky:ly,j)) &
                                       + MATMUL(kpsi(ix1:ix2,ix5:ix6,i,j),psii(kx:lx,j))
         lpsi(ky:ly,j) = lpsi(ky:ly,j) + MATMUL(kpsi(ix3:ix4,ix1:ix2,i,j),psii(kz:lz,j)) &
                                       + MATMUL(kpsi(ix3:ix4,ix3:ix4,i,j),psii(ky:ly,j)) &
                                       + MATMUL(kpsi(ix3:ix4,ix5:ix6,i,j),psii(kx:lx,j))
         lpsi(kx:lx,j) = lpsi(kx:lx,j) + MATMUL(kpsi(ix5:ix6,ix1:ix2,i,j),psii(kz:lz,j)) &
                                       + MATMUL(kpsi(ix5:ix6,ix3:ix4,i,j),psii(ky:ly,j)) &
                                       + MATMUL(kpsi(ix5:ix6,ix5:ix6,i,j),psii(kx:lx,j))
      END DO
   END DO
END IF
lpsi = psio(:,:,g) - lpsi

! Compute the residual
rphi = bphi - lphi
rpsi = bpsi - lpsi

RETURN
END SUBROUTINE residual
