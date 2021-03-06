SUBROUTINE residual(rphi,rpsi)

!-------------------------------------------------------------
!
!  Compute residual of the fine grid calculation
!
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER :: i, j, k, ix1, ix2, ix3, ix4, ix5, ix6, kz, lz, ky, ly, kx, lx, jndx
REAL*8, DIMENSION(bcs,8) :: dpsi
REAL*8, DIMENSION(neq), INTENT(OUT) :: rphi
REAL*8, DIMENSION(bcs,8), INTENT(OUT) :: rpsi
real*8, dimension(neq) :: bphi, lphi

! Get difference in psi-in iterates
dpsi = psii - oldpsi

! Compute the residuals
! phi
rphi = 0.0
IF (tpose == 1) THEN
   DO j = 1, 8
      DO i = 1, neq
         rphi(i) = rphi(i) + DOT_PRODUCT(kmat(:,i,j),dpsi(:,j))
      END DO
   END DO
ELSE
   DO j = 1, 8
      rphi = rphi + MATMUL(kmat(:,:,j),dpsi(:,j))
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
            rpsi(k,j) = rpsi(k,j) + DOT_PRODUCT(kpsi(ix1:ix2,jndx,i,j),dpsi(kz:lz,j)) &
                                  + DOT_PRODUCT(kpsi(ix3:ix4,jndx,i,j),dpsi(ky:ly,j)) &
                                  + DOT_PRODUCT(kpsi(ix5:ix6,jndx,i,j),dpsi(kx:lx,j))
         END DO
         DO k = ky, ly
            jndx = jndx + 1
            rpsi(k,j) = rpsi(k,j) + DOT_PRODUCT(kpsi(ix1:ix2,jndx,i,j),dpsi(kz:lz,j)) &
                                  + DOT_PRODUCT(kpsi(ix3:ix4,jndx,i,j),dpsi(ky:ly,j)) &
                                  + DOT_PRODUCT(kpsi(ix5:ix6,jndx,i,j),dpsi(kx:lx,j))
         END DO
         DO k = kx, lx
            jndx = jndx + 1
            rpsi(k,j) = rpsi(k,j) + DOT_PRODUCT(kpsi(ix1:ix2,jndx,i,j),dpsi(kz:lz,j)) &
                                  + DOT_PRODUCT(kpsi(ix3:ix4,jndx,i,j),dpsi(ky:ly,j)) &
                                  + DOT_PRODUCT(kpsi(ix5:ix6,jndx,i,j),dpsi(kx:lx,j))
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
         rpsi(kz:lz,j) = rpsi(kz:lz,j) + MATMUL(kpsi(ix1:ix2,ix1:ix2,i,j),dpsi(kz:lz,j)) &
                                       + MATMUL(kpsi(ix1:ix2,ix3:ix4,i,j),dpsi(ky:ly,j)) &
                                       + MATMUL(kpsi(ix1:ix2,ix5:ix6,i,j),dpsi(kx:lx,j))
         rpsi(ky:ly,j) = rpsi(ky:ly,j) + MATMUL(kpsi(ix3:ix4,ix1:ix2,i,j),dpsi(kz:lz,j)) &
                                       + MATMUL(kpsi(ix3:ix4,ix3:ix4,i,j),dpsi(ky:ly,j)) &
                                       + MATMUL(kpsi(ix3:ix4,ix5:ix6,i,j),dpsi(kx:lx,j))
         rpsi(kx:lx,j) = rpsi(kx:lx,j) + MATMUL(kpsi(ix5:ix6,ix1:ix2,i,j),dpsi(kz:lz,j)) &
                                       + MATMUL(kpsi(ix5:ix6,ix3:ix4,i,j),dpsi(ky:ly,j)) &
                                       + MATMUL(kpsi(ix5:ix6,ix5:ix6,i,j),dpsi(kx:lx,j))
      END DO
   END DO
END IF

RETURN
END SUBROUTINE residual
