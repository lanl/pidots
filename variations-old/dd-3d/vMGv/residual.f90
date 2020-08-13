SUBROUTINE residual(v,rphi,rpsi)

!-------------------------------------------------------------
!
!  Compute residual of the fine grid calculation
!
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER, INTENT(IN) :: v
INTEGER :: i, j, k, ix1, ix2, ix3, ix4, ix5, ix6, kz, lz, ky, ly, kx, lx, jndx
REAL*8, DIMENSION(neq) :: dphi
REAL*8, DIMENSION(neq), INTENT(OUT) :: rphi
REAL*8, DIMENSION(bcs,8) :: dpsi
REAL*8, DIMENSION(bcs,8), INTENT(OUT) :: rpsi

! Get difference in psi-in iterates
dpsi = psii(:,:,v) - oldpsi
dphi = phi(:,v) - oldphi

! Compute the residuals
! phi
rphi = 0.0
DO j = 1, 8
   DO i = 1, neq
      rphi(i) = rphi(i) + DOT_PRODUCT(kmat(:,i,j,v),dpsi(:,j))
   END DO
END DO

! psi
rpsi = 0.0
ix1 = 1
ix2 = xys
ix3 = ix2 + 1
ix4 = ix2 + xzs
ix5 = ix4 + 1
ix6 = ix4 + yzs
DO j = 1, 8
   DO i = 1, bcs
!      rpsi(i,j) = rpsi(i,j) + DOT_PRODUCT(jpsi(:,i,j,v),dphi)
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
         rpsi(k,j) = rpsi(k,j) + DOT_PRODUCT(kpsi(ix1:ix2,jndx,i,j,v),dpsi(kz:lz,j)) &
                               + DOT_PRODUCT(kpsi(ix3:ix4,jndx,i,j,v),dpsi(ky:ly,j)) &
                               + DOT_PRODUCT(kpsi(ix5:ix6,jndx,i,j,v),dpsi(kx:lx,j))
      END DO
      DO k = ky, ly
         jndx = jndx + 1
         rpsi(k,j) = rpsi(k,j) + DOT_PRODUCT(kpsi(ix1:ix2,jndx,i,j,v),dpsi(kz:lz,j)) &
                               + DOT_PRODUCT(kpsi(ix3:ix4,jndx,i,j,v),dpsi(ky:ly,j)) &
                               + DOT_PRODUCT(kpsi(ix5:ix6,jndx,i,j,v),dpsi(kx:lx,j))
      END DO
      DO k = kx, lx
         jndx = jndx + 1
         rpsi(k,j) = rpsi(k,j) + DOT_PRODUCT(kpsi(ix1:ix2,jndx,i,j,v),dpsi(kz:lz,j)) &
                               + DOT_PRODUCT(kpsi(ix3:ix4,jndx,i,j,v),dpsi(ky:ly,j)) &
                               + DOT_PRODUCT(kpsi(ix5:ix6,jndx,i,j,v),dpsi(kx:lx,j))
      END DO
   END DO
END DO

RETURN
END SUBROUTINE residual
