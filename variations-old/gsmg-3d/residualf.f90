SUBROUTINE residualf(fpi,foldpsi,frphi,frpsi,fkmt,fkps)

!-------------------------------------------------------------
!
!  Compute residual of the fine grid calculation
!
!-------------------------------------------------------------

USE invar
IMPLICIT NONE
INTEGER :: i, j, k, ix1, ix2, ix3, ix4, ix5, ix6, kz, lz, ky, ly, kx, lx, jndx
REAL*8, DIMENSION(neq), INTENT(OUT) :: frphi
REAL*8, DIMENSION(bcs,8), INTENT(IN) :: fpi, foldpsi
REAL*8, DIMENSION(bcs,8) :: dpsi
REAL*8, DIMENSION(bcs,8), INTENT(OUT) :: frpsi
REAL*8, DIMENSION(bcs,neq,8), INTENT(IN) :: fkmt
REAL*8, DIMENSION(bcs2,bcs2,apo,8), INTENT(IN) :: fkps

! Get difference in psi-in iterates
dpsi = fpi - foldpsi

! Compute the residuals
! phi
frphi = 0.0
DO j = 1, 8
   DO i = 1, neq
      frphi(i) = frphi(i) + DOT_PRODUCT(fkmt(:,i,j),dpsi(:,j))
   END DO
END DO

! psi
frpsi = 0.0
ix1 = 1
ix2 = xys
ix3 = ix2 + 1
ix4 = ix2 + xzs
ix5 = ix4 + 1
ix6 = ix4 + yzs
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
         frpsi(k,j) = frpsi(k,j) + DOT_PRODUCT(fkps(ix1:ix2,jndx,i,j),dpsi(kz:lz,j)) &
                                 + DOT_PRODUCT(fkps(ix3:ix4,jndx,i,j),dpsi(ky:ly,j)) &
                                 + DOT_PRODUCT(fkps(ix5:ix6,jndx,i,j),dpsi(kx:lx,j))
      END DO
      DO k = ky, ly
         jndx = jndx + 1
         frpsi(k,j) = frpsi(k,j) + DOT_PRODUCT(fkps(ix1:ix2,jndx,i,j),dpsi(kz:lz,j)) &
                                 + DOT_PRODUCT(fkps(ix3:ix4,jndx,i,j),dpsi(ky:ly,j)) &
                                 + DOT_PRODUCT(fkps(ix5:ix6,jndx,i,j),dpsi(kx:lx,j))
      END DO
      DO k = kx, lx
         jndx = jndx + 1
         frpsi(k,j) = frpsi(k,j) + DOT_PRODUCT(fkps(ix1:ix2,jndx,i,j),dpsi(kz:lz,j)) &
                                 + DOT_PRODUCT(fkps(ix3:ix4,jndx,i,j),dpsi(ky:ly,j)) &
                                 + DOT_PRODUCT(fkps(ix5:ix6,jndx,i,j),dpsi(kx:lx,j))
      END DO
   END DO
END DO

RETURN
END SUBROUTINE residualf
