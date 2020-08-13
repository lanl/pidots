SUBROUTINE residual(psii,oldpsi,rphi,rpsi,kmat,kpsi)

!-------------------------------------------------------------
!
!  Compute residual of the fine grid calculation
!
!-------------------------------------------------------------

USE invar
IMPLICIT NONE
INTEGER :: i, j, k, ix1, ix2, ix3, ix4, ix5, ix6, kz, lz, ky, ly, kx, lx, jndx, sdi
REAL*8, DIMENSION(neq,nsdp), INTENT(OUT) :: rphi
REAL*8, DIMENSION(bcs,8,nsdp) :: dpsi
REAL*8, DIMENSION(bcs,8,nsdp), INTENT(IN) :: psii, oldpsi
REAL*8, DIMENSION(bcs,8,nsdp), INTENT(OUT) :: rpsi
REAL*8, DIMENSION(bcs,neq,8,nsdp), INTENT(IN) :: kmat
REAL*8, DIMENSION(bcs2,bcs2,apo,8,nsdp), INTENT(IN) :: kpsi

! Get difference in psi-in iterates
dpsi = psii - oldpsi

! Compute the residuals
! phi
rphi = 0.0
DO sdi = 1, nsdp
   DO j = 1, 8
      DO i = 1, neq
         rphi(i,sdi) = rphi(i,sdi) + DOT_PRODUCT(kmat(:,i,j,sdi),dpsi(:,j,sdi))
      END DO
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
DO sdi = 1, nsdp
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
            rpsi(k,j,sdi) = rpsi(k,j,sdi) + DOT_PRODUCT(kpsi(ix1:ix2,jndx,i,j,sdi),dpsi(kz:lz,j,sdi)) &
                                          + DOT_PRODUCT(kpsi(ix3:ix4,jndx,i,j,sdi),dpsi(ky:ly,j,sdi)) &
                                          + DOT_PRODUCT(kpsi(ix5:ix6,jndx,i,j,sdi),dpsi(kx:lx,j,sdi))
         END DO
         DO k = ky, ly
            jndx = jndx + 1
            rpsi(k,j,sdi) = rpsi(k,j,sdi) + DOT_PRODUCT(kpsi(ix1:ix2,jndx,i,j,sdi),dpsi(kz:lz,j,sdi)) &
                                          + DOT_PRODUCT(kpsi(ix3:ix4,jndx,i,j,sdi),dpsi(ky:ly,j,sdi)) &
                                          + DOT_PRODUCT(kpsi(ix5:ix6,jndx,i,j,sdi),dpsi(kx:lx,j,sdi))
         END DO
         DO k = kx, lx
            jndx = jndx + 1
            rpsi(k,j,sdi) = rpsi(k,j,sdi) + DOT_PRODUCT(kpsi(ix1:ix2,jndx,i,j,sdi),dpsi(kz:lz,j,sdi)) &
                                          + DOT_PRODUCT(kpsi(ix3:ix4,jndx,i,j,sdi),dpsi(ky:ly,j,sdi)) &
                                          + DOT_PRODUCT(kpsi(ix5:ix6,jndx,i,j,sdi),dpsi(kx:lx,j,sdi))
         END DO
      END DO
   END DO
END DO

RETURN
END SUBROUTINE residual
