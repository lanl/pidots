SUBROUTINE idot(g)

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
INTEGER :: i, j, ky, ly, kx, lx, ix1, ix2, ix3, ix4
REAL*8, DIMENSION(ordsq*nx*ny) :: q

! Solve for the new phi with the sum of sv with the product of kmat and psi
phi(:,g) = sv + MATMUL(kmat(:,:,1),psi1) + MATMUL(kmat(:,:,2),psi2) &
              + MATMUL(kmat(:,:,3),psi3) + MATMUL(kmat(:,:,4),psi4)
cnvf(g) = 1

! Compute the new angular flux values at the boundaries from jpsi and kpsi
q = phi(:,g) + src

! Set up the loop independent indices
ix1 = 1
ix2 = nx*order
ix3 = ix2 + 1
ix4 = ix2 + ny*order
DO j = 1, 4
   psio(:,j,g) = MATMUL(jpsi(:,:,j),q)
   DO i = 1, apo
      ky = (i-1)*nx*order + 1
      ly = i*nx*order
      kx = apo*nx*order + (i-1)*ny*order + 1
      lx = apo*nx*order + i*ny*order
      IF (j == 1) THEN
         psio(ky:ly,j,g) = psio(ky:ly,j,g) + MATMUL(kpsi(ix1:ix2,ix1:ix2,i,j),psi1(ky:ly)) &
                                           + MATMUL(kpsi(ix1:ix2,ix3:ix4,i,j),psi1(kx:lx))
         psio(kx:lx,j,g) = psio(kx:lx,j,g) + MATMUL(kpsi(ix3:ix4,ix1:ix2,i,j),psi1(ky:ly)) &
                                           + MATMUL(kpsi(ix3:ix4,ix3:ix4,i,j),psi1(kx:lx))
      ELSE IF (j == 2) THEN
         psio(ky:ly,j,g) = psio(ky:ly,j,g) + MATMUL(kpsi(ix1:ix2,ix1:ix2,i,j),psi2(ky:ly)) &
                                           + MATMUL(kpsi(ix1:ix2,ix3:ix4,i,j),psi2(kx:lx))
         psio(kx:lx,j,g) = psio(kx:lx,j,g) + MATMUL(kpsi(ix3:ix4,ix1:ix2,i,j),psi2(ky:ly)) &
                                           + MATMUL(kpsi(ix3:ix4,ix3:ix4,i,j),psi2(kx:lx))
      ELSE IF (j == 3) THEN
         psio(ky:ly,j,g) = psio(ky:ly,j,g) + MATMUL(kpsi(ix1:ix2,ix1:ix2,i,j),psi3(ky:ly)) &
                                           + MATMUL(kpsi(ix1:ix2,ix3:ix4,i,j),psi3(kx:lx))
         psio(kx:lx,j,g) = psio(kx:lx,j,g) + MATMUL(kpsi(ix3:ix4,ix1:ix2,i,j),psi3(ky:ly)) &
                                           + MATMUL(kpsi(ix3:ix4,ix3:ix4,i,j),psi3(kx:lx))
      ELSE
         psio(ky:ly,j,g) = psio(ky:ly,j,g) + MATMUL(kpsi(ix1:ix2,ix1:ix2,i,j),psi4(ky:ly)) &
                                           + MATMUL(kpsi(ix1:ix2,ix3:ix4,i,j),psi4(kx:lx))
         psio(kx:lx,j,g) = psio(kx:lx,j,g) + MATMUL(kpsi(ix3:ix4,ix1:ix2,i,j),psi4(ky:ly)) &
                                           + MATMUL(kpsi(ix3:ix4,ix3:ix4,i,j),psi4(kx:lx))
      END IF
   END DO
END DO

RETURN
END SUBROUTINE idot
