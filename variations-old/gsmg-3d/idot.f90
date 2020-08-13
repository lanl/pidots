SUBROUTINE idot(phi,sv,psii,psio,av,jmat,kmat,jpsi,kpsi,piv)

!-------------------------------------------------------------
!
!  Integral discrete ordinates transport
!
!  Add sv to kmat*psii to get the new phi
!  Compute psio with jpsi*(phi+src) + kpsi*psii
!
!-------------------------------------------------------------

USE invar
IMPLICIT NONE
INTEGER :: i, j, k, ix1, ix2, ix3, ix4, ix5, ix6, kz, lz, ky, ly, kx, lx, jndx, info
INTEGER, DIMENSION(neq), INTENT(IN) :: piv
REAL*8, DIMENSION(neq), INTENT(OUT) :: phi
REAL*8, DIMENSION(neq), INTENT(IN) :: sv
REAL*8, DIMENSION(bcs,8), INTENT(OUT) :: psio
REAL*8, DIMENSION(bcs,8), INTENT(IN) :: psii, av
REAL*8, DIMENSION(neq,neq), INTENT(IN) :: jmat
REAL*8, DIMENSION(bcs,neq,8), INTENT(IN) :: kmat
REAL*8, DIMENSION(neq,bcs,8), INTENT(IN) :: jpsi
REAL*8, DIMENSION(bcs2,bcs2,apo,8), INTENT(IN) :: kpsi

! Always compute the solution directly
! Always use sv and product of kmat and psii to either solve for phi, or form RHS
phi = sv
DO j = 1, 8
   DO i = 1, neq
      phi(i) = phi(i) + DOT_PRODUCT(kmat(:,i,j),psii(:,j))
   END DO
END DO

CALL dgetrs('T',neq,1,jmat,neq,piv,phi,neq,info)
IF (info /= 0) THEN
   WRITE (8,'(//,1X,A)') "ERROR: dgetrs cannot solve."
   STOP
END IF

! Set up the loop independent variables
ix1 = 1
ix2 = xys
ix3 = ix2 + 1
ix4 = ix2 + xzs
ix5 = ix4 + 1
ix6 = ix4 + yzs
DO j = 1, 8
   DO i = 1, bcs
      psio(i,j) = DOT_PRODUCT(jpsi(:,i,j),phi)
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
         psio(k,j) = psio(k,j) + DOT_PRODUCT(kpsi(ix1:ix2,jndx,i,j),psii(kz:lz,j)) &
                               + DOT_PRODUCT(kpsi(ix3:ix4,jndx,i,j),psii(ky:ly,j)) &
                               + DOT_PRODUCT(kpsi(ix5:ix6,jndx,i,j),psii(kx:lx,j))
      END DO
      ! ybc's in psio
      DO k = ky, ly
         jndx = jndx + 1
         psio(k,j) = psio(k,j) + DOT_PRODUCT(kpsi(ix1:ix2,jndx,i,j),psii(kz:lz,j)) &
                               + DOT_PRODUCT(kpsi(ix3:ix4,jndx,i,j),psii(ky:ly,j)) &
                               + DOT_PRODUCT(kpsi(ix5:ix6,jndx,i,j),psii(kx:lx,j))
      END DO
      ! xbc's in psio
      DO k = kx, lx
         jndx = jndx + 1
         psio(k,j) = psio(k,j) + DOT_PRODUCT(kpsi(ix1:ix2,jndx,i,j),psii(kz:lz,j)) &
                               + DOT_PRODUCT(kpsi(ix3:ix4,jndx,i,j),psii(ky:ly,j)) &
                               + DOT_PRODUCT(kpsi(ix5:ix6,jndx,i,j),psii(kx:lx,j))
      END DO
   END DO
END DO

! Add in RHS vector
psio = psio + av

RETURN
END SUBROUTINE idot
