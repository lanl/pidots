SUBROUTINE idot(sdi,bit)

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
INTEGER, INTENT(IN) :: sdi, bit
INTEGER :: i, j, k, ix1, ix2, ix3, ix4, ix5, ix6, kz, lz, ky, ly, kx, lx, jndx, info
REAL*8, DIMENSION(neq) :: q

! Always compute the solution directly
! Always use sv and product of kmat and psii to either solve for phi, or form RHS
phi(:,sdi) = sv(:,sdi)
DO j = 1, 8
   DO i = 1, neq
      phi(i,sdi) = phi(i,sdi) + DOT_PRODUCT(kmat(:,i,j,sdi),psii(:,j,sdi))
   END DO
END DO

CALL dgetrs('T',neq,1,jmat(:,:,sdi),neq,piv(:,sdi),phi(:,sdi),neq,info)
IF (info /= 0) THEN
   WRITE (8,'(//,1X,A)') "ERROR: dgetrs cannot solve."
   STOP
END IF

! Compute the new angular flux values at the boundaries from jpsi and kpsi
q = phi(:,sdi) + src(:,sdi)

! Set up the loop independent variables
ix1 = 1
ix2 = xys
ix3 = ix2 + 1
ix4 = ix2 + xzs
ix5 = ix4 + 1
ix6 = ix4 + yzs
DO j = 1, 8
   DO i = 1, bcs
      psio(i,j,sdi) = DOT_PRODUCT(jpsi(:,i,j,sdi),q)
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
         psio(k,j,sdi) = psio(k,j,sdi) + DOT_PRODUCT(kpsi(ix1:ix2,jndx,i,j,sdi),psii(kz:lz,j,sdi)) &
                                       + DOT_PRODUCT(kpsi(ix3:ix4,jndx,i,j,sdi),psii(ky:ly,j,sdi)) &
                                       + DOT_PRODUCT(kpsi(ix5:ix6,jndx,i,j,sdi),psii(kx:lx,j,sdi))
      END DO
      ! ybc's in psio
      DO k = ky, ly
         jndx = jndx + 1
         psio(k,j,sdi) = psio(k,j,sdi) + DOT_PRODUCT(kpsi(ix1:ix2,jndx,i,j,sdi),psii(kz:lz,j,sdi)) &
                                       + DOT_PRODUCT(kpsi(ix3:ix4,jndx,i,j,sdi),psii(ky:ly,j,sdi)) &
                                       + DOT_PRODUCT(kpsi(ix5:ix6,jndx,i,j,sdi),psii(kx:lx,j,sdi))
      END DO
      ! xbc's in psio
      DO k = kx, lx
         jndx = jndx + 1
         psio(k,j,sdi) = psio(k,j,sdi) + DOT_PRODUCT(kpsi(ix1:ix2,jndx,i,j,sdi),psii(kz:lz,j,sdi)) &
                                       + DOT_PRODUCT(kpsi(ix3:ix4,jndx,i,j,sdi),psii(ky:ly,j,sdi)) &
                                       + DOT_PRODUCT(kpsi(ix5:ix6,jndx,i,j,sdi),psii(kx:lx,j,sdi))
      END DO
   END DO
END DO

RETURN
END SUBROUTINE idot
