SUBROUTINE idot(v,bit)

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
INTEGER, INTENT(IN) :: v, bit
INTEGER :: i, j, k, ix1, ix2, ix3, ix4, ix5, ix6, kz, lz, ky, ly, kx, lx, jndx, info

! Compute the solution directly
! Always use sv and product of kmat and psii to form RHS
phi(:,v) = sv(:,v)
DO j = 1, 8
   DO i = 1, neq
      phi(i,v) = phi(i,v) + DOT_PRODUCT(kmat(:,i,j,v),psii(:,j,v))
   END DO
END DO

! Use LAPACK solvers if matrix was factored only
CALL dgetrs('T',neq,1,jmat(:,:,v),neq,piv(:,v),phi(:,v),neq,info)
IF (info /= 0) THEN
   WRITE (8,'(//,1X,A)') "ERROR: dgetrs cannot solve."
   STOP
END IF

! Compute the new angular flux values at the boundaries from jpsi and kpsi
! Set up the loop independent variables
ix1 = 1
ix2 = xys
ix3 = ix2 + 1
ix4 = ix2 + xzs
ix5 = ix4 + 1
ix6 = ix4 + yzs
! Solve for psio
DO j = 1, 8
   DO i = 1, bcs
      psio(i,j,v) = DOT_PRODUCT(jpsi(:,i,j,v),phi(:,v))
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
         psio(k,j,v) = psio(k,j,v) + DOT_PRODUCT(kpsi(ix1:ix2,jndx,i,j,v),psii(kz:lz,j,v)) &
                                   + DOT_PRODUCT(kpsi(ix3:ix4,jndx,i,j,v),psii(ky:ly,j,v)) &
                                   + DOT_PRODUCT(kpsi(ix5:ix6,jndx,i,j,v),psii(kx:lx,j,v))
      END DO
      ! ybc's in psio
      DO k = ky, ly
         jndx = jndx + 1
         psio(k,j,v) = psio(k,j,v) + DOT_PRODUCT(kpsi(ix1:ix2,jndx,i,j,v),psii(kz:lz,j,v)) &
                                   + DOT_PRODUCT(kpsi(ix3:ix4,jndx,i,j,v),psii(ky:ly,j,v)) &
                                   + DOT_PRODUCT(kpsi(ix5:ix6,jndx,i,j,v),psii(kx:lx,j,v))
      END DO
      ! xbc's in psio
      DO k = kx, lx
         jndx = jndx + 1
         psio(k,j,v) = psio(k,j,v) + DOT_PRODUCT(kpsi(ix1:ix2,jndx,i,j,v),psii(kz:lz,j,v)) &
                                   + DOT_PRODUCT(kpsi(ix3:ix4,jndx,i,j,v),psii(ky:ly,j,v)) &
                                   + DOT_PRODUCT(kpsi(ix5:ix6,jndx,i,j,v),psii(kx:lx,j,v))
      END DO
   END DO
END DO

! Add in RHS vector
psio(:,:,v) = psio(:,:,v) + av(:,:,v)

RETURN
END SUBROUTINE idot
