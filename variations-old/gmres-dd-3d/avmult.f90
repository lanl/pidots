SUBROUTINE avmult(sv,ov,iv,ws,wa)

!-------------------------------------------------------------
!
!  GMRES: Av(:,j)
!  Matrix vector multiplication performed with ITMM operators.
!  Follows same steps as IDOT from the block Jacobi code.
!
!  Also serves the purpose of Ax (ITMM operators multiplied by
!   solution vectors)
!
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER :: i, j, k, ix1, ix2, ix3, ix4, ix5, ix6, kz, lz, ky, ly, kx, lx, jndx
REAL*8, DIMENSION(neq), INTENT(IN) :: sv
REAL*8, DIMENSION(neq), INTENT(OUT) :: ws
REAL*8, DIMENSION(bcs,8), INTENT(IN) :: ov, iv
REAL*8, DIMENSION(bcs,8), INTENT(OUT) :: wa

! Initialize
ws = 0.0
wa = 0.0

! Compute the scalar flux portion, ws: jmat, kmat multiplications
IF (tpose == 1) THEN
   ! jmat
   DO i = 1, neq
      ws(i) = DOT_PRODUCT(jmat(:,i),sv)
   END DO
   ! kmat
   ! Subtraction necessary because RHS terms moved to LHS of equation
   DO j = 1, 8
      DO i = 1, neq
         ws(i) = ws(i) - DOT_PRODUCT(kmat(:,i,j),iv(:,j))
      END DO
   END DO
ELSE
   ! jmat
   ws = MATMUL(jmat,sv)
   ! kmat
   ws = ws - (MATMUL(kmat(:,:,1),iv(:,1)) + MATMUL(kmat(:,:,2),iv(:,2)) + MATMUL(kmat(:,:,3),iv(:,3)) &
            + MATMUL(kmat(:,:,4),iv(:,4)) + MATMUL(kmat(:,:,5),iv(:,5)) + MATMUL(kmat(:,:,6),iv(:,6)) &
            + MATMUL(kmat(:,:,7),iv(:,7)) + MATMUL(kmat(:,:,8),iv(:,8)))
END IF

! Compute the angular flux portion, wa: jpsi, kpsi multiplications
! Set up the loop independent variables
ix1 = 1
ix2 = xys
ix3 = ix2 + 1
ix4 = ix2 + xzs
ix5 = ix4 + 1
ix6 = ix4 + yzs
! Find wa depending on 'tpose' flag
IF (tpose == 1) THEN
   DO j = 1, 8
      DO i = 1, bcs
         wa(i,j) = DOT_PRODUCT(jpsi(:,i,j),sv)
      END DO
      DO i = 1, apo
         kz = (i-1)*xys + 1
         lz = i*xys
         ky = apo*xys + (i-1)*xzs + 1
         ly = apo*xys + i*xzs
         kx = apo*(xys+xzs) + (i-1)*yzs + 1
         lx = apo*(xys+xzs) + i*yzs
         jndx = 0
         ! zbc's in wa
         DO k = kz, lz
            jndx = jndx + 1
            wa(k,j) = wa(k,j) + DOT_PRODUCT(kpsi(ix1:ix2,jndx,i,j),iv(kz:lz,j)) &
                              + DOT_PRODUCT(kpsi(ix3:ix4,jndx,i,j),iv(ky:ly,j)) &
                              + DOT_PRODUCT(kpsi(ix5:ix6,jndx,i,j),iv(kx:lx,j))
         END DO
         ! ybc's in wa
         DO k = ky, ly
            jndx = jndx + 1
            wa(k,j) = wa(k,j) + DOT_PRODUCT(kpsi(ix1:ix2,jndx,i,j),iv(kz:lz,j)) &
                              + DOT_PRODUCT(kpsi(ix3:ix4,jndx,i,j),iv(ky:ly,j)) &
                              + DOT_PRODUCT(kpsi(ix5:ix6,jndx,i,j),iv(kx:lx,j))
         END DO
         ! xbc's in wa
         DO k = kx, lx
            jndx = jndx + 1
            wa(k,j) = wa(k,j) + DOT_PRODUCT(kpsi(ix1:ix2,jndx,i,j),iv(kz:lz,j)) &
                              + DOT_PRODUCT(kpsi(ix3:ix4,jndx,i,j),iv(ky:ly,j)) &
                              + DOT_PRODUCT(kpsi(ix5:ix6,jndx,i,j),iv(kx:lx,j))
         END DO
      END DO
   END DO
ELSE
   DO j = 1, 8
      wa(:,j) = MATMUL(jpsi(:,:,j),sv)
      DO i = 1, apo
         kz = (i-1)*xys + 1
         lz = i*xys
         ky = apo*xys + (i-1)*xzs + 1
         ly = apo*xys + i*xzs
         kx = apo*(xys+xzs) + (i-1)*yzs + 1
         lx = apo*(xys+xzs) + i*yzs
         ! zbc's in wa
         wa(kz:lz,j) = wa(kz:lz,j) + MATMUL(kpsi(ix1:ix2,ix1:ix2,i,j),iv(kz:lz,j)) &
                                   + MATMUL(kpsi(ix1:ix2,ix3:ix4,i,j),iv(ky:ly,j)) &
                                   + MATMUL(kpsi(ix1:ix2,ix5:ix6,i,j),iv(kx:lx,j))
         ! ybc's in wa
         wa(ky:ly,j) = wa(ky:ly,j) + MATMUL(kpsi(ix3:ix4,ix1:ix2,i,j),iv(kz:lz,j)) &
                                   + MATMUL(kpsi(ix3:ix4,ix3:ix4,i,j),iv(ky:ly,j)) &
                                   + MATMUL(kpsi(ix3:ix4,ix5:ix6,i,j),iv(kx:lx,j))
         ! xbc's in wa
         wa(kx:lx,j) = wa(kx:lx,j) + MATMUL(kpsi(ix5:ix6,ix1:ix2,i,j),iv(kz:lz,j)) &
                                   + MATMUL(kpsi(ix5:ix6,ix3:ix4,i,j),iv(ky:ly,j)) &
                                   + MATMUL(kpsi(ix5:ix6,ix5:ix6,i,j),iv(kx:lx,j))
      END DO
   END DO
END IF

! When moving terms to LHS, everything subtracted from solution vector (1 on
! diagonal in the psio equations). Must be reflected in wa.
wa = ov - wa

RETURN
END SUBROUTINE avmult
