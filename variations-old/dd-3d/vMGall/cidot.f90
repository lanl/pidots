SUBROUTINE cidot(rsv)

!-------------------------------------------------------------
!
!  Coarse Integral discrete ordinates transport
!
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER :: i, j, k, ix1, ix2, ix3, ix4, ix5, ix6, kz, lz, ky, ly, kx, lx, jndx, info
REAL*8, DIMENSION(cneq) :: q
REAL*8, DIMENSION(cneq), INTENT(IN) :: rsv

! Compute the solution directly
IF (tpose == 1) THEN
   f = rsv
   DO j = 1, 8
      DO i = 1, cneq
         f(i) = f(i) + DOT_PRODUCT(ckmat(:,i,j),cpsii(:,j))
      END DO
   END DO
ELSE
   f = rsv + MATMUL(ckmat(:,:,1),cpsii(:,1)) + MATMUL(ckmat(:,:,2),cpsii(:,2)) + MATMUL(ckmat(:,:,3),cpsii(:,3)) &
           + MATMUL(ckmat(:,:,4),cpsii(:,4)) + MATMUL(ckmat(:,:,5),cpsii(:,5)) + MATMUL(ckmat(:,:,6),cpsii(:,6)) &
           + MATMUL(ckmat(:,:,7),cpsii(:,7)) + MATMUL(ckmat(:,:,8),cpsii(:,8))
END IF

! Use LAPACK solvers if matrix was factored only
IF (tpose == 1) THEN
   IF (sym == 1) THEN
      CALL dpotrs('U',cneq,1,cjmat,cneq,f,cneq,info)
      IF (info /= 0) THEN
         WRITE (8,'(//,1X,A)') "ERROR: dpotrs cannot solve."
         STOP
      END IF
   ELSE
      CALL dgetrs('T',cneq,1,cjmat,cneq,cpiv,f,cneq,info)
      IF (info /= 0) THEN
         WRITE (8,'(//,1X,A)') "ERROR: dgetrs cannot solve."
         STOP
      END IF
   END IF
ELSE
   IF (sym == 1) THEN
      CALL dpotrs('U',cneq,1,cjmat,cneq,f,cneq,info)
      IF (info /= 0) THEN
         WRITE (8,'(//,1X,A)') "ERROR: dpotrs cannot solve."
         STOP
      END IF
   ELSE
      CALL dgetrs('N',cneq,1,cjmat,cneq,cpiv,f,cneq,info)
      IF (info /= 0) THEN
         WRITE (8,'(//,1X,A)') "ERROR: dgetrs cannot solve."
         STOP
      END IF
   END IF
END IF

! Compute the new angular flux values at the boundaries from jpsi and kpsi
q = f

! Set up the loop independent variables
ix1 = 1
ix2 = cxys
ix3 = ix2 + 1
ix4 = ix2 + cxzs
ix5 = ix4 + 1
ix6 = ix4 + cyzs
! Find psio depending on 'tpose' flag
IF (tpose == 1) THEN
   DO j = 1, 8
      DO i = 1, cbcs
         cpsio(i,j) = DOT_PRODUCT(cjpsi(:,i,j),q)
      END DO
      DO i = 1, apo
         kz = (i-1)*cxys + 1
         lz = i*cxys
         ky = apo*cxys + (i-1)*cxzs + 1
         ly = apo*cxys + i*cxzs
         kx = apo*(cxys+cxzs) + (i-1)*cyzs + 1
         lx = apo*(cxys+cxzs) + i*cyzs
         jndx = 0
         ! zbc's in psio
         DO k = kz, lz
            jndx = jndx + 1
            cpsio(k,j) = cpsio(k,j) + DOT_PRODUCT(ckpsi(ix1:ix2,jndx,i,j),cpsii(kz:lz,j)) &
                                    + DOT_PRODUCT(ckpsi(ix3:ix4,jndx,i,j),cpsii(ky:ly,j)) &
                                    + DOT_PRODUCT(ckpsi(ix5:ix6,jndx,i,j),cpsii(kx:lx,j))
         END DO
         ! ybc's in psio
         DO k = ky, ly
            jndx = jndx + 1
            cpsio(k,j) = cpsio(k,j) + DOT_PRODUCT(ckpsi(ix1:ix2,jndx,i,j),cpsii(kz:lz,j)) &
                                    + DOT_PRODUCT(ckpsi(ix3:ix4,jndx,i,j),cpsii(ky:ly,j)) &
                                    + DOT_PRODUCT(ckpsi(ix5:ix6,jndx,i,j),cpsii(kx:lx,j))
         END DO
         ! xbc's in psio
         DO k = kx, lx
            jndx = jndx + 1
            cpsio(k,j) = cpsio(k,j) + DOT_PRODUCT(ckpsi(ix1:ix2,jndx,i,j),cpsii(kz:lz,j)) &
                                    + DOT_PRODUCT(ckpsi(ix3:ix4,jndx,i,j),cpsii(ky:ly,j)) &
                                    + DOT_PRODUCT(ckpsi(ix5:ix6,jndx,i,j),cpsii(kx:lx,j))
         END DO
      END DO
   END DO
ELSE
   DO j = 1, 8
      cpsio(:,j) = MATMUL(cjpsi(:,:,j),q)
      DO i = 1, apo
         kz = (i-1)*cxys + 1
         lz = i*cxys
         ky = apo*cxys + (i-1)*cxzs + 1
         ly = apo*cxys + i*cxzs
         kx = apo*(cxys+cxzs) + (i-1)*cyzs + 1
         lx = apo*(cxys+cxzs) + i*cyzs
         ! zbc's in psio
         cpsio(kz:lz,j) = cpsio(kz:lz,j) + MATMUL(ckpsi(ix1:ix2,ix1:ix2,i,j),cpsii(kz:lz,j)) &
                                         + MATMUL(ckpsi(ix1:ix2,ix3:ix4,i,j),cpsii(ky:ly,j)) &
                                         + MATMUL(ckpsi(ix1:ix2,ix5:ix6,i,j),cpsii(kx:lx,j))
         ! ybc's in psio
         cpsio(ky:ly,j) = cpsio(ky:ly,j) + MATMUL(ckpsi(ix3:ix4,ix1:ix2,i,j),cpsii(kz:lz,j)) &
                                         + MATMUL(ckpsi(ix3:ix4,ix3:ix4,i,j),cpsii(ky:ly,j)) &
                                         + MATMUL(ckpsi(ix3:ix4,ix5:ix6,i,j),cpsii(kx:lx,j))
         ! xbc's in psio
         cpsio(kx:lx,j) = cpsio(kx:lx,j) + MATMUL(ckpsi(ix5:ix6,ix1:ix2,i,j),cpsii(kz:lz,j)) &
                                         + MATMUL(ckpsi(ix5:ix6,ix3:ix4,i,j),cpsii(ky:ly,j)) &
                                         + MATMUL(ckpsi(ix5:ix6,ix5:ix6,i,j),cpsii(kx:lx,j))
      END DO
   END DO
END IF
! Lastly add in the RHS
cpsio = cpsio + rrpsi

RETURN
END SUBROUTINE cidot
