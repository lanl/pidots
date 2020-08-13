SUBROUTINE restrict(rphi,rsv,rpsi,rpave)

!-------------------------------------------------------------
!
!  Restrict the fine grid data for use on coarse grid
!
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER :: o, n, k, j, i, jeq, ieq, ierr
REAL*8 :: tmp
REAL*8, DIMENSION(neq), INTENT(IN) :: rphi
REAL*8, DIMENSION(cneq), INTENT(OUT) :: rsv
REAL*8, DIMENSION(bcs,8), INTENT(IN) :: rpsi
REAL*8, DIMENSION(3*apo,8), INTENT(OUT) :: rpave

INCLUDE 'mpif.h'

! phi
! Get average of rphi
rsv = 0.0
tmp = SUM(rphi)/neq
CALL MPI_GATHER(tmp,1,MPI_DOUBLE_PRECISION,rsv,1,MPI_DOUBLE_PRECISION,root,allcomm,ierr)

!! Get the average phi as well and store ratios for interpolation later
tmp = SUM(phi)/neq
phir = phi/tmp

!print*, irank, rphi
!phir = rphi/tmp

!print *, irank, phir

!if (nrank == root) print*, rphi, rsv

! psi
rpave = 0.0
DO o = 1, 8
   DO n = 1, apo
      ieq = n
      DO j = 1, ny
         DO i = 1, nx
            jeq = (n-1)*xys + (j-1)*nx + i
            rpave(ieq,o) = rpave(ieq,o) + rpsi(jeq,o)
         END DO
      END DO
      rpave(ieq,o) = rpave(ieq,o)/xys
   END DO
   DO n = 1, apo
      ieq = apo + n
      DO k = 1, nz
         DO i = 1, nx
            jeq = apo*xys + (n-1)*xzs + (k-1)*nx + i
            rpave(ieq,o) = rpave(ieq,o) + rpsi(jeq,o)
         END DO
      END DO
      rpave(ieq,o) = rpave(ieq,o)/xzs
   END DO
   DO n = 1, apo
      ieq = 2*apo + n
      DO k = 1, nz
         DO j = 1, ny
            jeq = apo*(xys+xzs) + (n-1)*yzs + (k-1)*ny + j
            rpave(ieq,o) = rpave(ieq,o) + rpsi(jeq,o)
         END DO
      END DO
      rpave(ieq,o) = rpave(ieq,o)/yzs
   END DO
END DO
CALL rpsigthr(rpave)


RETURN
END SUBROUTINE restrict
