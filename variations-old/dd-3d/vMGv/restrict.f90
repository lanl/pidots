SUBROUTINE restrict(v,sp,rphi,rpsi, vit)

!-------------------------------------------------------------
!
!  Restrict the fine grid data for use on coarse grid
!
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER, INTENT(IN) :: v, sp, vit
INTEGER :: o, n, k, j, i, jeq, ieq, it, jt, kt
REAL*8 :: tmp
REAL*8, DIMENSION(neq), INTENT(IN) :: rphi
REAL*8, DIMENSION(cneq) :: rsv, avg
REAL*8, DIMENSION(bcs,8), INTENT(IN) :: rpsi
REAL*8, DIMENSION(cbcs,8) :: rpave

! phi
! Get average of rphi of 8 adjacent cells
tmp = 8.0
rsv = 0.0
DO k = 1, nz
   kt = (k+1)/2
   DO j = 1, ny
      jt = (j+1)/2
      DO i = 1, nx
         it = (i+1)/2
         ieq = (k-1)*xys + (j-1)*nx + i
         jeq = (kt-1)*xys/4 + (jt-1)*nx/2 + it
         rsv(jeq) = rsv(jeq) + rphi(ieq)/tmp
      END DO
   END DO
END DO
! Call for a communication that collects rsv on processors that will go to next v-cycle stage
CALL rsvgthr(v,sp,rsv)

!if (vit ==3 .and. v==2) print*, nrank, sv(:,3)

!print*, xrank, yrank, zrank, nrank, "out of rsvgthr"

! Do same for to sum and get an average for phi
avg = 0.0
DO k = 1, nz
   kt = (k+1)/2
   DO j = 1, ny
      jt = (j+1)/2
      DO i = 1, nx
         it = (i+1)/2
         ieq = (k-1)*xys + (j-1)*nx + i
         jeq = (kt-1)*xys/4 + (jt-1)*nx/2 + it
         avg(jeq) = avg(jeq) + phi(ieq,v)/tmp
      END DO
   END DO
END DO
! Use the avg vector and phi to get a ratio that's used for interpolation
DO k = 1, nz
   kt = (k+1)/2
   DO j = 1, ny
      jt = (j+1)/2
      DO i = 1, nx
         it = (i+1)/2
         ieq = (k-1)*xys + (j-1)*nx + i
         jeq = (kt-1)*xys/4 + (jt-1)*nx/2 + it
         phir(ieq,v) = phi(ieq,v)/avg(jeq)
      END DO
   END DO
END DO

! psi
tmp = 4.0
rpave = 0.0
DO o = 1, 8
   DO n = 1, apo
      DO j = 1, ny
         jt = (j+1)/2
         DO i = 1, nx
            it = (i+1)/2
            ieq = (n-1)*xys + (j-1)*nx + i
            jeq = (n-1)*xys/4 + (jt-1)*nx/2 + it
            rpave(jeq,o) = rpave(jeq,o) + rpsi(ieq,o)/tmp
         END DO
      END DO
   END DO
   DO n = 1, apo
      DO k = 1, nz
         kt = (k+1)/2
         DO i = 1, nx
            it = (i+1)/2
            ieq = apo*xys + (n-1)*xzs + (k-1)*nx + i
            jeq = apo*xys/4 + (n-1)*xzs/4 + (kt-1)*nx/2 + it
            rpave(jeq,o) = rpave(jeq,o) + rpsi(ieq,o)/tmp
         END DO
      END DO
   END DO
   DO n = 1, apo
      DO k = 1, nz
         kt = (k+1)/2
         DO j = 1, ny
            jt = (j+1)/2
            ieq = apo*(xys+xzs) + (n-1)*yzs + (k-1)*ny + j
            jeq = apo*(xys+xzs)/4 + (n-1)*yzs/4 + (kt-1)*ny/2 + jt
            rpave(jeq,o) = rpave(jeq,o) + rpsi(ieq,o)/tmp
         END DO
      END DO
   END DO
END DO
CALL rpsigthr(v,sp,rpave)

RETURN
END SUBROUTINE restrict
