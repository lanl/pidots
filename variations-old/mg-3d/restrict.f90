SUBROUTINE restrict(v,spx,spy,spz,kc,jc,ic,cneq,cbcs,crsf,sdcf,rphi,rpsi)

!-------------------------------------------------------------
!
!  Restrict the fine grid data for use on coarse grid
!
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER, INTENT(IN) :: v, spx, spy, spz, kc, jc, ic, cneq, cbcs, crsf
INTEGER :: o, n, k, j, i, jeq, ieq, it, jt, kt
INTEGER :: ka, ja, ia, ijc, ikc, jkc, ix0, jx0, ix1, jx1, ix2, jx2
REAL*8, INTENT(IN) :: sdcf
REAL*8 :: tmp
REAL*8, DIMENSION(neq), INTENT(IN) :: rphi
REAL*8, DIMENSION(cneq) :: rsv, avg
REAL*8, DIMENSION(bcs,8), INTENT(IN) :: rpsi
REAL*8, DIMENSION(cbcs,8) :: rpave

! Use kc/jc,ic to determine other coarsening parameters
ka = kc - 1
ja = jc - 1
ia = ic - 1
ijc = ic*jc
ikc = ic*kc
jkc = jc*kc

! phi
! Get average of rphi of fine grid cells
tmp = REAL(kc*jc*ic)
rsv = 0.0
DO k = 1, nz
   kt = (k+ka)/kc
   ix1 = (k-1)*xys
   jx1 = (kt-1)*xys/ijc
   DO j = 1, ny
      jt = (j+ja)/jc
      ix2 = ix1 + (j-1)*nx
      jx2 = jx1 + (jt-1)*nx/ic
      DO i = 1, nx
         it = (i+ia)/ic
         ieq = ix2 + i
         jeq = jx2 + it
         rsv(jeq) = rsv(jeq) + rphi(ieq)/tmp
      END DO
   END DO
END DO
! Call for a communication that collects rsv on processors that will go to next v-cycle stage
CALL rsvgthr(v,spx,spy,spz,kc,jc,ic,cneq,crsf,sdcf,rsv)

! Do same for to sum and get an average for phi
avg = 0.0
DO k = 1, nz
   kt = (k+ka)/kc
   ix1 = (k-1)*xys
   jx1 = (kt-1)*xys/ijc
   DO j = 1, ny
      jt = (j+ja)/jc
      ix2 = ix1 + (j-1)*nx
      jx2 = jx1 + (jt-1)*nx/ic
      DO i = 1, nx
         it = (i+ia)/ic
         ieq = ix2 + i
         jeq = jx2 + it
         avg(jeq) = avg(jeq) + phi(ieq,v)/tmp
      END DO
   END DO
END DO
! Use the avg vector and phi to get a ratio that's used for interpolation
DO k = 1, nz
   kt = (k+ka)/kc
   ix1 = (k-1)*xys
   jx1 = (kt-1)*xys/ijc
   DO j = 1, ny
      jt = (j+ja)/jc
      ix2 = ix1 + (j-1)*nx
      jx2 = jx1 + (jt-1)*nx/ic
      DO i = 1, nx
         it = (i+ia)/ic
         ieq = ix2 + i
         jeq = jx2 + it
         phir(ieq,v) = phi(ieq,v)/avg(jeq)
      END DO
   END DO
END DO

! psi
rpave = 0.0
DO o = 1, 8
   tmp = REAL(ijc)
   DO n = 1, apo
      ix1 = (n-1)*xys
      jx1 = (n-1)*xys/ijc
      DO j = 1, ny
         jt = (j+ja)/jc
         ix2 = ix1 + (j-1)*nx
         jx2 = jx1 + (jt-1)*nx/ic
         DO i = 1, nx
            it = (i+ia)/ic
            ieq = ix2 + i
            jeq = jx2 + it
            rpave(jeq,o) = rpave(jeq,o) + rpsi(ieq,o)/tmp
         END DO
      END DO
   END DO
   tmp = REAL(ikc)
   ix0 = apo*xys
   jx0 = apo*xys/ijc
   DO n = 1, apo
      ix1 = ix0 + (n-1)*xzs
      jx1 = jx0 + (n-1)*xzs/ikc
      DO k = 1, nz
         kt = (k+ka)/kc
         ix2 = ix1 + (k-1)*nx
         jx2 = jx1 + (kt-1)*nx/ic
         DO i = 1, nx
            it = (i+ia)/ic
            ieq = ix2 + i
            jeq = jx2 + it
            rpave(jeq,o) = rpave(jeq,o) + rpsi(ieq,o)/tmp
         END DO
      END DO
   END DO
   tmp = REAL(jkc)
   ix0 = apo*(xys+xzs)
   jx0 = apo*(xys/ijc + xzs/ikc)
   DO n = 1, apo
      ix1 = ix0 + (n-1)*yzs
      jx1 = jx0 + (n-1)*yzs/jkc
      DO k = 1, nz
         kt = (k+ka)/kc
         ix2 = ix1 + (k-1)*ny
         jx2 = jx1 + (kt-1)*ny/jc
         DO j = 1, ny
            jt = (j+ja)/jc
            ieq = ix2 + j
            jeq = jx2 + jt
            rpave(jeq,o) = rpave(jeq,o) + rpsi(ieq,o)/tmp
         END DO
      END DO
   END DO
END DO
CALL rpsigthr(v,spx,spy,spz,kc,jc,ic,cbcs,crsf,sdcf,rpave)

RETURN
END SUBROUTINE restrict
