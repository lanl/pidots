SUBROUTINE interp(v,spx,spy,spz,kc,jc,ic,cneq)

!------------------------------------------------------------------
!
! Coarse grid processes send correction to fine grid processes.
! Then interpolate the correction using phir.
!
!------------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER, INTENT(IN) :: v, spx, spy, spz, kc, jc, ic, cneq
INTEGER :: ka, ja, ia, ijc
INTEGER :: kt, jt, it, k, j, i, ix1, jx1, ix2, jx2, ieq, jeq
REAL*8, DIMENSION(cneq) :: rsv

! Use kc/jc/ic to determine other needed coarsening parameters
ka = kc - 1
ja = jc - 1
ia = ic - 1
ijc = ic*jc

! First get the correction sent to the fine grid processes
CALL corrscat(v,spx,spy,spz,kc,jc,ic,cneq,rsv)

! Now interpolate
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
         phi(ieq,v) = phi(ieq,v) + phir(ieq,v)*rsv(jeq)
      END DO
   END DO
END DO

RETURN
END SUBROUTINE interp
