SUBROUTINE interpf(v,spx,spy,spz,kc,jc,ic,cneq,f1,f2,fr)

!------------------------------------------------------------------
!
! Coarse grid processes send correction to fine grid processes.
! Then interpolate the correction using phir.
!
!------------------------------------------------------------------

USE invar
IMPLICIT NONE
INTEGER, INTENT(IN) :: v, spx, spy, spz, kc, jc, ic, cneq
INTEGER :: ka, ja, ia, ijc
INTEGER :: kt, jt, it, k, j, i, ix1, jx1, ix2, jx2, ieq, jeq
REAL*8, DIMENSION(cneq) :: rsv
REAL*8, DIMENSION(neq), INTENT(INOUT) :: f1
REAL*8, DIMENSION(neq), INTENT(IN) :: f2, fr

! Use kc/jc/ic to determine other needed coarsening parameters
ka = kc - 1
ja = jc - 1
ia = ic - 1
ijc = ic*jc

! First get the correction sent to the fine grid processes
CALL corrscatf(v,spx,spy,spz,kc,jc,ic,cneq,rsv,f2)

! Now interpolate
DO k = 1, sdnz
   kt = (k+ka)/kc
   ix1 = (k-1)*xys
   jx1 = (kt-1)*xys/ijc
   DO j = 1, sdny
      jt = (j+ja)/jc
      ix2 = ix1 + (j-1)*sdnx
      jx2 = jx1 + (jt-1)*sdnx/ic
      DO i = 1, sdnx
         it = (i+ia)/ic
         ieq = ix2 + i
         jeq = jx2 + it
         f1(ieq) = f1(ieq) + fr(ieq)*rsv(jeq)
      END DO
   END DO
END DO

RETURN
END SUBROUTINE interpf
