SUBROUTINE interp(kc,jc,ic,cneq,phi,f,phir)

!------------------------------------------------------------------
!
! Coarse grid processes send correction to fine grid processes.
! Then interpolate the correction using phir.
!
!------------------------------------------------------------------

USE invar
IMPLICIT NONE
INTEGER, INTENT(IN) :: kc, jc, ic, cneq
INTEGER :: ka, ja, ia, ijc, sdi
INTEGER :: kt, jt, it, k, j, i, ix1, jx1, ix2, jx2, ieq, jeq
REAL*8, DIMENSION(cneq,nsdp) :: rsv
REAL*8, DIMENSION(neq,nsdp), INTENT(INOUT) :: phi
REAL*8, DIMENSION(neq,nsdp), INTENT(IN) :: phir
REAL*8, DIMENSION(neq), INTENT(IN) :: f

! Use kc/jc/ic to determine other needed coarsening parameters
ka = kc - 1
ja = jc - 1
ia = ic - 1
ijc = ic*jc

! First get the correction sent to the fine grid processes
CALL corrscat(kc,jc,ic,cneq,rsv,f)

! Now interpolate
DO sdi = 1, nsdp
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
            phi(ieq,sdi) = phi(ieq,sdi) + phir(ieq,sdi)*rsv(jeq,sdi)
         END DO
      END DO
   END DO
END DO

RETURN
END SUBROUTINE interp
