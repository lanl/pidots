SUBROUTINE corrscat(kc,jc,ic,cneq,rsv,f)

!------------------------------------------------------------------
!
! Coarse grid processes send correction to fine grid processes.
!
!------------------------------------------------------------------

USE invar
IMPLICIT NONE
INTEGER, INTENT(IN) :: kc, jc, ic, cneq
INTEGER :: kk, jj, ii, kt, jt, it, k, j, i, sdi
INTEGER :: kcc, jcc, icc, ijc, ix1, ix2, jx1, jx2, ieq, jeq
REAL*8, DIMENSION(cneq,nsdp), INTENT(OUT) :: rsv
REAL*8, DIMENSION(neq), INTENT(IN) :: f

! Need zcf/ycf/xcf for communication
kcc = zcf(0)
jcc = ycf(0)
icc = xcf(0)
ijc = ic*jc

sdi = 0
DO kk = 1, kcc
   kt = (kk-1)*sdnz/kc
   DO jj = 1, jcc
      jt = (jj-1)*sdny/jc
      DO ii = 1, icc
         it = (ii-1)*sdnx/ic
         sdi = sdi + 1
         DO k = kt + 1, kt + sdnz/kc
            ix1 = (k-kt-1)*xys/ijc
            jx1 = (k-1)*xys
            IF (kcc > sdnz) jx1 = ((k-1+kcc/sdnz)*sdnz/kcc-1)*xys
            DO j = jt + 1, jt + sdny/jc
               ix2 = ix1 + (j-jt-1)*sdnx/ic
               jx2 = jx1 + (j-1)*sdnx
               IF (jcc > sdny) jx2 = jx1 + ((j-1+jcc/sdny)*sdny/jcc-1)*sdnx
               DO i = it + 1, it + sdnx/ic
                  ieq = ix2 + i-it
                  jeq = jx2 + i
                  IF (icc > sdnx) jeq = jx2 + (i-1+icc/sdnx)*sdnx/icc
                  rsv(ieq,sdi) = f(jeq)
               END DO
            END DO
         END DO
      END DO
   END DO
END DO

RETURN
END SUBROUTINE corrscat
