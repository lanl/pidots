SUBROUTINE restrict(kc,jc,ic,cneq,cbcs,sdcf,rphi,rpsi,phi,phir,fsv,fav)

!-------------------------------------------------------------
!
!  Restrict the fine grid data for use on coarse grid
!
!-------------------------------------------------------------

USE invar
IMPLICIT NONE
INTEGER, INTENT(IN) :: kc, jc, ic, cneq, cbcs
INTEGER :: o, n, k, j, i, jeq, ieq, it, jt, kt, sdi
INTEGER :: ka, ja, ia, ijc, ikc, jkc, ix0, jx0, ix1, jx1, ix2, jx2
REAL*8, INTENT(IN) :: sdcf
REAL*8 :: tmp
REAL*8, DIMENSION(neq,nsdp), INTENT(IN) :: phi, rphi
REAL*8, DIMENSION(neq,nsdp), INTENT(OUT) :: phir
REAL*8, DIMENSION(neq), INTENT(OUT) :: fsv
REAL*8, DIMENSION(cneq,nsdp) :: rsv, avg
REAL*8, DIMENSION(bcs,8,nsdp), INTENT(IN) :: rpsi
REAL*8, DIMENSION(bcs,8), INTENT(OUT) :: fav
REAL*8, DIMENSION(cbcs,8,nsdp) :: rpave

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
            rsv(jeq,sdi) = rsv(jeq,sdi) + rphi(ieq,sdi)/tmp
         END DO
      END DO
   END DO
END DO

! Call for a communication that collects rsv on processors that will go to next v-cycle stage
CALL rsvgthr(kc,jc,ic,cneq,sdcf,rsv,fsv)

! Do same for to sum and get an average for phi
avg = 0.0
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
            avg(jeq,sdi) = avg(jeq,sdi) + phi(ieq,sdi)/tmp
         END DO
      END DO
   END DO
END DO
! Use the avg vector and phi to get a ratio that's used for interpolation
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
            phir(ieq,sdi) = phi(ieq,sdi)/avg(jeq,sdi)
         END DO
      END DO
   END DO
END DO

! psi
rpave = 0.0
DO sdi = 1, nsdp
   DO o = 1, 8
      tmp = REAL(ijc)
      DO n = 1, apo
         ix1 = (n-1)*xys
         jx1 = (n-1)*xys/ijc
         DO j = 1, sdny
            jt = (j+ja)/jc
            ix2 = ix1 + (j-1)*sdnx
            jx2 = jx1 + (jt-1)*sdnx/ic
            DO i = 1, sdnx
               it = (i+ia)/ic
               ieq = ix2 + i
               jeq = jx2 + it
               rpave(jeq,o,sdi) = rpave(jeq,o,sdi) + rpsi(ieq,o,sdi)/tmp
            END DO
         END DO
      END DO
      tmp = REAL(ikc)
      ix0 = apo*xys
      jx0 = apo*xys/ijc
      DO n = 1, apo
         ix1 = ix0 + (n-1)*xzs
         jx1 = jx0 + (n-1)*xzs/ikc
         DO k = 1, sdnz
            kt = (k+ka)/kc
            ix2 = ix1 + (k-1)*sdnx
            jx2 = jx1 + (kt-1)*sdnx/ic
            DO i = 1, sdnx
               it = (i+ia)/ic
               ieq = ix2 + i
               jeq = jx2 + it
               rpave(jeq,o,sdi) = rpave(jeq,o,sdi) + rpsi(ieq,o,sdi)/tmp
            END DO
         END DO
      END DO
      tmp = REAL(jkc)
      ix0 = apo*(xys+xzs)
      jx0 = apo*(xys/ijc + xzs/ikc)
      DO n = 1, apo
         ix1 = ix0 + (n-1)*yzs
         jx1 = jx0 + (n-1)*yzs/jkc
         DO k = 1, sdnz
            kt = (k+ka)/kc
            ix2 = ix1 + (k-1)*sdny
            jx2 = jx1 + (kt-1)*sdny/jc
            DO j = 1, sdny
               jt = (j+ja)/jc
               ieq = ix2 + j
               jeq = jx2 + jt
               rpave(jeq,o,sdi) = rpave(jeq,o,sdi) + rpsi(ieq,o,sdi)/tmp
            END DO
         END DO
      END DO
   END DO
END DO
CALL rpsigthr(kc,jc,ic,cbcs,rpave,fav)

RETURN
END SUBROUTINE restrict
