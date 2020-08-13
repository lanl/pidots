SUBROUTINE rpsigthr(kc,jc,ic,cbcs,rpave,fav)

!-------------------------------------------------------------
!
! Gathers all the rpave vectors and puts them into rrpsi 
! 
!-------------------------------------------------------------

USE invar
IMPLICIT NONE
INTEGER, INTENT(IN) :: kc, jc, ic, cbcs
INTEGER :: i, j, k, it, jt, kt, n, o, ieq, jeq, sdi, kk, jj, ii
INTEGER :: ix0, jx0, ix1, jx1, ix2, jx2, ijc, ikc, jkc
INTEGER :: tz, ty, tx, txy, txz, tyz, tbcs, kcc, jcc, icc
REAL*8 :: xyf, xzf, yzf
REAL*8, DIMENSION(cbcs,8,nsdp), INTENT(IN) :: rpave
REAL*8, DIMENSION(:,:), ALLOCATABLE :: tv
REAL*8, DIMENSION(bcs,8), INTENT(OUT) :: fav

! Need zcf/ycf/xcf for communication
kcc = zcf(0)
jcc = ycf(0)
icc = xcf(0)

! Set up the temporary matrix
tz = sdnz
ty = sdny
tx = sdnx
IF (kcc > sdnz) tz = kcc
IF (jcc > sdny) ty = jcc
IF (icc > sdnx) tx = icc
txy = tx*ty
txz = tx*tz
tyz = ty*tz
tbcs = apo*(txy + txz + tyz)
ALLOCATE(tv(tbcs,8))
xyf = REAL(sdnx*sdny)/REAL(icc*jcc)
xzf = REAL(sdnx*sdnz)/REAL(icc*kcc)
yzf = REAL(sdny*sdnz)/REAL(jcc*kcc)
ijc = ic*jc
ikc = ic*kc
jkc = jc*kc

! Copy all the rpave to one large matrix tv
! Only copy if sub-domain boundary is part of the coarser sub-domain boundary
sdi = 0
DO kk = 1, kcc
   kt = (kk-1)*sdnz/kc
   DO jj = 1, jcc
      jt = (jj-1)*sdny/jc
      DO ii = 1, icc
         it = (ii-1)*sdnx/ic
         sdi = sdi + 1
         IF (kk == kcc) THEN
            DO n = 1, apo
               ix1 = (n-1)*xys/ijc
               jx1 = (n-1)*txy
               DO j = jt + 1, jt + sdny/jc
                  ix2 = ix1 + (j-jt-1)*sdnx/ic
                  jx2 = jx1 + (j-1)*tx
                  DO i = it + 1, it + sdnx/ic
                     ieq = ix2 + i-it
                     jeq = jx2 + i
                     tv(jeq,1) = rpave(ieq,1,sdi)
                     tv(jeq,2) = rpave(ieq,2,sdi)
                     tv(jeq,3) = rpave(ieq,3,sdi)
                     tv(jeq,4) = rpave(ieq,4,sdi)
                  END DO
               END DO
            END DO
         ELSE IF (kk == 1) THEN
            DO n = 1, apo
               ix1 = (n-1)*xys/ijc
               jx1 = (n-1)*txy
               DO j = jt + 1, jt + sdny/jc
                  ix2 = ix1 + (j-jt-1)*sdnx/ic
                  jx2 = jx1 + (j-1)*tx
                  DO i = it + 1, it + sdnx/ic
                     ieq = ix2 + i-it
                     jeq = jx2 + i
                     tv(jeq,5) = rpave(ieq,5,sdi)
                     tv(jeq,6) = rpave(ieq,6,sdi)
                     tv(jeq,7) = rpave(ieq,7,sdi)
                     tv(jeq,8) = rpave(ieq,8,sdi)
                  END DO
               END DO
            END DO
         END IF
         IF (jj == jcc) THEN
            ix0 = apo*xys/ijc
            jx0 = apo*txy
            DO n = 1, apo
               ix1 = ix0 + (n-1)*xzs/ikc
               jx1 = jx0 + (n-1)*txz
               DO k = kt + 1, kt + sdnz/kc
                  ix2 = ix1 + (k-kt-1)*sdnx/ic
                  jx2 = jx1 + (k-1)*tx
                  DO i = it + 1, it + sdnx/ic
                     ieq = ix2 + i-it
                     jeq = jx2 + i
                     tv(jeq,1) = rpave(ieq,1,sdi)
                     tv(jeq,2) = rpave(ieq,2,sdi)
                     tv(jeq,5) = rpave(ieq,5,sdi)
                     tv(jeq,6) = rpave(ieq,6,sdi)
                  END DO
               END DO
            END DO
         ELSE IF (jj == 1) THEN
            ix0 = apo*xys/ijc
            jx0 = apo*txy
            DO n = 1, apo
               ix1 = ix0 + (n-1)*xzs/ikc
               jx1 = jx0 + (n-1)*txz
               DO k = kt + 1, kt + sdnz/kc
                  ix2 = ix1 + (k-kt-1)*sdnx/ic
                  jx2 = jx1 + (k-1)*tx
                  DO i = it + 1, it + sdnx/ic
                     ieq = ix2 + i-it
                     jeq = jx2 + i
                     tv(jeq,3) = rpave(ieq,3,sdi)
                     tv(jeq,4) = rpave(ieq,4,sdi)
                     tv(jeq,7) = rpave(ieq,7,sdi)
                     tv(jeq,8) = rpave(ieq,8,sdi)
                  END DO
               END DO
            END DO
         END IF
         IF (ii == icc) THEN
            ix0 = apo*(xys/ijc + xzs/ikc)
            jx0 = apo*(txy + txz)
            DO n = 1, apo
               ix1 = ix0 + (n-1)*yzs/jkc
               jx1 = jx0 + (n-1)*tyz
               DO k = kt + 1, kt + sdnz/kc
                  ix2 = ix1 + (k-kt-1)*sdny/jc
                  jx2 = jx1 + (k-1)*ty
                  DO j = jt + 1, jt + sdny/jc
                     ieq = ix2 + j-jt
                     jeq = jx2 + j
                     tv(jeq,1) = rpave(ieq,1,sdi)
                     tv(jeq,4) = rpave(ieq,4,sdi)
                     tv(jeq,5) = rpave(ieq,5,sdi)
                     tv(jeq,8) = rpave(ieq,8,sdi)
                  END DO
               END DO
            END DO
         ELSE IF (ii == 1) THEN
            ix0 = apo*(xys/ijc + xzs/ikc)
            jx0 = apo*(txy + txz)
            DO n = 1, apo
               ix1 = ix0 + (n-1)*yzs/jkc
               jx1 = jx0 + (n-1)*tyz
               DO k = kt + 1, kt + sdnz/kc
                  ix2 = ix1 + (k-kt-1)*sdny/jc
                  jx2 = jx1 + (k-1)*ty
                  DO j = jt + 1, jt + sdny/jc
                     ieq = ix2 + j-jt
                     jeq = jx2 + j
                     tv(jeq,2) = rpave(ieq,2,sdi)
                     tv(jeq,3) = rpave(ieq,3,sdi)
                     tv(jeq,6) = rpave(ieq,6,sdi)
                     tv(jeq,7) = rpave(ieq,7,sdi)
                  END DO
               END DO
            END DO
         END IF
      END DO
   END DO
END DO

! Reinitialize av
fav = 0.0
! Now compute average for combined cells/sub-domains
! Place into fav array for use in next v-cycle stage
DO n = 1, apo
   ix1 = (n-1)*txy
   jx1 = (n-1)*xys
   DO j = 1, ty
      ix2 = ix1 + (j-1)*tx
      jx2 = jx1 + (j-1)*sdnx
      IF (jcc > sdny) jx2 = jx1 + ((j-1+jcc/sdny)*sdny/jcc-1)*sdnx
      DO i = 1, tx
         ieq = ix2 + i
         jeq = jx2 + i
         IF (icc > sdnx) jeq = jx2 + (i-1+icc/sdnx)*sdnx/icc
         DO o = 1, 8
            fav(jeq,o) = fav(jeq,o) + tv(ieq,o)*xyf
         END DO
      END DO
   END DO
END DO
ix0 = apo*txy
jx0 = apo*xys
DO n = 1, apo
   ix1 = ix0 + (n-1)*txz
   jx1 = jx0 + (n-1)*xzs
   DO k = 1, tz
      ix2 = ix1 + (k-1)*tx
      jx2 = jx1 + (k-1)*sdnx
      IF (kcc > sdnz) jx2 = jx1 + ((k-1+kcc/sdnz)*sdnz/kcc-1)*sdnx
      DO i = 1, tx
         ieq = ix2 + i
         jeq = jx2 + i
         IF (icc > sdnx) jeq = jx2 + (i-1+icc/sdnx)*sdnx/icc
         DO o = 1, 8
            fav(jeq,o) = fav(jeq,o) + tv(ieq,o)*xzf
         END DO
      END DO
   END DO
END DO
ix0 = apo*(txy + txz)
jx0 = apo*(xys + xzs)
DO n = 1, apo
   ix1 = ix0 + (n-1)*tyz
   jx1 = jx0 + (n-1)*yzs
   DO k = 1, tz
      ix2 = ix1 + (k-1)*ty
      jx2 = jx1 + (k-1)*sdny
      IF (kcc > sdnz) jx2 = jx1 + ((k-1+kcc/sdnz)*sdnz/kcc-1)*sdny
      DO j = 1, ty
         ieq = ix2 + j
         jeq = jx2 + j
         IF (jcc > sdny) jeq = jx2 + (j-1+jcc/sdny)*sdny/jcc
         DO o = 1, 8
            fav(jeq,o) = fav(jeq,o) + tv(ieq,o)*yzf
         END DO
      END DO
   END DO
END DO

DEALLOCATE(tv)

RETURN
END SUBROUTINE rpsigthr
