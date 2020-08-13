SUBROUTINE rsvgthr(kc,jc,ic,cneq,sdcf,rsv,fsv)

!-------------------------------------------------------------
!
! Gathers rsv vectors and puts solution into sv(:,v+1)
! 
!-------------------------------------------------------------

USE invar
IMPLICIT NONE
INTEGER, INTENT(IN) :: kc, jc, ic, cneq
INTEGER :: kt, jt, it, k, j, i, ix1, jx1, ix2, jx2, ieq, sdi
INTEGER :: kcc, jcc, icc, ijc, tz, ty, tx, ii, jj, kk
REAL*8, INTENT(IN) :: sdcf
REAL*8, DIMENSION(cneq,nsdp), INTENT(IN) :: rsv
REAL*8, DIMENSION(neq), INTENT(OUT) :: fsv
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: tv

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
ALLOCATE(tv(tx,ty,tz))
ijc = ic*jc

! Copy all the rsv's to one large matrix tv
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
            DO j = jt + 1, jt + sdny/jc
               ix2 = ix1 + (j-jt-1)*sdnx/ic
               DO i = it + 1, it + sdnx/ic
                  ieq = ix2 + i-it
                  tv(i,j,k) = rsv(ieq,sdi)
               END DO
            END DO
         END DO
      END DO
   END DO
END DO

! Reinitialize fsv using tv
fsv = 0.0
! Now compute average for combined cells/sub-domains
DO k = 1, tz
   ix1 = (k-1)*xys
   IF (kcc > sdnz) ix1 = ((k-1+kcc/sdnz)*sdnz/kcc-1)*xys
   DO j = 1, ty
      ix2 = ix1 + (j-1)*sdnx
      IF (jcc > sdny) ix2 = ix1 + ((j-1+jcc/sdny)*sdny/jcc-1)*sdnx
      DO i = 1, tx
         ieq = ix2 + i
         IF (icc > sdnx) ieq = ix2 + (i-1+icc/sdnx)*sdnx/icc
         fsv(ieq) = fsv(ieq) + tv(i,j,k)*sdcf
      END DO
   END DO
END DO

DEALLOCATE(tv)

RETURN
END SUBROUTINE rsvgthr
