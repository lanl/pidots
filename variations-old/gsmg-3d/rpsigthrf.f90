SUBROUTINE rpsigthrf(v,spx,spy,spz,kc,jc,ic,cbcs,crsf,sdcf,rpave,fav)

!-------------------------------------------------------------
!
! Gathers all the rpave vectors and puts them into av 
! 
!-------------------------------------------------------------

USE invar
IMPLICIT NONE
INTEGER, INTENT(IN) :: v, spx, spy, spz, kc, jc, ic, cbcs, crsf
INTEGER :: spx2, spy2, spz2, i, j, k, xt, yt, zt, it, jt, kt, n, o, ieq, jeq
INTEGER :: ix0, jx0, ix1, jx1, ix2, jx2, kaa, jaa, iaa, ijc, ikc, jkc
INTEGER :: tz, ty, tx, txy, txz, tyz, tbcs, kcc, jcc, icc
INTEGER :: pk, pj, pi, trank, tag, ierr
INTEGER, DIMENSION(3) :: coord, istat
REAL*8, INTENT(IN) :: sdcf
REAL*8 :: xyf, xzf, yzf
REAL*8, DIMENSION(cbcs,8) :: tmpsi
REAL*8, DIMENSION(cbcs,8), INTENT(IN) :: rpave
REAL*8, DIMENSION(bcs,8), INTENT(OUT) :: fav
REAL*8, DIMENSION(:,:), ALLOCATABLE :: tv

INCLUDE 'mpif.h'

! Need zcf/ycf/xcf for communication
kcc = zcf(v)
jcc = ycf(v)
icc = xcf(v)
kaa = kcc - 1
jaa = jcc - 1
iaa = icc - 1

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

! Setup for parallel
spx2 = spx*icc
spy2 = spy*jcc
spz2 = spz*kcc
xt = MOD(xrank,spx2)
yt = MOD(yrank,spy2)
zt = MOD(zrank,spz2)

! If not part of the next stage, need to send
IF (xt /= 0 .OR. yt /= 0 .OR. zt /= 0) THEN
   ! Determine recipient with coordinates
   coord(1) = zrank - zt
   coord(2) = yrank - yt
   coord(3) = xrank - xt
   CALL MPI_CART_RANK(allcomm,coord,trank,ierr)
   tag = 200 + nrank
   CALL MPI_SEND(rpave,cbcs*8,MPI_DOUBLE_PRECISION,trank,tag,allcomm,ierr)

! If part of the next stage, need to receive
ELSE
   DO pk = zrank, zrank+kaa*spz, spz
      coord(1) = pk
      kt = (pk - zrank)*sdnz/(spz*kc)
      DO pj = yrank, yrank+jaa*spy, spy
         coord(2) = pj
         jt = (pj - yrank)*sdny/(spy*jc)
         DO pi = xrank, xrank+iaa*spx, spx
            coord(3) = pi
            it = (pi - xrank)*sdnx/(spx*ic)
            CALL MPI_CART_RANK(allcomm,coord,trank,ierr)
            IF (trank /= nrank) THEN
               ! Receive from another processor and copy into av of the next stage
               tag = 200 + trank
               CALL MPI_RECV(tmpsi,cbcs*8,MPI_DOUBLE_PRECISION,trank,tag,allcomm,istat,ierr)
               IF (pk == zrank+kaa*spz) THEN
                  DO n = 1, apo
                     ix1 = (n-1)*xys/ijc
                     jx1 = (n-1)*txy
                     DO j = jt + 1, jt + sdny/jc
                        ix2 = ix1 + (j-jt-1)*sdnx/ic
                        jx2 = jx1 + (j-1)*tx
                        DO i = it + 1, it + sdnx/ic
                           ieq = ix2 + i-it
                           jeq = jx2 + i
                           tv(jeq,1) = tmpsi(ieq,1)
                           tv(jeq,2) = tmpsi(ieq,2)
                           tv(jeq,3) = tmpsi(ieq,3)
                           tv(jeq,4) = tmpsi(ieq,4)
                        END DO
                     END DO
                  END DO
               ELSE IF (pk == zrank) THEN
                  DO n = 1, apo
                     ix1 = (n-1)*xys/ijc
                     jx1 = (n-1)*txy
                     DO j = jt + 1, jt + sdny/jc
                        ix2 = ix1 + (j-jt-1)*sdnx/ic
                        jx2 = jx1 + (j-1)*tx
                        DO i = it + 1, it + sdnx/ic
                           ieq = ix2 + i-it
                           jeq = jx2 + i
                           tv(jeq,5) = tmpsi(ieq,5)
                           tv(jeq,6) = tmpsi(ieq,6)
                           tv(jeq,7) = tmpsi(ieq,7)
                           tv(jeq,8) = tmpsi(ieq,8)
                        END DO
                     END DO
                  END DO
               END IF
               IF (pj == yrank+jaa*spy) THEN
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
                           tv(jeq,1) = tmpsi(ieq,1)
                           tv(jeq,2) = tmpsi(ieq,2)
                           tv(jeq,5) = tmpsi(ieq,5)
                           tv(jeq,6) = tmpsi(ieq,6)
                        END DO
                     END DO
                  END DO
               ELSE IF (pj == yrank) THEN
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
                           tv(jeq,3) = tmpsi(ieq,3)
                           tv(jeq,4) = tmpsi(ieq,4)
                           tv(jeq,7) = tmpsi(ieq,7)
                           tv(jeq,8) = tmpsi(ieq,8)
                        END DO
                     END DO
                  END DO
               END IF
               IF (pi == xrank+iaa*spx) THEN
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
                           tv(jeq,1) = tmpsi(ieq,1)
                           tv(jeq,4) = tmpsi(ieq,4)
                           tv(jeq,5) = tmpsi(ieq,5)
                           tv(jeq,8) = tmpsi(ieq,8)
                        END DO
                     END DO
                  END DO
               ELSE IF (pi == xrank) THEN
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
                           tv(jeq,2) = tmpsi(ieq,2)
                           tv(jeq,3) = tmpsi(ieq,3)
                           tv(jeq,6) = tmpsi(ieq,6)
                           tv(jeq,7) = tmpsi(ieq,7)
                        END DO
                     END DO
                  END DO
               END IF
            ELSE
               ! Know that it=jt=kt=0
               DO n = 1, apo
                  ix1 = (n-1)*xys/ijc
                  jx1 = (n-1)*txy
                  DO j = 1, sdny/jc
                     ix2 = ix1 + (j-1)*sdnx/ic
                     jx2 = jx1 + (j-1)*tx
                     DO i = 1, sdnx/ic
                        ieq = ix2 + i
                        jeq = jx2 + i
                        tv(jeq,5) = rpave(ieq,5)
                        tv(jeq,6) = rpave(ieq,6)
                        tv(jeq,7) = rpave(ieq,7)
                        tv(jeq,8) = rpave(ieq,8)
                     END DO
                  END DO
               END DO
               ix0 = apo*xys/ijc
               jx0 = apo*txy
               DO n = 1, apo
                  ix1 = ix0 + (n-1)*xzs/ikc
                  jx1 = jx0 + (n-1)*txz
                  DO k = 1, sdnz/kc
                     ix2 = ix1 + (k-1)*sdnx/ic
                     jx2 = jx1 + (k-1)*tx
                     DO i = 1, sdnx/ic
                        ieq = ix2 + i
                        jeq = jx2 + i
                        tv(jeq,3) = rpave(ieq,3)
                        tv(jeq,4) = rpave(ieq,4)
                        tv(jeq,7) = rpave(ieq,7)
                        tv(jeq,8) = rpave(ieq,8)
                     END DO
                  END DO
               END DO
               ix0 = apo*(xys/ijc + xzs/ikc)
               jx0 = apo*(txy + txz)
               DO n = 1, apo
                  ix1 = ix0 + (n-1)*yzs/jkc
                  jx1 = jx0 + (n-1)*tyz
                  DO k = 1, sdnz/kc
                     ix2 = ix1 + (k-1)*sdny/jc
                     jx2 = jx1 + (k-1)*ty
                     DO j = 1, sdny/jc
                        ieq = ix2 + j
                        jeq = jx2 + j
                        tv(jeq,2) = rpave(ieq,2)
                        tv(jeq,3) = rpave(ieq,3)
                        tv(jeq,6) = rpave(ieq,6)
                        tv(jeq,7) = rpave(ieq,7)
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
END IF

DEALLOCATE(tv)

RETURN
END SUBROUTINE rpsigthrf
