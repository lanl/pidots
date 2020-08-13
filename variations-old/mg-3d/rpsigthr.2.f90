SUBROUTINE rpsigthr(v,spx,spy,spz,kc,jc,ic,cbcs,rpave)

!-------------------------------------------------------------
!
! Gathers all the rpave vectors and puts them into av 
! 
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER, INTENT(IN) :: v, spx, spy, spz, kc, jc, ic, cbcs
INTEGER :: vv, spx2, spy2, spz2, i, j, k, xt, yt, zt, it, jt, kt, n, ieq, jeq
INTEGER :: ix0, jx0, ix1, jx1, ix2, jx2, ka, ja, ia, ijc, ikc, jkc
INTEGER :: pk, pj, pi, trank, tag, ierr
INTEGER, DIMENSION(3) :: coord, istat
REAL*8, DIMENSION(cbcs,8) :: tmpsi
REAL*8, DIMENSION(cbcs,8), INTENT(IN) :: rpave

INCLUDE 'mpif.h'

! Use kc/jc/ic to determine other coarsening parameters
ka = kc - 1
ja = jc - 1
ia = ic - 1
ijc = ic*jc
ikc = ic*kc
jkc = jc*kc

! Setup for parallel
vv = v+1
spx2 = spx*ic
spy2 = spy*jc
spz2 = spz*kc
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
   DO pk = zrank, zrank+ka*spz, spz
      coord(1) = pk
      kt = (pk - zrank)*nz/(spz*kc)
      DO pj = yrank, yrank+ja*spy, spy
         coord(2) = pj
         jt = (pj - yrank)*ny/(spy*jc)
         DO pi = xrank, xrank+ia*spx, spx
            coord(3) = pi
            it = (pi - xrank)*nx/(spx*ic)
            CALL MPI_CART_RANK(allcomm,coord,trank,ierr)
            IF (trank /= nrank) THEN
               ! Receive from another processor and copy into av of the next stage
               tag = 200 + trank
               CALL MPI_RECV(tmpsi,cbcs*8,MPI_DOUBLE_PRECISION,trank,tag,allcomm,istat,ierr)
               IF (pk == zrank+ka*spz) THEN
                  DO n = 1, apo
                     ix1 = (n-1)*xys/ijc
                     jx1 = (n-1)*xys
                     DO j = jt + 1, jt + ny/jc
                        ix2 = ix1 + (j-jt-1)*nx/ic
                        jx2 = jx1 + (j-1)*nx
                        DO i = it + 1, it + nx/ic
                           ieq = ix2 + i-it
                           jeq = jx2 + i
                           av(jeq,1,vv) = tmpsi(ieq,1)
                           av(jeq,2,vv) = tmpsi(ieq,2)
                           av(jeq,3,vv) = tmpsi(ieq,3)
                           av(jeq,4,vv) = tmpsi(ieq,4)
                        END DO
                     END DO
                  END DO
               ELSE IF (pk == zrank) THEN
                  DO n = 1, apo
                     ix1 = (n-1)*xys/ijc
                     jx1 = (n-1)*xys
                     DO j = jt + 1, jt + ny/jc
                        ix2 = ix1 + (j-jt-1)*nx/ic
                        jx2 = jx1 + (j-1)*nx
                        DO i = it + 1, it + nx/ic
                           ieq = ix2 + i-it
                           jeq = jx2 + i
                           av(jeq,5,vv) = tmpsi(ieq,5)
                           av(jeq,6,vv) = tmpsi(ieq,6)
                           av(jeq,7,vv) = tmpsi(ieq,7)
                           av(jeq,8,vv) = tmpsi(ieq,8)
                        END DO
                     END DO
                  END DO
               END IF
               IF (pj == yrank+ja*spy) THEN
                  ix0 = apo*xys/ijc
                  jx0 = apo*xys
                  DO n = 1, apo
                     ix1 = ix0 + (n-1)*xzs/ikc
                     jx1 = jx0 + (n-1)*xzs
                     DO k = kt + 1, kt + nz/kc
                        ix2 = ix1 + (k-kt-1)*nx/ic
                        jx2 = jx1 + (k-1)*nx
                        DO i = it + 1, it + nx/ic
                           ieq = ix2 + i-it
                           jeq = jx2 + i
                           av(jeq,1,vv) = tmpsi(ieq,1)
                           av(jeq,2,vv) = tmpsi(ieq,2)
                           av(jeq,5,vv) = tmpsi(ieq,5)
                           av(jeq,6,vv) = tmpsi(ieq,6)
                        END DO
                     END DO
                  END DO
               ELSE IF (pj == yrank) THEN
                  ix0 = apo*xys/ijc
                  jx0 = apo*xys
                  DO n = 1, apo
                     ix1 = ix0 + (n-1)*xzs/ikc
                     jx1 = jx0 + (n-1)*xzs
                     DO k = kt + 1, kt + nz/kc
                        ix2 = ix1 + (k-kt-1)*nx/ic
                        jx2 = jx1 + (k-1)*nx
                        DO i = it + 1, it + nx/ic
                           ieq = ix2 + i-it
                           jeq = jx2 + i
                           av(jeq,3,vv) = tmpsi(ieq,3)
                           av(jeq,4,vv) = tmpsi(ieq,4)
                           av(jeq,7,vv) = tmpsi(ieq,7)
                           av(jeq,8,vv) = tmpsi(ieq,8)
                        END DO
                     END DO
                  END DO
               END IF
               IF (pi == xrank+ia*spx) THEN
                  ix0 = apo*(xys/ijc + xzs/ikc)
                  jx0 = apo*(xys + xzs)
                  DO n = 1, apo
                     ix1 = ix0 + (n-1)*yzs/jkc
                     jx1 = jx0 + (n-1)*yzs
                     DO k = kt + 1, kt + nz/kc
                        ix2 = ix1 + (k-kt-1)*ny/jc
                        jx2 = jx1 + (k-1)*ny
                        DO j = jt + 1, jt + ny/jc
                           ieq = ix2 + j-jt
                           jeq = jx2 + j
                           av(jeq,1,vv) = tmpsi(ieq,1)
                           av(jeq,4,vv) = tmpsi(ieq,4)
                           av(jeq,5,vv) = tmpsi(ieq,5)
                           av(jeq,8,vv) = tmpsi(ieq,8)
                        END DO
                     END DO
                  END DO
               ELSE IF (pi == xrank) THEN
                  ix0 = apo*(xys/ijc + xzs/ikc)
                  jx0 = apo*(xys + xzs)
                  DO n = 1, apo
                     ix1 = ix0 + (n-1)*yzs/jkc
                     jx1 = jx0 + (n-1)*yzs
                     DO k = kt + 1, kt + nz/kc
                        ix2 = ix1 + (k-kt-1)*ny/jc
                        jx2 = jx1 + (k-1)*ny
                        DO j = jt + 1, jt + ny/jc
                           ieq = ix2 + j-jt
                           jeq = jx2 + j
                           av(jeq,2,vv) = tmpsi(ieq,2)
                           av(jeq,3,vv) = tmpsi(ieq,3)
                           av(jeq,6,vv) = tmpsi(ieq,6)
                           av(jeq,7,vv) = tmpsi(ieq,7)
                        END DO
                     END DO
                  END DO
               END IF
            ELSE
               ! Know that it=jt=kt=0
               DO n = 1, apo
                  ix1 = (n-1)*xys/ijc
                  jx1 = (n-1)*xys
                  DO j = 1, ny/jc
                     ix2 = ix1 + (j-1)*nx/ic
                     jx2 = jx1 + (j-1)*nx
                     DO i = 1, nx/ic
                        ieq = ix2 + i
                        jeq = jx2 + i
                        av(jeq,5,vv) = rpave(ieq,5)
                        av(jeq,6,vv) = rpave(ieq,6)
                        av(jeq,7,vv) = rpave(ieq,7)
                        av(jeq,8,vv) = rpave(ieq,8)
                     END DO
                  END DO
               END DO
               ix0 = apo*xys/ijc
               jx0 = apo*xys
               DO n = 1, apo
                  ix1 = ix0 + (n-1)*xzs/ikc
                  jx1 = jx0 + (n-1)*xzs
                  DO k = 1, nz/kc
                     ix2 = ix1 + (k-1)*nx/ic
                     jx2 = jx1 + (k-1)*nx
                     DO i = 1, nx/ic
                        ieq = ix2 + i
                        jeq = jx2 + i
                        av(jeq,3,vv) = rpave(ieq,3)
                        av(jeq,4,vv) = rpave(ieq,4)
                        av(jeq,7,vv) = rpave(ieq,7)
                        av(jeq,8,vv) = rpave(ieq,8)
                     END DO
                  END DO
               END DO
               ix0 = apo*(xys/ijc + xzs/ikc)
               jx0 = apo*(xys + xzs)
               DO n = 1, apo
                  ix1 = ix0 + (n-1)*yzs/jkc
                  jx1 = jx0 + (n-1)*yzs
                  DO k = 1, nz/kc
                     ix2 = ix1 + (k-1)*ny/jc
                     jx2 = jx1 + (k-1)*ny
                     DO j = 1, ny/jc
                        ieq = ix2 + j
                        jeq = jx2 + j
                        av(jeq,2,vv) = rpave(ieq,2)
                        av(jeq,3,vv) = rpave(ieq,3)
                        av(jeq,6,vv) = rpave(ieq,6)
                        av(jeq,7,vv) = rpave(ieq,7)
                     END DO
                  END DO
               END DO
            END IF
         END DO
      END DO
   END DO
END IF

RETURN
END SUBROUTINE rpsigthr
