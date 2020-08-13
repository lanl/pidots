SUBROUTINE rpsigthr(v,sp,rpave)

!-------------------------------------------------------------
!
! Gathers all the rpave vectors and puts them into rrpsi 
! 
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER, INTENT(IN) :: v, sp
INTEGER :: vv, spp, i, j, k, xt, yt, zt, it, jt, kt, n, ieq, jeq
INTEGER :: pk, pj, pi, trank, tag, ierr
INTEGER, DIMENSION(3) :: coord, istat
REAL*8, DIMENSION(cbcs,8) :: tmpsi
REAL*8, DIMENSION(cbcs,8), INTENT(IN) :: rpave

INCLUDE 'mpif.h'

! Setup for parallel
vv = v+1
spp = 2**(vv-1)
xt = MOD(xrank,spp)
yt = MOD(yrank,spp)
zt = MOD(zrank,spp)

! If not part of the next stage, need to send
IF (xt /= 0 .OR. yt /= 0 .OR. zt /= 0) THEN
   ! Determine recipient with coordinates
   IF (zt /= 0) THEN
      coord(1) = zrank - sp
   ELSE
      coord(1) = zrank
   END IF
   IF (yt /= 0) THEN
      coord(2) = yrank - sp
   ELSE
      coord(2) = yrank
   END IF
   IF (xt /= 0) THEN
      coord(3) = xrank - sp
   ELSE
      coord(3) = xrank
   END IF
   CALL MPI_CART_RANK(allcomm,coord,trank,ierr)
   tag = 200 + nrank
   CALL MPI_SEND(rpave,cbcs*8,MPI_DOUBLE_PRECISION,trank,tag,allcomm,ierr)

! If part of the next stage, need to receive
ELSE
   DO pk = zrank, zrank+sp, sp
      coord(1) = pk
      kt = (pk - zrank)/sp
      DO pj = yrank, yrank+sp, sp
         coord(2) = pj
         jt = (pj - yrank)/sp
         DO pi = xrank, xrank+sp, sp
            coord(3) = pi
            it = (pi - xrank)/sp
            CALL MPI_CART_RANK(allcomm,coord,trank,ierr)
            IF (trank /= nrank) THEN
               ! Receive from another processor and copy into av of the next stage
               tag = 200 + trank
               CALL MPI_RECV(tmpsi,cbcs*8,MPI_DOUBLE_PRECISION,trank,tag,allcomm,istat,ierr)
               IF (kt == 1) THEN
                  DO n = 1, apo
                     DO j = 1, ny/2
                        DO i = 1, nx/2
                           ieq = (n-1)*xys/4 + (j-1)*nx/2 + i
                           jeq = (n-1)*xys + jt*xys/2 + (j-1)*nx + it*nx/2 + i
                           av(jeq,1,vv) = tmpsi(ieq,1)
                           av(jeq,2,vv) = tmpsi(ieq,2)
                           av(jeq,3,vv) = tmpsi(ieq,3)
                           av(jeq,4,vv) = tmpsi(ieq,4)
                        END DO
                     END DO
                  END DO
               ELSE
                  DO n = 1, apo
                     DO j = 1, ny/2
                        DO i = 1, nx/2
                           ieq = (n-1)*xys/4 + (j-1)*nx/2 + i
                           jeq = (n-1)*xys + jt*xys/2 + (j-1)*nx + it*nx/2 + i
                           av(jeq,5,vv) = tmpsi(ieq,5)
                           av(jeq,6,vv) = tmpsi(ieq,6)
                           av(jeq,7,vv) = tmpsi(ieq,7)
                           av(jeq,8,vv) = tmpsi(ieq,8)
                        END DO
                     END DO
                  END DO
               END IF
               IF (jt == 1) THEN
                  DO n = 1, apo
                     DO k = 1, nz/2
                        DO i = 1, nx/2
                           ieq = apo*xys/4 + (n-1)*xzs/4 + (k-1)*nx/2 + i
                           jeq = apo*xys + (n-1)*xzs + kt*xzs/2 + (k-1)*nx + it*nx/2 + i
                           av(jeq,1,vv) = tmpsi(ieq,1)
                           av(jeq,2,vv) = tmpsi(ieq,2)
                           av(jeq,5,vv) = tmpsi(ieq,5)
                           av(jeq,6,vv) = tmpsi(ieq,6)
                        END DO
                     END DO
                  END DO
               ELSE
                  DO n = 1, apo
                     DO k = 1, nz/2
                        DO i = 1, nx/2
                           ieq = apo*xys/4 + (n-1)*xzs/4 + (k-1)*nx/2 + i
                           jeq = apo*xys + (n-1)*xzs + kt*xzs/2 + (k-1)*nx + it*nx/2 + i
                           av(jeq,3,vv) = tmpsi(ieq,3)
                           av(jeq,4,vv) = tmpsi(ieq,4)
                           av(jeq,7,vv) = tmpsi(ieq,7)
                           av(jeq,8,vv) = tmpsi(ieq,8)
                        END DO
                     END DO
                  END DO
               END IF
               IF (it == 1) THEN
                  DO n = 1, apo
                     DO k = 1, nz/2
                        DO j = 1, ny/2
                           ieq = apo*(xys+xzs)/4 + (n-1)*yzs/4 + (k-1)*ny/2 + j
                           jeq = apo*(xys+xzs) + (n-1)*yzs + kt*yzs/2 + (k-1)*ny + jt*ny/2 + j
                           av(jeq,1,vv) = tmpsi(ieq,1)
                           av(jeq,4,vv) = tmpsi(ieq,4)
                           av(jeq,5,vv) = tmpsi(ieq,5)
                           av(jeq,8,vv) = tmpsi(ieq,8)
                        END DO
                     END DO
                  END DO
               ELSE
                  DO n = 1, apo
                     DO k = 1, nz/2
                        DO j = 1, ny/2
                           ieq = apo*(xys+xzs)/4 + (n-1)*yzs/4 + (k-1)*ny/2 + j
                           jeq = apo*(xys+xzs) + (n-1)*yzs + kt*yzs/2 + (k-1)*ny + jt*ny/2 + j
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
                  DO j = 1, ny/2
                     DO i = 1, nx/2
                        ieq = (n-1)*xys/4 + (j-1)*nx/2 + i
                        jeq = (n-1)*xys + (j-1)*nx + i
                        av(jeq,5,vv) = rpave(ieq,5)
                        av(jeq,6,vv) = rpave(ieq,6)
                        av(jeq,7,vv) = rpave(ieq,7)
                        av(jeq,8,vv) = rpave(ieq,8)
                     END DO
                  END DO
               END DO
               DO n = 1, apo
                  DO k = 1, nz/2
                     DO i = 1, nx/2
                        ieq = apo*xys/4 + (n-1)*xzs/4 + (k-1)*nx/2 + i
                        jeq = apo*xys + (n-1)*xzs + (k-1)*nx + i
                        av(jeq,3,vv) = rpave(ieq,3)
                        av(jeq,4,vv) = rpave(ieq,4)
                        av(jeq,7,vv) = rpave(ieq,7)
                        av(jeq,8,vv) = rpave(ieq,8)
                     END DO
                  END DO
               END DO
               DO n = 1, apo
                  DO k = 1, nz/2
                     DO j = 1, ny/2
                        ieq = apo*(xys+xzs)/4 + (n-1)*yzs/4 + (k-1)*ny/2 + j
                        jeq = apo*(xys+xzs) + (n-1)*yzs + (k-1)*ny + j
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
