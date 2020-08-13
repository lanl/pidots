SUBROUTINE rsvgthr(v,spx,spy,spz,kc,jc,ic,cneq,rsv)

!-------------------------------------------------------------
!
! Gathers rsv vectors and puts solution into sv(:,v+1)
! 
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER, INTENT(IN) :: v, spx, spy, spz, kc, jc, ic, cneq
INTEGER :: vv, spx2, spy2, spz2, tag, pk, pj, pi, trank, ieq, jeq, ierr
INTEGER :: xt, yt, zt, kt, jt, it, k, j, i, ix1, jx1, ix2, jx2
INTEGER :: ka, ja, ia, ijc
INTEGER, DIMENSION(3) :: coord, istat
REAL*8, DIMENSION(cneq) :: tmpsv
REAL*8, DIMENSION(cneq), INTENT(IN) :: rsv

INCLUDE 'mpif.h'

! Use kc/jc/ic to determine other coarsening parameters
ka = kc - 1
ja = jc - 1
ia = ic - 1
ijc = ic*jc

! Setup for parallel
vv = v+1
spx2 = spx*ic
spy2 = spy*jc
spz2 = spz*kc
xt = MOD(xrank,spx2)
yt = MOD(yrank,spy2)
zt = MOD(zrank,spz2)

! If not part of next stage, need to send
IF (xt /= 0 .OR. yt /= 0 .OR. zt /= 0) THEN
   ! Determine recipient with coordinates
   coord(1) = zrank - zt
   coord(2) = yrank - yt
   coord(3) = xrank - xt
   CALL MPI_CART_RANK(allcomm,coord,trank,ierr)
   tag = 100 + nrank
   CALL MPI_SEND(rsv,cneq,MPI_DOUBLE_PRECISION,trank,tag,allcomm,ierr)

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
               ! Receive from another processor and copy into sv of the next stage
               tag = 100 + trank
               CALL MPI_RECV(tmpsv,cneq,MPI_DOUBLE_PRECISION,trank,tag,allcomm,istat,ierr)
               DO k = kt + 1, kt + nz/kc
                  ix1 = (k-kt-1)*xys/ijc
                  jx1 = (k-1)*xys
                  DO j = jt + 1, jt + ny/jc
                     ix2 = ix1 + (j-jt-1)*nx/ic
                     jx2 = jx1 + (j-1)*nx
                     DO i = it + 1, it + nx/ic
                        ieq = ix2 + i-it
                        jeq = jx2 + i
                        sv(jeq,vv) = tmpsv(ieq)
                     END DO
                  END DO
               END DO
            ELSE
               ! Process going to next stage copies its own rsv into sv
               DO k = 1, nz/kc
                  ix1 = (k-1)*xys/ijc
                  jx1 = (k-1)*xys
                  DO j = 1, ny/jc
                     ix2 = ix1 + (j-1)*nx/ic
                     jx2 = jx1 + (j-1)*nx
                     DO i = 1, nx/ic
                        ieq = ix2 + i
                        jeq = jx2 + i
                        sv(jeq,vv) = rsv(ieq)
                     END DO
                  END DO
               END DO
            END IF
         END DO
      END DO
   END DO
END IF

RETURN
END SUBROUTINE rsvgthr
