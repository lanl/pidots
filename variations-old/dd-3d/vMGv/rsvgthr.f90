SUBROUTINE rsvgthr(v,sp,rsv)

!-------------------------------------------------------------
!
! Gathers rsv vectors and puts solution into sv(:,v+1)
! 
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER, INTENT(IN) :: v, sp
INTEGER :: vv, spp, tag, pk, pj, pi, trank, ieq, jeq, ierr
INTEGER :: xt, yt, zt, kt, jt, it, k, j, i
INTEGER, DIMENSION(3) :: coord, istat
REAL*8, DIMENSION(cneq) :: tmpsv
REAL*8, DIMENSION(cneq), INTENT(IN) :: rsv

INCLUDE 'mpif.h'

! Setup for parallel
vv = v+1
spp = 2**(vv-1)
xt = MOD(xrank,spp)
yt = MOD(yrank,spp)
zt = MOD(zrank,spp)

! If not part of next stage, need to send
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
   tag = 100 + nrank
   CALL MPI_SEND(rsv,cneq,MPI_DOUBLE_PRECISION,trank,tag,allcomm,ierr)

!print *, xrank, yrank, zrank, nrank, "rsv sent to", trank


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
               ! Receive from another processor and copy into sv of the next stage
               tag = 100 + trank
               CALL MPI_RECV(tmpsv,cneq,MPI_DOUBLE_PRECISION,trank,tag,allcomm,istat,ierr)


!print *, xrank, yrank, zrank, nrank, "rsv received"

               DO k = 1, nz/2
                  DO j = 1, ny/2
                     DO i = 1, nx/2
                        ieq = (k-1)*xys/4 + (j-1)*nx/2 + i
                        jeq = kt*neq/2 + (k-1)*xys + jt*xys/2 + (j-1)*nx + it*nx/2 + i
                        sv(jeq,vv) = tmpsv(ieq)
                     END DO
                  END DO
               END DO
            ELSE
               ! Process going to next stage copies its own rsv into sv
               DO k = 1, nz/2
                  DO j = 1, ny/2
                     DO i = 1, nx/2
                        ieq = (k-1)*xys/4 + (j-1)*nx/2 + i
                        jeq = (k-1)*xys + (j-1)*nx + i
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
