SUBROUTINE corrscat(v,sp,rsv, vit)

!------------------------------------------------------------------
!
! Coarse grid processes send correction to fine grid processes.
!
!------------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER, INTENT(IN) :: v, sp, vit
INTEGER :: vv, spp, xt, yt, zt, pk, pj, pi, kt, jt, it, k, j, i, ieq, jeq
INTEGER :: trank, tag, ierr
INTEGER, DIMENSION(3) :: coord, istat
REAL*8, DIMENSION(cneq), INTENT(OUT) :: rsv
REAL*8, DIMENSION(cneq) :: tsv

INCLUDE 'mpif.h'

vv  = v+1
spp = 2**(vv-1)
xt  = MOD(xrank,spp)
yt  = MOD(yrank,spp)
zt  = MOD(zrank,spp)

! If part of the previous stage, need to send
IF (xt == 0 .AND. yt == 0 .AND. zt == 0) THEN
   DO pk = zrank, zrank+sp, sp
      coord(1) = pk
      kt = (pk - zrank)*nz/(sp*2)
      DO pj = yrank, yrank+sp, sp
         coord(2) = pj
         jt = (pj - yrank)*ny/(sp*2)
         DO pi = xrank, xrank+sp, sp
            coord(3) = pi
            it = (pi - xrank)*nx/(sp*2)
            CALL MPI_CART_RANK(allcomm,coord,trank,ierr)
            IF (trank /= nrank) THEN
               DO k = kt + 1, kt + nz/2
                  DO j = jt + 1, jt + ny/2
                     DO i = it + 1, it + nx/2
                        ieq = (k-1)*xys + (j-1)*nx + i
                        jeq = (k-kt-1)*xys/4 + (j-jt-1)*nx/2 + i-it
                        tsv(jeq) = phi(ieq,vv)
                     END DO
                  END DO
               END DO

!if (vit == 3 .and. v==1)  print *, nrank, trank, rsv

               tag = 300 + trank
               CALL MPI_SEND(tsv,cneq,MPI_DOUBLE_PRECISION,trank,tag,allcomm,ierr)
            ELSE
               DO k = 1, nz/2
                  DO j = 1, ny/2
                     DO i = 1, nx/2
                        ieq = (k-1)*xys + (j-1)*nx + i
                        jeq = (k-1)*xys/4 + (j-1)*nx/2 + i
                        rsv(jeq) = phi(ieq,vv)
                     END DO
                  END DO
               END DO
            END IF


!if (nrank == 0 .and. v==1) print*, "rsv", vit, trank, rsv



         END DO
      END DO
   END DO
ELSE

   ! Fine grid processes receive the coarse grid's correction
   ! Determine where it came from
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

   tag = 300 + nrank
   CALL MPI_RECV(rsv,cneq,MPI_DOUBLE_PRECISION,trank,tag,allcomm,istat,ierr)
END IF

RETURN
END SUBROUTINE corrscat
