SUBROUTINE rsvgthr(v,spx,spy,spz,kc,jc,ic,cneq,crsf,sdcf,rsv)

!-------------------------------------------------------------
!
! Gathers rsv vectors and puts solution into sv(:,v+1)
! 
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER, INTENT(IN) :: v, spx, spy, spz, kc, jc, ic, cneq, crsf
INTEGER :: vv, spx2, spy2, spz2, tag, pk, pj, pi, trank, ieq, jeq, ierr
INTEGER :: xt, yt, zt, kt, jt, it, k, j, i, ix1, jx1, ix2, jx2
INTEGER :: kcc, jcc, icc, kaa, jaa, iaa, ijc, tz, ty, tx
INTEGER, DIMENSION(3) :: coord, istat
REAL*8, INTENT(IN) :: sdcf
REAL*8, DIMENSION(cneq) :: tmpsv
REAL*8, DIMENSION(cneq), INTENT(IN) :: rsv
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: tv

INCLUDE 'mpif.h'

! Need zcf/ycf/xcf for communication
kcc = zcf(v)
jcc = ycf(v)
icc = xcf(v)
kaa = kcc - 1
jaa = jcc - 1
iaa = icc - 1

! Set up the temporary matrix
tz = nz
ty = ny
tx = nx
IF (kcc > nz) tz = kcc
IF (jcc > ny) ty = jcc
IF (icc > nx) tx = icc
ALLOCATE(tv(tx,ty,tz))
ijc = ic*jc

! Setup for parallel
vv = v+1
spx2 = spx*icc
spy2 = spy*jcc
spz2 = spz*kcc
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
   DO pk = zrank, zrank+kaa*spz, spz
      coord(1) = pk
      kt = (pk - zrank)*nz/(spz*kc)
      DO pj = yrank, yrank+jaa*spy, spy
         coord(2) = pj
         jt = (pj - yrank)*ny/(spy*jc)
         DO pi = xrank, xrank+iaa*spx, spx
            coord(3) = pi
            it = (pi - xrank)*nx/(spx*ic)
            CALL MPI_CART_RANK(allcomm,coord,trank,ierr)
            IF (trank /= nrank) THEN
               ! Receive from another processor and copy into sv of the next stage
               tag = 100 + trank
               CALL MPI_RECV(tmpsv,cneq,MPI_DOUBLE_PRECISION,trank,tag,allcomm,istat,ierr)
               DO k = kt + 1, kt + nz/kc
                  ix1 = (k-kt-1)*xys/ijc
                  DO j = jt + 1, jt + ny/jc
                     ix2 = ix1 + (j-jt-1)*nx/ic
                     DO i = it + 1, it + nx/ic
                        ieq = ix2 + i-it
                        tv(i,j,k) = tmpsv(ieq)
                     END DO
                  END DO
               END DO
            ELSE
               DO k = 1, nz/kc
                  ix1 = (k-1)*xys/ijc
                  DO j = 1, ny/jc
                     ix2 = ix1 + (j-1)*nx/ic
                     DO i = 1, nx/ic
                        ieq = ix2 + i
                        tv(i,j,k) = rsv(ieq)
                     END DO
                  END DO
               END DO
            END IF
         END DO
      END DO
   END DO
   ! Reinitialize sv
   sv(:,vv) = 0.0
   ! Now compute average for combined cells/sub-domains
   DO k = 1, tz
      ix1 = (k-1)*xys
      IF (kcc > nz) ix1 = ((k-1+kcc/nz)*nz/kcc-1)*xys
      DO j = 1, ty
         ix2 = ix1 + (j-1)*nx
         IF (jcc > ny) ix2 = ix1 + ((j-1+jcc/ny)*ny/jcc-1)*nx
         DO i = 1, tx
            ieq = ix2 + i
            IF (icc > nx) ieq = ix2 + (i-1+icc/nx)*nx/icc
            sv(ieq,vv) = sv(ieq,vv) + tv(i,j,k)*sdcf
         END DO
      END DO
   END DO
END IF

DEALLOCATE(tv)

RETURN
END SUBROUTINE rsvgthr
