SUBROUTINE corrscat(v,spx,spy,spz,kc,jc,ic,cneq,rsv)

!------------------------------------------------------------------
!
! Coarse grid processes send correction to fine grid processes.
!
!------------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER, INTENT(IN) :: v, spx, spy, spz, kc, jc, ic, cneq
INTEGER :: vv, spx2, spy2, spz2, xt, yt, zt, pk, pj, pi, kt, jt, it, k, j, i, ieq
INTEGER :: kaa, jaa, iaa, ijc, ix1, ix2, trank, tag, ierr
INTEGER :: tz, ty, tx, kcc, jcc, icc
INTEGER, DIMENSION(3) :: coord, istat
REAL*8, DIMENSION(cneq), INTENT(OUT) :: rsv
REAL*8, DIMENSION(cneq) :: tsv
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

vv  = v+1
spx2 = spx*icc
spy2 = spy*jcc
spz2 = spz*kcc
xt  = MOD(xrank,spx2)
yt  = MOD(yrank,spy2)
zt  = MOD(zrank,spz2)

! If part of the previous stage, need to send
IF (xt == 0 .AND. yt == 0 .AND. zt == 0) THEN
   ! Copy phi multiplied by phir2 to tv. Will then send tv elements to pack in tsv
   DO k = 1, tz
      ix1 = (k-1)*xys
      IF (kcc > nz) ix1 = ((k-1+kcc/nz)*nz/kcc-1)*xys
      DO j = 1, ty
         ix2 = ix1 + (j-1)*nx
         IF (jcc > ny) ix2 = ix1 + ((j-1+jcc/ny)*ny/jcc-1)*nx
         DO i = 1, tx
            ieq = ix2 + i
            IF (icc > nx) ieq = ix2 + (i-1+icc/nx)*nx/icc
            tv(i,j,k) = phi(ieq,vv)      !phir2(i,j,k,vv)*phi(ieq,vv)
         END DO
      END DO
   END DO
   ! Now copy tv into tsv and send
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
               DO k = kt + 1, kt + nz/kc
                  ix1 = (k-kt-1)*xys/ijc
                  DO j = jt + 1, jt + ny/jc
                     ix2 = ix1 + (j-jt-1)*nx/ic
                     DO i = it + 1, it + nx/ic
                        ieq = ix2 + i-it
                        tsv(ieq) = tv(i,j,k)
                     END DO
                  END DO
               END DO
               tag = 300 + trank
               CALL MPI_SEND(tsv,cneq,MPI_DOUBLE_PRECISION,trank,tag,allcomm,ierr)
            ELSE
               DO k = 1, nz/kc
                  ix1 = (k-1)*xys/ijc
                  DO j = 1, ny/jc
                     ix2 = ix1 + (j-1)*nx/ic
                     DO i = 1, nx/ic
                        ieq = ix2 + i
                        rsv(ieq) = tv(i,j,k)
                     END DO
                  END DO
               END DO
            END IF
         END DO
      END DO
   END DO
ELSE

   ! Fine grid processes receive the coarse grid's correction
   ! Determine where it came from
   coord(1) = zrank - zt
   coord(2) = yrank - yt
   coord(3) = xrank - xt
   CALL MPI_CART_RANK(allcomm,coord,trank,ierr)
   tag = 300 + nrank
   CALL MPI_RECV(rsv,cneq,MPI_DOUBLE_PRECISION,trank,tag,allcomm,istat,ierr)
END IF

DEALLOCATE(tv)

RETURN
END SUBROUTINE corrscat
