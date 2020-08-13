SUBROUTINE corrscatf(v,spx,spy,spz,kc,jc,ic,cneq,rsv,f2)

!------------------------------------------------------------------
!
! Coarse grid processes send correction to fine grid processes.
!
!------------------------------------------------------------------

USE invar
IMPLICIT NONE
INTEGER, INTENT(IN) :: v, spx, spy, spz, kc, jc, ic, cneq
INTEGER :: spx2, spy2, spz2, xt, yt, zt, pk, pj, pi, kt, jt, it, k, j, i, ieq
INTEGER :: kaa, jaa, iaa, ijc, ix1, ix2, trank, tag, ierr
INTEGER :: kcc, jcc, icc, jx1, jx2, jeq
INTEGER, DIMENSION(3) :: coord, istat
REAL*8, DIMENSION(cneq), INTENT(OUT) :: rsv
REAL*8, DIMENSION(cneq) :: tsv
REAL*8, DIMENSION(neq), INTENT(IN) :: f2

INCLUDE 'mpif.h'

! Need zcf/ycf/xcf for communication
kcc = zcf(v)
jcc = ycf(v)
icc = xcf(v)
kaa = kcc - 1
jaa = jcc - 1
iaa = icc - 1
ijc = ic*jc

! Parameters to determine send or receive
spx2 = spx*icc
spy2 = spy*jcc
spz2 = spz*kcc
xt  = MOD(xrank,spx2)
yt  = MOD(yrank,spy2)
zt  = MOD(zrank,spz2)

! If part of the previous stage, need to send
IF (xt == 0 .AND. yt == 0 .AND. zt == 0) THEN
   ! Now copy tv into tsv and send
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
               DO k = kt + 1, kt + sdnz/kc
                  ix1 = (k-kt-1)*xys/ijc
                  jx1 = (k-1)*xys
                  IF (kcc > sdnz) jx1 = ((k-1+kcc/sdnz)*sdnz/kcc-1)*xys
                  DO j = jt + 1, jt + sdny/jc
                     ix2 = ix1 + (j-jt-1)*sdnx/ic
                     jx2 = jx1 + (j-1)*sdnx
                     IF (jcc > sdny) jx2 = jx1 + ((j-1+jcc/sdny)*sdny/jcc-1)*sdnx
                     DO i = it + 1, it + sdnx/ic
                        ieq = ix2 + i-it
                        jeq = jx2 + i
                        IF (icc > sdnx) jeq = jx2 + (i-1+icc/sdnx)*sdnx/icc
                        tsv(ieq) = f2(jeq)
                     END DO
                  END DO
               END DO
               tag = 300 + trank
               CALL MPI_SEND(tsv,cneq,MPI_DOUBLE_PRECISION,trank,tag,allcomm,ierr)
            ELSE
               DO k = 1, sdnz/kc
                  ix1 = (k-1)*xys/ijc
                  jx1 = (k-1)*xys
                  IF (kcc > sdnz) jx1 = ((k-1+kcc/sdnz)*sdnz/kcc-1)*xys
                  DO j = 1, sdny/jc
                     ix2 = ix1 + (j-1)*sdnx/ic
                     jx2 = jx1 + (j-1)*sdnx
                     IF (jcc > sdny) jx2 = jx1 + ((j-1+jcc/sdny)*sdny/jcc-1)*sdnx
                     DO i = 1, sdnx/ic
                        ieq = ix2 + i
                        jeq = jx2 + i
                        IF (icc > sdnx) jeq = jx2 + (i-1+icc/sdnx)*sdnx/icc
                        rsv(ieq) = f2(jeq)
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

RETURN
END SUBROUTINE corrscatf
