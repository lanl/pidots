SUBROUTINE pgsblk(bit,itmx,cnvf,phi,phiold,psii,psio,oldpsi)

!-------------------------------------------------------------
!
!  Parallel Gauss-Seidel  -  Black
!
!  Check the convergence for all sub-domains in the system
!   If the system is not converged pass the black angular
!   flux data back to the red sub-domains. Update angular flux
!   at reflective edges.
!  If convergence is met, say so and exit the PGS loop
!
!-------------------------------------------------------------

USE invar
IMPLICIT NONE
INTEGER, INTENT(IN) :: bit, itmx
INTEGER, INTENT(OUT) :: cnvf
INTEGER :: i, j, k, i1, sdi, sdii, sdio, reset
INTEGER :: myztag, myytag, myxtag, znxt, zprv, ynxt, yprv, xnxt, xprv
INTEGER :: znxtag, zprtag, ynxtag, yprtag, xnxtag, xprtag, neqz, neqy, neqx, ieq, ierr
INTEGER, DIMENSION(3) :: istat
REAL*8 :: df, dfmx
REAL*8, DIMENSION(nsdp) :: dfsd
REAL*8, DIMENSION(neq,nsdp), INTENT(IN) :: phi
REAL*8, DIMENSION(neq,nsdp), INTENT(OUT) :: phiold
REAL*8, DIMENSION(bcs,8,nsdp), INTENT(INOUT) :: psii, oldpsi
REAL*8, DIMENSION(bcs,8,nsdp), INTENT(OUT) :: psio

INCLUDE 'mpif.h'

! Set up some tags
myztag = 100 + zrank
myytag = 100 + yrank
myxtag = 100 + xrank
znxt = zrank + 1
znxtag = 100 + znxt
zprv = zrank - 1
zprtag = 100 + zprv
ynxt = yrank + 1
ynxtag = 100 + ynxt
yprv = yrank - 1
yprtag = 100 + yprv
xnxt = xrank + 1
xnxtag = 100 + xnxt
xprv = xrank - 1
xprtag = 100 + xprv

! Check convergence of scalar flux
! Each process checks its own fluxes, then perform MPI Reduction
IF (MINVAL(ABS(phiold)) >= tolr) THEN
   DO i = 1, nsdp
      dfsd(i) = MAXVAL(ABS((phi(:,i) - phiold(:,i))/phiold(:,i)))
   END DO
ELSE
   DO i = 1, nsdp
      dfsd(i) = MAXVAL(ABS(phi(:,i) - phiold(:,i)))
   END DO
END IF
df = MAXVAL(dfsd)

! Get the overall max difference
CALL MPI_ALLREDUCE(df,dfmx,1,MPI_DOUBLE_PRECISION,MPI_MAX,allcomm,ierr)

! Only print iteration information if requested
IF (irank == root) THEN
   IF (itp == 1) WRITE(8,111) bit, dfmx
END IF
111 FORMAT(2X,'It ', I5,' Dfmx ', ES11.3)

! Set cnvf to 1 if converged
IF (dfmx < err) cnvf = 1

! Set the previous iterate of the flux equal to the current
phiold = phi

! Get a copy of the old psii for computing the residual for red sub-domains
IF (bit == itmx) THEN
   sdi = 0
   DO k = 1, nzsd
      DO j = 1, nysd
         i1 = MOD(j+k,2) + 1
         DO i = i1, nxsd, 2
            sdi = (k-1)*nysd*nxsd + (j-1)*nxsd + i
            oldpsi(:,:,sdi) = psii(:,:,sdi)
         END DO
      END DO
   END DO
END IF

! Send psio to neighbors. Inner blocks have 24 communications (8 octants to
! 3 neigbors each). Then set the incoming data into psii

!--------------------------------------------------------------------------------
! Perform the copies 'up z' from black to red, internally
neqz = apo*xys
DO k = 1, nzsd-1
   DO j = 1, nysd
      i1 = MOD(j+k+1,2) + 1
      DO i = i1, nxsd, 2
         sdio = (k-1)*nysd*nxsd + (j-1)*nxsd + i
         sdii = sdio + nysd*nxsd
         psii(1:neqz,1,sdii) = psio(1:neqz,1,sdio)
         psii(1:neqz,2,sdii) = psio(1:neqz,2,sdio)
         psii(1:neqz,3,sdii) = psio(1:neqz,3,sdio)
         psii(1:neqz,4,sdii) = psio(1:neqz,4,sdio)
      END DO
   END DO
END DO
! Now perform the sends 'up z' from black to red across processors
reset = (nzsd-1)*nysd*nxsd
DO j = 1, nysd
   i1 = MOD(j+1,2) + 1
   DO i = i1, nxsd, 2
      sdii = (j-1)*nxsd + i
      sdio = sdii + reset
      IF (zrank/=npz-1) CALL MPI_SEND(psio(1,1,sdio),neqz,MPI_DOUBLE_PRECISION,znxt,znxtag,zcomm,ierr)
      IF (zrank/=0) CALL MPI_RECV(psii(1,1,sdii),neqz,MPI_DOUBLE_PRECISION,zprv,myztag,zcomm,istat,ierr)
      IF (zrank/=npz-1) CALL MPI_SEND(psio(1,2,sdio),neqz,MPI_DOUBLE_PRECISION,znxt,znxtag,zcomm,ierr)
      IF (zrank/=0) CALL MPI_RECV(psii(1,2,sdii),neqz,MPI_DOUBLE_PRECISION,zprv,myztag,zcomm,istat,ierr)
      IF (zrank/=npz-1) CALL MPI_SEND(psio(1,3,sdio),neqz,MPI_DOUBLE_PRECISION,znxt,znxtag,zcomm,ierr)
      IF (zrank/=0) CALL MPI_RECV(psii(1,3,sdii),neqz,MPI_DOUBLE_PRECISION,zprv,myztag,zcomm,istat,ierr)
      IF (zrank/=npz-1) CALL MPI_SEND(psio(1,4,sdio),neqz,MPI_DOUBLE_PRECISION,znxt,znxtag,zcomm,ierr)
      IF (zrank/=0) CALL MPI_RECV(psii(1,4,sdii),neqz,MPI_DOUBLE_PRECISION,zprv,myztag,zcomm,istat,ierr)
   END DO
END DO
! Perform the copies 'down z' from black to red, internally
DO k = 2, nzsd
   DO j = 1, nysd
      i1 = MOD(j+k+1,2) + 1
      DO i = i1, nxsd, 2
         sdio = (k-1)*nysd*nxsd + (j-1)*nxsd + i
         sdii = sdio - nysd*nxsd
         psii(1:neqz,5,sdii) = psio(1:neqz,5,sdio)
         psii(1:neqz,6,sdii) = psio(1:neqz,6,sdio)
         psii(1:neqz,7,sdii) = psio(1:neqz,7,sdio)
         psii(1:neqz,8,sdii) = psio(1:neqz,8,sdio)
      END DO
   END DO
END DO
! Now perform the sends 'down z' from black to red across processors
DO j = 1, nysd
   i1 = MOD(j,2) + 1
   DO i = i1, nxsd, 2
      sdio = (j-1)*nxsd + i
      sdii = sdio + reset
      IF (zrank/=0) CALL MPI_SEND(psio(1,5,sdio),neqz,MPI_DOUBLE_PRECISION,zprv,zprtag,zcomm,ierr)
      IF (zrank/=npz-1) CALL MPI_RECV(psii(1,5,sdii),neqz,MPI_DOUBLE_PRECISION,znxt,myztag,zcomm,istat,ierr)
      IF (zrank/=0) CALL MPI_SEND(psio(1,6,sdio),neqz,MPI_DOUBLE_PRECISION,zprv,zprtag,zcomm,ierr)
      IF (zrank/=npz-1) CALL MPI_RECV(psii(1,6,sdii),neqz,MPI_DOUBLE_PRECISION,znxt,myztag,zcomm,istat,ierr)
      IF (zrank/=0) CALL MPI_SEND(psio(1,7,sdio),neqz,MPI_DOUBLE_PRECISION,zprv,zprtag,zcomm,ierr)
      IF (zrank/=npz-1) CALL MPI_RECV(psii(1,7,sdii),neqz,MPI_DOUBLE_PRECISION,znxt,myztag,zcomm,istat,ierr)
      IF (zrank/=0) CALL MPI_SEND(psio(1,8,sdio),neqz,MPI_DOUBLE_PRECISION,zprv,zprtag,zcomm,ierr)
      IF (zrank/=npz-1) CALL MPI_RECV(psii(1,8,sdii),neqz,MPI_DOUBLE_PRECISION,znxt,myztag,zcomm,istat,ierr)
   END DO
END DO

!----------------------------------------------------------------------------------------------------
! Perform the copies 'up y' from black to red, internally
ieq = neqz + 1
neqy = apo*xzs
DO k = 1, nzsd
   DO j = 1, nysd-1
      i1 = MOD(j+k+1,2) + 1
      DO i = i1, nxsd, 2
         sdio = (k-1)*nysd*nxsd + (j-1)*nxsd + i
         sdii = sdio + nxsd
         psii(ieq:neqz+neqy,1,sdii) = psio(ieq:neqz+neqy,1,sdio)
         psii(ieq:neqz+neqy,2,sdii) = psio(ieq:neqz+neqy,2,sdio)
         psii(ieq:neqz+neqy,5,sdii) = psio(ieq:neqz+neqy,5,sdio)
         psii(ieq:neqz+neqy,6,sdii) = psio(ieq:neqz+neqy,6,sdio)
      END DO
   END DO
END DO
! Now perform the sends 'up y' from black to red across processors
reset = (nysd-1)*nxsd
DO k = 1, nzsd
   i1 = MOD(k+1,2) + 1
   DO i = i1, nxsd, 2
      sdii = (k-1)*nysd*nxsd + i
      sdio = sdii + reset
      IF (yrank/=npy-1) CALL MPI_SEND(psio(ieq,1,sdio),neqy,MPI_DOUBLE_PRECISION,ynxt,ynxtag,ycomm,ierr)
      IF (yrank/=0) CALL MPI_RECV(psii(ieq,1,sdii),neqy,MPI_DOUBLE_PRECISION,yprv,myytag,ycomm,istat,ierr)
      IF (yrank/=npy-1) CALL MPI_SEND(psio(ieq,2,sdio),neqy,MPI_DOUBLE_PRECISION,ynxt,ynxtag,ycomm,ierr)
      IF (yrank/=0) CALL MPI_RECV(psii(ieq,2,sdii),neqy,MPI_DOUBLE_PRECISION,yprv,myytag,ycomm,istat,ierr)
      IF (yrank/=npy-1) CALL MPI_SEND(psio(ieq,5,sdio),neqy,MPI_DOUBLE_PRECISION,ynxt,ynxtag,ycomm,ierr)
      IF (yrank/=0) CALL MPI_RECV(psii(ieq,5,sdii),neqy,MPI_DOUBLE_PRECISION,yprv,myytag,ycomm,istat,ierr)
      IF (yrank/=npy-1) CALL MPI_SEND(psio(ieq,6,sdio),neqy,MPI_DOUBLE_PRECISION,ynxt,ynxtag,ycomm,ierr)
      IF (yrank/=0) CALL MPI_RECV(psii(ieq,6,sdii),neqy,MPI_DOUBLE_PRECISION,yprv,myytag,ycomm,istat,ierr)
   END DO
END DO
! Perform the copies 'down y' from black to red, internally
DO k = 1, nzsd
   DO j = 2, nysd
      i1 = MOD(j+k+1,2) + 1
      DO i = i1, nxsd, 2
         sdio = (k-1)*nysd*nxsd + (j-1)*nxsd + i
         sdii = sdio - nxsd
         psii(ieq:neqz+neqy,3,sdii) = psio(ieq:neqz+neqy,3,sdio)
         psii(ieq:neqz+neqy,4,sdii) = psio(ieq:neqz+neqy,4,sdio)
         psii(ieq:neqz+neqy,7,sdii) = psio(ieq:neqz+neqy,7,sdio)
         psii(ieq:neqz+neqy,8,sdii) = psio(ieq:neqz+neqy,8,sdio)
      END DO
   END DO
END DO
! Now perform the sends 'down y' from black to red across processors
reset = (nysd-1)*nxsd
DO k = 1, nzsd
   i1 = MOD(k,2) + 1
   DO i = i1, nxsd, 2
      sdio = (k-1)*nysd*nxsd + i
      sdii = sdio + reset
      IF (yrank/=0) CALL MPI_SEND(psio(ieq,3,sdio),neqy,MPI_DOUBLE_PRECISION,yprv,yprtag,ycomm,ierr)
      IF (yrank/=npy-1) CALL MPI_RECV(psii(ieq,3,sdii),neqy,MPI_DOUBLE_PRECISION,ynxt,myytag,ycomm,istat,ierr)
      IF (yrank/=0) CALL MPI_SEND(psio(ieq,4,sdio),neqy,MPI_DOUBLE_PRECISION,yprv,yprtag,ycomm,ierr)
      IF (yrank/=npy-1) CALL MPI_RECV(psii(ieq,4,sdii),neqy,MPI_DOUBLE_PRECISION,ynxt,myytag,ycomm,istat,ierr)
      IF (yrank/=0) CALL MPI_SEND(psio(ieq,7,sdio),neqy,MPI_DOUBLE_PRECISION,yprv,yprtag,ycomm,ierr)
      IF (yrank/=npy-1) CALL MPI_RECV(psii(ieq,7,sdii),neqy,MPI_DOUBLE_PRECISION,ynxt,myytag,ycomm,istat,ierr)
      IF (yrank/=0) CALL MPI_SEND(psio(ieq,8,sdio),neqy,MPI_DOUBLE_PRECISION,yprv,yprtag,ycomm,ierr)
      IF (yrank/=npy-1) CALL MPI_RECV(psii(ieq,8,sdii),neqy,MPI_DOUBLE_PRECISION,ynxt,myytag,ycomm,istat,ierr)
   END DO
END DO

!--------------------------------------------------------------------------------------------
! Perform the copies 'up x' from black to red, internally
ieq = ieq + neqy
neqx = apo*yzs
DO k = 1, nzsd
   DO j = 1, nysd
      i1 = MOD(j+k+1,2) + 1
      DO i = i1, nxsd-1, 2
         sdio = (k-1)*nysd*nxsd + (j-1)*nxsd + i
         sdii = sdio + 1
         psii(ieq:bcs,1,sdii) = psio(ieq:bcs,1,sdio)
         psii(ieq:bcs,4,sdii) = psio(ieq:bcs,4,sdio)
         psii(ieq:bcs,5,sdii) = psio(ieq:bcs,5,sdio)
         psii(ieq:bcs,8,sdii) = psio(ieq:bcs,8,sdio)
      END DO
   END DO
END DO
! Now perform the sends 'up x' from black to red across processors
reset = nxsd-1
DO k = 1, nzsd
   i1 = MOD(k+1,2) + 1
   DO j = i1, nysd, 2
      sdio = (k-1)*nysd*nxsd + j*nxsd
      sdii = sdio - reset
      IF (xrank/=npx-1) CALL MPI_SEND(psio(ieq,1,sdio),neqx,MPI_DOUBLE_PRECISION,xnxt,xnxtag,xcomm,ierr)
      IF (xrank/=0) CALL MPI_RECV(psii(ieq,1,sdii),neqx,MPI_DOUBLE_PRECISION,xprv,myxtag,xcomm,istat,ierr)
      IF (xrank/=npx-1) CALL MPI_SEND(psio(ieq,4,sdio),neqx,MPI_DOUBLE_PRECISION,xnxt,xnxtag,xcomm,ierr)
      IF (xrank/=0) CALL MPI_RECV(psii(ieq,4,sdii),neqx,MPI_DOUBLE_PRECISION,xprv,myxtag,xcomm,istat,ierr)
      IF (xrank/=npx-1) CALL MPI_SEND(psio(ieq,5,sdio),neqx,MPI_DOUBLE_PRECISION,xnxt,xnxtag,xcomm,ierr)
      IF (xrank/=0) CALL MPI_RECV(psii(ieq,5,sdii),neqx,MPI_DOUBLE_PRECISION,xprv,myxtag,xcomm,istat,ierr)
      IF (xrank/=npx-1) CALL MPI_SEND(psio(ieq,8,sdio),neqx,MPI_DOUBLE_PRECISION,xnxt,xnxtag,xcomm,ierr)
      IF (xrank/=0) CALL MPI_RECV(psii(ieq,8,sdii),neqx,MPI_DOUBLE_PRECISION,xprv,myxtag,xcomm,istat,ierr)
   END DO
END DO
! Perform the copies 'down x' from black to red, internally
DO k = 1, nzsd
   DO j = 1, nysd
      i1 = MOD(j+k,2) + 2
      DO i = i1, nxsd, 2
         sdio = (k-1)*nysd*nxsd + (j-1)*nxsd + i
         sdii = sdio - 1
         psii(ieq:bcs,2,sdii) = psio(ieq:bcs,2,sdio)
         psii(ieq:bcs,3,sdii) = psio(ieq:bcs,3,sdio)
         psii(ieq:bcs,6,sdii) = psio(ieq:bcs,6,sdio)
         psii(ieq:bcs,7,sdii) = psio(ieq:bcs,7,sdio)
      END DO
   END DO
END DO
! Now perform the sends 'down x' from black to red across processors
DO k = 1, nzsd
   i1 = MOD(k,2) + 1
   DO j = i1, nysd, 2
      sdii = (k-1)*nysd*nxsd + j*nxsd
      sdio = sdii - reset
      IF (xrank/=0) CALL MPI_SEND(psio(ieq,2,sdio),neqx,MPI_DOUBLE_PRECISION,xprv,xprtag,xcomm,ierr)
      IF (xrank/=npx-1) CALL MPI_RECV(psii(ieq,2,sdii),neqx,MPI_DOUBLE_PRECISION,xnxt,myxtag,xcomm,istat,ierr)
      IF (xrank/=0) CALL MPI_SEND(psio(ieq,3,sdio),neqx,MPI_DOUBLE_PRECISION,xprv,xprtag,xcomm,ierr)
      IF (xrank/=npx-1) CALL MPI_RECV(psii(ieq,3,sdii),neqx,MPI_DOUBLE_PRECISION,xnxt,myxtag,xcomm,istat,ierr)
      IF (xrank/=0) CALL MPI_SEND(psio(ieq,6,sdio),neqx,MPI_DOUBLE_PRECISION,xprv,xprtag,xcomm,ierr)
      IF (xrank/=npx-1) CALL MPI_RECV(psii(ieq,6,sdii),neqx,MPI_DOUBLE_PRECISION,xnxt,myxtag,xcomm,istat,ierr)
      IF (xrank/=0) CALL MPI_SEND(psio(ieq,7,sdio),neqx,MPI_DOUBLE_PRECISION,xprv,xprtag,xcomm,ierr)
      IF (xrank/=npx-1) CALL MPI_RECV(psii(ieq,7,sdii),neqx,MPI_DOUBLE_PRECISION,xnxt,myxtag,xcomm,istat,ierr)
   END DO
END DO

RETURN
END SUBROUTINE pgsblk
