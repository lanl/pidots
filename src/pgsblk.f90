SUBROUTINE pgsblk(bit,tsolve)

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
USE solvar
IMPLICIT NONE
INTEGER, INTENT(IN) :: bit
INTEGER :: i, j, k, i1, sdii, sdio, reset, tmpi
integer, dimension(1) :: ti
INTEGER :: myztag, myytag, myxtag, znxt, zprv, ynxt, yprv, xnxt, xprv
INTEGER :: znxtag, zprtag, ynxtag, yprtag, xnxtag, xprtag, neqz, neqy, neqx, ieq, ierr
INTEGER, DIMENSION(3) :: istat
REAL*8 :: tmpr
REAL*8, DIMENSION(nsdp) :: dfsd, tmpr1, tmpr2
real*8, dimension(2) :: df, dfmx, tr1, tr2
REAL*8, INTENT(OUT) :: tsolve

REAL*8, DIMENSION(:,:,:), ALLOCATABLE:: transfer_array

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
      ti = MAXLOC(ABS((phi(:,i) - phiold(:,i))/phiold(:,i)))
      tmpr1(i) = ABS(phi(ti(1),i) - phiold(ti(1),i))
      tmpr2(i) = ABS(phiold(ti(1),i))
   END DO
ELSE
   DO i = 1, nsdp
      dfsd(i) = MAXVAL(ABS(phi(:,i) - phiold(:,i)))
   END DO
END IF
df(1) = MAXVAL(dfsd)
df(2) = nrank
ti = MAXLOC(dfsd)
tr1(1) = tmpr1(ti(1))
tr1(2) = tmpr2(ti(1))

! Get the overall max difference
CALL MPI_ALLREDUCE(df,dfmx,1,MPI_2DOUBLE_PRECISION,MPI_MAXloc,allcomm,ierr)

tmpr = dfmx(1)
tmpi = dfmx(2)

! determine process with greatest, then get its flux, delta-flux
IF (nrank == tmpi .and. tmpi /= 0) THEN
   CALL MPI_SEND(tr1,2,MPI_DOUBLE_PRECISION,0,1,allcomm,ierr)
END IF
IF (nrank == 0 .and. tmpi/=0) THEN
   CALL MPI_RECV(tr2,2,MPI_DOUBLE_PRECISION,tmpi,1,allcomm,istat,ierr)
END IF
if (tmpi == 0) then
   tr2 = tr1
end if

IF (tmpr > err .AND. bit < bitmx) THEN
   ! Only print iteration information if requested
   IF (irank == root) THEN
      IF (itp == 1) WRITE(8,111) bit, tmpr, tr2(1), tr2(2)
   END IF
   ! Set the previous iterate of the flux equal to the current
   phiold = phi

   ! Send psio to neighbors. Inner blocks have 24 communications (8 octants to
   ! 3 neigbors each). Then set the incoming data into psii

   !!------------------------------------------------------------------------------
   !!------------------------------------------------------------------------------
   !!UP-Z--------------------------------------------------------------------------
   !!------------------------------------------------------------------------------
   !!------------------------------------------------------------------------------

   !!Equation count is angles_per_octant * number of cells per face of processor subdomain
   neqz = apo*xys

   !!Alocate transfer array to be neqz by number of face subdomains
   ALLOCATE(transfer_array(neqz, nysd, nxsd))
   transfer_array=0

   !!------------------------------------------------------------------------------
   !!UP-Z, INTERNAL----------------------------------------------------------------
   !!------------------------------------------------------------------------------
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

   !!------------------------------------------------------------------------------
   !!UP-Z, OCTANT 1----------------------------------------------------------------
   !!------------------------------------------------------------------------------
   reset = (nzsd-1)*nysd*nxsd
   DO j = 1, nysd
      i1 = MOD(j+1,2) + 1
      DO i = i1, nxsd, 2
         sdii = (j-1)*nxsd + i
         sdio = sdii + reset
         IF (zrank/=npz-1) THEN
           transfer_array(1:neqz,j,i) = psio(1:neqz,1,sdio)
           !!CALL MPI_SEND(psio(1:neqz,1,sdio),neqz,MPI_DOUBLE_PRECISION,znxt,znxtag,zcomm,ierr)
         END IF
       !!IF (zrank/=0) CALL MPI_RECV(psii(1:neqz,1,sdii),neqz,MPI_DOUBLE_PRECISION,zprv,myztag,zcomm,istat,ierr)
     END DO
   END DO

   !!Send all subdomain data
   IF (zrank/=npz-1) CALL MPI_SEND(transfer_array,neqz*nysd*nxsd,MPI_DOUBLE_PRECISION,znxt,znxtag,zcomm,ierr)
   !!Recieve all subdomain data
   IF (zrank/=0)     CALL MPI_RECV(transfer_array,neqz*nysd*nxsd,MPI_DOUBLE_PRECISION,zprv,myztag,zcomm,istat,ierr)

   DO j = 1, nysd
      i1 = MOD(j+1,2) + 1
      DO i = i1, nxsd, 2
        sdii = (j-1)*nxsd + i
        sdio = sdii + reset
        IF (zrank/=0) psii(1:neqz,1,sdii) = transfer_array(1:neqz,j,i)
      END DO
    END DO
    !!------------------------------------------------------------------------------
    !!UP-Z, OCTANT 2----------------------------------------------------------------
    !!------------------------------------------------------------------------------

    DO j = 1, nysd
       i1 = MOD(j+1,2) + 1
       DO i = i1, nxsd, 2
          sdii = (j-1)*nxsd + i
          sdio = sdii + reset
          ! Octant 2
          IF (zrank/=npz-1) THEN
            transfer_array(1:neqz,j,i) = psio(1:neqz,2,sdio)
            !!CALL MPI_SEND(psio(1,2,sdio),neqz,MPI_DOUBLE_PRECISION,znxt,znxtag,zcomm,ierr)
          END IF
        !!IF (zrank/=0) CALL MPI_RECV(psii(1,2,sdii),neqz,MPI_DOUBLE_PRECISION,zprv,myztag,zcomm,istat,ierr)
      END DO
    END DO

    !!Send all subdomain data
    IF (zrank/=npz-1) CALL MPI_SEND(transfer_array,neqz*nysd*nxsd,MPI_DOUBLE_PRECISION,znxt,znxtag,zcomm,ierr)
    !!Recieve all subdomain data
    IF (zrank/=0)     CALL MPI_RECV(transfer_array,neqz*nysd*nxsd,MPI_DOUBLE_PRECISION,zprv,myztag,zcomm,istat,ierr)

    DO j = 1, nysd
       i1 = MOD(j+1,2) + 1
       DO i = i1, nxsd, 2
          sdii = (j-1)*nxsd + i
          sdio = sdii + reset
          ! Octant 2
          IF (zrank .NE. 0) psii(1:neqz,2,sdii) = transfer_array(1:neqz,j,i)
          !!IF (zrank/=0) CALL MPI_RECV(psii(1,2,sdii),neqz,MPI_DOUBLE_PRECISION,zprv,myztag,zcomm,istat,ierr)
      END DO
    END DO


    !!------------------------------------------------------------------------------
    !!UP-Z, OCTANT 3----------------------------------------------------------------
    !!------------------------------------------------------------------------------

    DO j = 1, nysd
       i1 = MOD(j+1,2) + 1
       DO i = i1, nxsd, 2
          sdii = (j-1)*nxsd + i
          sdio = sdii + reset
          ! Octant 3
          IF (zrank/=npz-1) THEN
            transfer_array(1:neqz,j,i) = psio(1:neqz,3,sdio)
            !!CALL MPI_SEND(psio(1,3,sdio),neqz,MPI_DOUBLE_PRECISION,znxt,znxtag,zcomm,ierr)
          END IF
        !!IF (zrank/=0) CALL MPI_RECV(psii(1,3,sdii),neqz,MPI_DOUBLE_PRECISION,zprv,myztag,zcomm,istat,ierr)
      END DO
    END DO

    !!Send all subdomain data
    IF (zrank/=npz-1) CALL MPI_SEND(transfer_array,neqz*nysd*nxsd,MPI_DOUBLE_PRECISION,znxt,znxtag,zcomm,ierr)
    !!Recieve all subdomain data
    IF (zrank/=0)     CALL MPI_RECV(transfer_array,neqz*nysd*nxsd,MPI_DOUBLE_PRECISION,zprv,myztag,zcomm,istat,ierr)

    DO j = 1, nysd
       i1 = MOD(j+1,2) + 1
       DO i = i1, nxsd, 2
          sdii = (j-1)*nxsd + i
          sdio = sdii + reset
          IF (zrank .NE. 0) psii(1:neqz,3,sdii) = transfer_array(1:neqz,j,i)
      END DO
    END DO

    !!------------------------------------------------------------------------------
    !!UP-Z, OCTANT 4----------------------------------------------------------------
    !!------------------------------------------------------------------------------

    DO j = 1, nysd
       i1 = MOD(j+1,2) + 1
       DO i = i1, nxsd, 2
          sdii = (j-1)*nxsd + i
          sdio = sdii + reset
          IF (zrank/=npz-1) THEN
            transfer_array(1:neqz,j,i) = psio(1:neqz,4,sdio)
            !!CALL MPI_SEND(psio(1,3,sdio),neqz,MPI_DOUBLE_PRECISION,znxt,znxtag,zcomm,ierr)
          END IF
        !!IF (zrank/=0) CALL MPI_RECV(psii(1,3,sdii),neqz,MPI_DOUBLE_PRECISION,zprv,myztag,zcomm,istat,ierr)
      END DO
    END DO

    !!Send all subdomain data
    IF (zrank/=npz-1) CALL MPI_SEND(transfer_array,neqz*nysd*nxsd,MPI_DOUBLE_PRECISION,znxt,znxtag,zcomm,ierr)
    !!Recieve all subdomain data
    IF (zrank/=0)     CALL MPI_RECV(transfer_array,neqz*nysd*nxsd,MPI_DOUBLE_PRECISION,zprv,myztag,zcomm,istat,ierr)

    DO j = 1, nysd
       i1 = MOD(j+1,2) + 1
       DO i = i1, nxsd, 2
          sdii = (j-1)*nxsd + i
          sdio = sdii + reset
          IF (zrank .NE. 0) psii(1:neqz,4,sdii) = transfer_array(1:neqz,j,i)
      END DO
    END DO

  !!------------------------------------------------------------------------------
  !!------------------------------------------------------------------------------
  !!DN-Z--------------------------------------------------------------------------
  !!------------------------------------------------------------------------------
  !!------------------------------------------------------------------------------

  transfer_array = 0

  !!------------------------------------------------------------------------------
  !!DN-Z, INTERNAL----------------------------------------------------------------
  !!------------------------------------------------------------------------------

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
   !!------------------------------------------------------------------------------
   !!DN-Z, OCTANT 5----------------------------------------------------------------
   !!------------------------------------------------------------------------------

   DO j = 1, nysd
      i1 = MOD(j,2) + 1
      DO i = i1, nxsd, 2
         sdio = (j-1)*nxsd + i
         sdii = sdio + reset
         IF (zrank/=0) THEN
           transfer_array(1:neqz,j,i) = psio(1:neqz,5,sdio)
           !!CALL MPI_SEND(psio(1,5,sdio),neqz,MPI_DOUBLE_PRECISION,zprv,zprtag,zcomm,ierr)
         END IF
       !!IF (zrank/=npz-1) CALL MPI_RECV(psii(1,5,sdii),neqz,MPI_DOUBLE_PRECISION,znxt,myztag,zcomm,istat,ierr)
     END DO
   END DO

   !!Send all subdomain data
   IF (zrank .NE. 0)     CALL MPI_SEND(transfer_array,neqz*nysd*nxsd,MPI_DOUBLE_PRECISION,zprv,zprtag,zcomm,ierr)
   !!Recieve all subdomain data
   IF (zrank .NE. npz-1) CALL MPI_RECV(transfer_array,neqz*nysd*nxsd,MPI_DOUBLE_PRECISION,znxt,myztag,zcomm,istat,ierr)

   DO j = 1, nysd
      i1 = MOD(j,2) + 1
      DO i = i1, nxsd, 2
         sdio = (j-1)*nxsd + i
         sdii = sdio + reset
         IF(zrank .NE. npz-1) psii(1:neqz,5,sdii) = transfer_array(1:neqz, j, i)
     END DO
   END DO

   !!------------------------------------------------------------------------------
   !!DN-Z, OCTANT 6----------------------------------------------------------------
   !!------------------------------------------------------------------------------

   DO j = 1, nysd
      i1 = MOD(j,2) + 1
      DO i = i1, nxsd, 2
         sdio = (j-1)*nxsd + i
         sdii = sdio + reset
         IF (zrank/=0) THEN
           transfer_array(1:neqz,j,i) = psio(1:neqz,6,sdio)
           !!CALL MPI_SEND(psio(1,5,sdio),neqz,MPI_DOUBLE_PRECISION,zprv,zprtag,zcomm,ierr)
         END IF
       !!IF (zrank/=npz-1) CALL MPI_RECV(psii(1,5,sdii),neqz,MPI_DOUBLE_PRECISION,znxt,myztag,zcomm,istat,ierr)
     END DO
   END DO

   !!Send all subdomain data
   IF (zrank .NE. 0)     CALL MPI_SEND(transfer_array,neqz*nysd*nxsd,MPI_DOUBLE_PRECISION,zprv,zprtag,zcomm,ierr)
   !!Recieve all subdomain data
   IF (zrank .NE. npz-1) CALL MPI_RECV(transfer_array,neqz*nysd*nxsd,MPI_DOUBLE_PRECISION,znxt,myztag,zcomm,istat,ierr)

   DO j = 1, nysd
      i1 = MOD(j,2) + 1
      DO i = i1, nxsd, 2
         sdio = (j-1)*nxsd + i
         sdii = sdio + reset
         IF(zrank .NE. npz-1) psii(1:neqz,6,sdii) = transfer_array(1:neqz, j, i)
     END DO
   END DO

   !!------------------------------------------------------------------------------
   !!DN-Z, OCTANT 7----------------------------------------------------------------
   !!------------------------------------------------------------------------------

   DO j = 1, nysd
      i1 = MOD(j,2) + 1
      DO i = i1, nxsd, 2
         sdio = (j-1)*nxsd + i
         sdii = sdio + reset
         IF (zrank/=0) THEN
           transfer_array(1:neqz,j,i) = psio(1:neqz,7,sdio)
           !!CALL MPI_SEND(psio(1,5,sdio),neqz,MPI_DOUBLE_PRECISION,zprv,zprtag,zcomm,ierr)
         END IF
       !!IF (zrank/=npz-1) CALL MPI_RECV(psii(1,5,sdii),neqz,MPI_DOUBLE_PRECISION,znxt,myztag,zcomm,istat,ierr)
     END DO
   END DO

   !!Send all subdomain data
   IF (zrank .NE. 0)     CALL MPI_SEND(transfer_array,neqz*nysd*nxsd,MPI_DOUBLE_PRECISION,zprv,zprtag,zcomm,ierr)
   !!Recieve all subdomain data
   IF (zrank .NE. npz-1) CALL MPI_RECV(transfer_array,neqz*nysd*nxsd,MPI_DOUBLE_PRECISION,znxt,myztag,zcomm,istat,ierr)

   DO j = 1, nysd
      i1 = MOD(j,2) + 1
      DO i = i1, nxsd, 2
         sdio = (j-1)*nxsd + i
         sdii = sdio + reset
         IF(zrank .NE. npz-1) psii(1:neqz,7,sdii) = transfer_array(1:neqz, j, i)
     END DO
   END DO

   !!------------------------------------------------------------------------------
   !!DN-Z, OCTANT 8----------------------------------------------------------------
   !!------------------------------------------------------------------------------

   DO j = 1, nysd
      i1 = MOD(j,2) + 1
      DO i = i1, nxsd, 2
         sdio = (j-1)*nxsd + i
         sdii = sdio + reset
         IF (zrank/=0) THEN
           transfer_array(1:neqz,j,i) = psio(1:neqz,8,sdio)
           !!CALL MPI_SEND(psio(1,5,sdio),neqz,MPI_DOUBLE_PRECISION,zprv,zprtag,zcomm,ierr)
         END IF
       !!IF (zrank/=npz-1) CALL MPI_RECV(psii(1,5,sdii),neqz,MPI_DOUBLE_PRECISION,znxt,myztag,zcomm,istat,ierr)
     END DO
   END DO

   !!Send all subdomain data
   IF (zrank .NE. 0)     CALL MPI_SEND(transfer_array,neqz*nysd*nxsd,MPI_DOUBLE_PRECISION,zprv,zprtag,zcomm,ierr)
   !!Recieve all subdomain data
   IF (zrank .NE. npz-1) CALL MPI_RECV(transfer_array,neqz*nysd*nxsd,MPI_DOUBLE_PRECISION,znxt,myztag,zcomm,istat,ierr)

   DO j = 1, nysd
      i1 = MOD(j,2) + 1
      DO i = i1, nxsd, 2
         sdio = (j-1)*nxsd + i
         sdii = sdio + reset
         IF(zrank .NE. npz-1) psii(1:neqz,8,sdii) = transfer_array(1:neqz, j, i)
     END DO
   END DO

   DEALLOCATE(transfer_array)

   ! Reflective boundary conditions
   IF (zrank == 0 .AND. bc(1) == 1) THEN               ! Reflective at low z
      DO i = 1, nxsd*nysd                              ! All z = 1 sub-domains (red and black)
         psii(1:neqz,1,i) = psio(1:neqz,5,i)
         psii(1:neqz,2,i) = psio(1:neqz,6,i)
         psii(1:neqz,3,i) = psio(1:neqz,7,i)
         psii(1:neqz,4,i) = psio(1:neqz,8,i)
      END DO
   END IF
   IF (zrank == npz-1 .AND. bc(2) == 1) THEN           ! Reflective at high z
      DO i = (nzsd-1)*nxsd*nysd + 1, nsdp                 ! All z = nzsd sub-domains
         psii(1:neqz,5,i) = psio(1:neqz,1,i)
         psii(1:neqz,6,i) = psio(1:neqz,2,i)
         psii(1:neqz,7,i) = psio(1:neqz,3,i)
         psii(1:neqz,8,i) = psio(1:neqz,4,i)
      END DO
   END IF

   !!------------------------------------------------------------------------------
   !!------------------------------------------------------------------------------
   !!UP-Y--------------------------------------------------------------------------
   !!------------------------------------------------------------------------------
   !!------------------------------------------------------------------------------
   !!Equation count is angles_per_octant * number of cells per face of processor subdomain
   ieq = neqz + 1
   neqy = apo*xzs

   !!Alocate transfer array to be neqz by number of face subdomains
   ALLOCATE(transfer_array(neqy, nzsd, nxsd))
   transfer_array=0

   !!------------------------------------------------------------------------------
   !!UP-Y, INTERNAL----------------------------------------------------------------
   !!------------------------------------------------------------------------------
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
   !!------------------------------------------------------------------------------
   !!UP-Y, OCTANT 1----------------------------------------------------------------
   !!------------------------------------------------------------------------------
   reset = (nysd-1)*nxsd
   DO k = 1, nzsd
      i1 = MOD(k+1,2) + 1
      DO i = i1, nxsd, 2
         sdii = (k-1)*nysd*nxsd + i
         sdio = sdii + reset
         IF (yrank/=npy-1) THEN
           transfer_array(1:neqy, k, i) = psio(ieq:ieq+neqy-1,1,sdio)
           !!CALL MPI_SEND(psio(ieq,1,sdio),neqy,MPI_DOUBLE_PRECISION,ynxt,ynxtag,ycomm,ierr)
         END IF
         !!IF (yrank/=0) CALL MPI_RECV(psii(ieq,1,sdii),neqy,MPI_DOUBLE_PRECISION,yprv,myytag,ycomm,istat,ierr)
       END DO
    END DO

   IF (yrank/=npy-1) CALL MPI_SEND(transfer_array,neqy*nzsd*nxsd,MPI_DOUBLE_PRECISION,ynxt,ynxtag,ycomm,ierr)
   IF (yrank/=0)     CALL MPI_RECV(transfer_array,neqy*nzsd*nxsd,MPI_DOUBLE_PRECISION,yprv,myytag,ycomm,istat,ierr)

    DO k = 1, nzsd
       i1 = MOD(k+1,2) + 1
       DO i = i1, nxsd, 2
          sdii = (k-1)*nysd*nxsd + i
          sdio = sdii + reset
          IF (yrank/=0) psii(ieq:ieq+neqy-1,1,sdii) = transfer_array(1:neqy, k, i)
          !!IF (yrank/=0) CALL MPI_RECV(psii(ieq,1,sdii),neqy,MPI_DOUBLE_PRECISION,yprv,myytag,ycomm,istat,ierr)
        END DO
     END DO
     !!------------------------------------------------------------------------------
     !!UP-Y, OCTANT 2----------------------------------------------------------------
     !!------------------------------------------------------------------------------
     DO k = 1, nzsd
        i1 = MOD(k+1,2) + 1
        DO i = i1, nxsd, 2
           sdii = (k-1)*nysd*nxsd + i
           sdio = sdii + reset
           IF (yrank/=npy-1) THEN
             transfer_array(1:neqy, k, i) = psio(ieq:ieq+neqy-1,2,sdio)
             !!CALL MPI_SEND(psio(ieq,1,sdio),neqy,MPI_DOUBLE_PRECISION,ynxt,ynxtag,ycomm,ierr)
           END IF
           !!IF (yrank/=0) CALL MPI_RECV(psii(ieq,1,sdii),neqy,MPI_DOUBLE_PRECISION,yprv,myytag,ycomm,istat,ierr)
         END DO
      END DO

     IF (yrank/=npy-1) CALL MPI_SEND(transfer_array,neqy*nzsd*nxsd,MPI_DOUBLE_PRECISION,ynxt,ynxtag,ycomm,ierr)
     IF (yrank/=0)     CALL MPI_RECV(transfer_array,neqy*nzsd*nxsd,MPI_DOUBLE_PRECISION,yprv,myytag,ycomm,istat,ierr)

      DO k = 1, nzsd
         i1 = MOD(k+1,2) + 1
         DO i = i1, nxsd, 2
            sdii = (k-1)*nysd*nxsd + i
            sdio = sdii + reset
            IF (yrank/=0) psii(ieq:ieq+neqy-1,2,sdii) = transfer_array(1:neqy, k, i)
            !!IF (yrank/=0) CALL MPI_RECV(psii(ieq,1,sdii),neqy,MPI_DOUBLE_PRECISION,yprv,myytag,ycomm,istat,ierr)
          END DO
       END DO
     !!------------------------------------------------------------------------------
     !!UP-Y, OCTANT 5----------------------------------------------------------------
     !!------------------------------------------------------------------------------
     DO k = 1, nzsd
        i1 = MOD(k+1,2) + 1
        DO i = i1, nxsd, 2
           sdii = (k-1)*nysd*nxsd + i
           sdio = sdii + reset
           IF (yrank/=npy-1) THEN
             transfer_array(1:neqy, k, i) = psio(ieq:ieq+neqy-1,5,sdio)
             !!CALL MPI_SEND(psio(ieq,1,sdio),neqy,MPI_DOUBLE_PRECISION,ynxt,ynxtag,ycomm,ierr)
           END IF
           !!IF (yrank/=0) CALL MPI_RECV(psii(ieq,1,sdii),neqy,MPI_DOUBLE_PRECISION,yprv,myytag,ycomm,istat,ierr)
         END DO
      END DO

     IF (yrank/=npy-1) CALL MPI_SEND(transfer_array,neqy*nzsd*nxsd,MPI_DOUBLE_PRECISION,ynxt,ynxtag,ycomm,ierr)
     IF (yrank/=0)     CALL MPI_RECV(transfer_array,neqy*nzsd*nxsd,MPI_DOUBLE_PRECISION,yprv,myytag,ycomm,istat,ierr)

      DO k = 1, nzsd
         i1 = MOD(k+1,2) + 1
         DO i = i1, nxsd, 2
            sdii = (k-1)*nysd*nxsd + i
            sdio = sdii + reset
            IF (yrank/=0) psii(ieq:ieq+neqy-1,5,sdii) = transfer_array(1:neqy, k, i)
            !!IF (yrank/=0) CALL MPI_RECV(psii(ieq,1,sdii),neqy,MPI_DOUBLE_PRECISION,yprv,myytag,ycomm,istat,ierr)
          END DO
       END DO
     !!------------------------------------------------------------------------------
     !!UP-Y, OCTANT 6----------------------------------------------------------------
     !!------------------------------------------------------------------------------
     DO k = 1, nzsd
        i1 = MOD(k+1,2) + 1
        DO i = i1, nxsd, 2
           sdii = (k-1)*nysd*nxsd + i
           sdio = sdii + reset
           IF (yrank/=npy-1) THEN
             transfer_array(1:neqy, k, i) = psio(ieq:ieq+neqy-1,6,sdio)
             !!CALL MPI_SEND(psio(ieq,1,sdio),neqy,MPI_DOUBLE_PRECISION,ynxt,ynxtag,ycomm,ierr)
           END IF
           !!IF (yrank/=0) CALL MPI_RECV(psii(ieq,1,sdii),neqy,MPI_DOUBLE_PRECISION,yprv,myytag,ycomm,istat,ierr)
         END DO
      END DO

     IF (yrank/=npy-1) CALL MPI_SEND(transfer_array,neqy*nzsd*nxsd,MPI_DOUBLE_PRECISION,ynxt,ynxtag,ycomm,ierr)
     IF (yrank/=0)     CALL MPI_RECV(transfer_array,neqy*nzsd*nxsd,MPI_DOUBLE_PRECISION,yprv,myytag,ycomm,istat,ierr)

      DO k = 1, nzsd
         i1 = MOD(k+1,2) + 1
         DO i = i1, nxsd, 2
            sdii = (k-1)*nysd*nxsd + i
            sdio = sdii + reset
            IF (yrank/=0) psii(ieq:ieq+neqy-1,6,sdii) = transfer_array(1:neqy, k, i)
            !!IF (yrank/=0) CALL MPI_RECV(psii(ieq,1,sdii),neqy,MPI_DOUBLE_PRECISION,yprv,myytag,ycomm,istat,ierr)
          END DO
       END DO
   !!------------------------------------------------------------------------------
   !!------------------------------------------------------------------------------
   !!DN-Y--------------------------------------------------------------------------
   !!------------------------------------------------------------------------------
   !!------------------------------------------------------------------------------

   !!------------------------------------------------------------------------------
   !!DN-Y, INTERNAL----------------------------------------------------------------
   !!------------------------------------------------------------------------------
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
   !!------------------------------------------------------------------------------
   !!DN-Y, OCTANT 3----------------------------------------------------------------
   !!------------------------------------------------------------------------------
   reset = (nysd-1)*nxsd
   DO k = 1, nzsd
      i1 = MOD(k,2) + 1
      DO i = i1, nxsd, 2
         sdio = (k-1)*nysd*nxsd + i
         sdii = sdio + reset
         IF (yrank/=0) THEN
           transfer_array(1:neqy, k,i) = psio(ieq:ieq+neqy-1,3,sdio)
           !!CALL MPI_SEND(psio(ieq,3,sdio),neqy,MPI_DOUBLE_PRECISION,yprv,yprtag,ycomm,ierr)
         END IF
       !!  IF (yrank/=npy-1) CALL MPI_RECV(psii(ieq,3,sdii),neqy,MPI_DOUBLE_PRECISION,ynxt,myytag,ycomm,istat,ierr)
       END DO
    END DO

    IF (yrank/=0)     CALL MPI_SEND(transfer_array,neqy*nzsd*nxsd,MPI_DOUBLE_PRECISION,yprv,yprtag,ycomm,ierr)
    IF (yrank/=npy-1) CALL MPI_RECV(transfer_array,neqy*nzsd*nxsd,MPI_DOUBLE_PRECISION,ynxt,myytag,ycomm,istat,ierr)

    DO k = 1, nzsd
       i1 = MOD(k,2) + 1
       DO i = i1, nxsd, 2
          sdio = (k-1)*nysd*nxsd + i
          sdii = sdio + reset
          IF (yrank/=npy-1) psii(ieq:ieq+neqy-1,3,sdii) = transfer_array(1:neqy, k,i)
          !!IF (yrank/=npy-1) CALL MPI_RECV(psii(ieq,3,sdii),neqy,MPI_DOUBLE_PRECISION,ynxt,myytag,ycomm,istat,ierr)
        END DO
     END DO
   !!------------------------------------------------------------------------------
   !!DN-Y, OCTANT 4----------------------------------------------------------------
   !!------------------------------------------------------------------------------
   DO k = 1, nzsd
      i1 = MOD(k,2) + 1
      DO i = i1, nxsd, 2
         sdio = (k-1)*nysd*nxsd + i
         sdii = sdio + reset
         IF (yrank/=0) THEN
           transfer_array(1:neqy, k,i) = psio(ieq:ieq+neqy-1,4,sdio)
           !!CALL MPI_SEND(psio(ieq,3,sdio),neqy,MPI_DOUBLE_PRECISION,yprv,yprtag,ycomm,ierr)
         END IF
       !!  IF (yrank/=npy-1) CALL MPI_RECV(psii(ieq,3,sdii),neqy,MPI_DOUBLE_PRECISION,ynxt,myytag,ycomm,istat,ierr)
       END DO
    END DO

    IF (yrank/=0)     CALL MPI_SEND(transfer_array,neqy*nzsd*nxsd,MPI_DOUBLE_PRECISION,yprv,yprtag,ycomm,ierr)
    IF (yrank/=npy-1) CALL MPI_RECV(transfer_array,neqy*nzsd*nxsd,MPI_DOUBLE_PRECISION,ynxt,myytag,ycomm,istat,ierr)

    DO k = 1, nzsd
       i1 = MOD(k,2) + 1
       DO i = i1, nxsd, 2
          sdio = (k-1)*nysd*nxsd + i
          sdii = sdio + reset
          IF (yrank/=npy-1) psii(ieq:ieq+neqy-1,4,sdii) = transfer_array(1:neqy, k,i)
          !!IF (yrank/=npy-1) CALL MPI_RECV(psii(ieq,3,sdii),neqy,MPI_DOUBLE_PRECISION,ynxt,myytag,ycomm,istat,ierr)
        END DO
     END DO
   !!------------------------------------------------------------------------------
   !!DN-Y, OCTANT 7----------------------------------------------------------------
   !!------------------------------------------------------------------------------
   DO k = 1, nzsd
      i1 = MOD(k,2) + 1
      DO i = i1, nxsd, 2
         sdio = (k-1)*nysd*nxsd + i
         sdii = sdio + reset
         IF (yrank/=0) THEN
           transfer_array(1:neqy, k,i) = psio(ieq:ieq+neqy-1,7,sdio)
           !!CALL MPI_SEND(psio(ieq,3,sdio),neqy,MPI_DOUBLE_PRECISION,yprv,yprtag,ycomm,ierr)
         END IF
       !!  IF (yrank/=npy-1) CALL MPI_RECV(psii(ieq,3,sdii),neqy,MPI_DOUBLE_PRECISION,ynxt,myytag,ycomm,istat,ierr)
       END DO
    END DO

    IF (yrank/=0)     CALL MPI_SEND(transfer_array,neqy*nzsd*nxsd,MPI_DOUBLE_PRECISION,yprv,yprtag,ycomm,ierr)
    IF (yrank/=npy-1) CALL MPI_RECV(transfer_array,neqy*nzsd*nxsd,MPI_DOUBLE_PRECISION,ynxt,myytag,ycomm,istat,ierr)

    DO k = 1, nzsd
       i1 = MOD(k,2) + 1
       DO i = i1, nxsd, 2
          sdio = (k-1)*nysd*nxsd + i
          sdii = sdio + reset
          IF (yrank/=npy-1) psii(ieq:ieq+neqy-1,7,sdii) = transfer_array(1:neqy, k,i)
          !!IF (yrank/=npy-1) CALL MPI_RECV(psii(ieq,3,sdii),neqy,MPI_DOUBLE_PRECISION,ynxt,myytag,ycomm,istat,ierr)
        END DO
     END DO
   !!------------------------------------------------------------------------------
   !!DN-Y, OCTANT 8----------------------------------------------------------------
   !!------------------------------------------------------------------------------
   DO k = 1, nzsd
      i1 = MOD(k,2) + 1
      DO i = i1, nxsd, 2
         sdio = (k-1)*nysd*nxsd + i
         sdii = sdio + reset
         IF (yrank/=0) THEN
           transfer_array(1:neqy, k,i) = psio(ieq:ieq+neqy-1,8,sdio)
           !!CALL MPI_SEND(psio(ieq,3,sdio),neqy,MPI_DOUBLE_PRECISION,yprv,yprtag,ycomm,ierr)
         END IF
       !!  IF (yrank/=npy-1) CALL MPI_RECV(psii(ieq,3,sdii),neqy,MPI_DOUBLE_PRECISION,ynxt,myytag,ycomm,istat,ierr)
       END DO
    END DO

    IF (yrank/=0)     CALL MPI_SEND(transfer_array,neqy*nzsd*nxsd,MPI_DOUBLE_PRECISION,yprv,yprtag,ycomm,ierr)
    IF (yrank/=npy-1) CALL MPI_RECV(transfer_array,neqy*nzsd*nxsd,MPI_DOUBLE_PRECISION,ynxt,myytag,ycomm,istat,ierr)

    DO k = 1, nzsd
       i1 = MOD(k,2) + 1
       DO i = i1, nxsd, 2
          sdio = (k-1)*nysd*nxsd + i
          sdii = sdio + reset
          IF (yrank/=npy-1) psii(ieq:ieq+neqy-1,8,sdii) = transfer_array(1:neqy, k,i)
          !!IF (yrank/=npy-1) CALL MPI_RECV(psii(ieq,3,sdii),neqy,MPI_DOUBLE_PRECISION,ynxt,myytag,ycomm,istat,ierr)
        END DO
     END DO

   DEALLOCATE(transfer_array)

   ! Reflective boundary conditions
   IF (yrank == 0 .AND. bc(3) == 1) THEN               ! Reflective at low y
      DO k = 1, nzsd                                   ! All y = 1 sub-domains
         sdii = (k-1)*nxsd*nysd
         DO i = 1, nxsd
            psii(ieq:neqz+neqy,1,sdii+i) = psio(ieq:neqz+neqy,4,sdii+i)
            psii(ieq:neqz+neqy,2,sdii+i) = psio(ieq:neqz+neqy,3,sdii+i)
            psii(ieq:neqz+neqy,5,sdii+i) = psio(ieq:neqz+neqy,8,sdii+i)
            psii(ieq:neqz+neqy,6,sdii+i) = psio(ieq:neqz+neqy,7,sdii+i)
         END DO
      END DO
   END IF
   IF (yrank == npy-1 .AND. bc(4) == 1) THEN            ! Reflective at high y
      DO k = 1, nzsd                                   ! All y = nysd sub-domains
         sdii = k*nxsd*nysd - nxsd
         DO i = 1, nxsd
            psii(ieq:neqz+neqy,4,sdii+i) = psio(ieq:neqz+neqy,1,sdii+i)
            psii(ieq:neqz+neqy,3,sdii+i) = psio(ieq:neqz+neqy,2,sdii+i)
            psii(ieq:neqz+neqy,8,sdii+i) = psio(ieq:neqz+neqy,5,sdii+i)
            psii(ieq:neqz+neqy,7,sdii+i) = psio(ieq:neqz+neqy,6,sdii+i)
         END DO
      END DO
   END IF

   !!------------------------------------------------------------------------------
   !!------------------------------------------------------------------------------
   !!UP-X--------------------------------------------------------------------------
   !!------------------------------------------------------------------------------
   !!------------------------------------------------------------------------------

   ieq = ieq + neqy
   neqx = apo*yzs

   !!Alocate transfer array to be neqz by number of face subdomains
   ALLOCATE(transfer_array(neqy, nzsd, nysd))
   transfer_array=0

   !!------------------------------------------------------------------------------
   !!UP-X, INTERNAL----------------------------------------------------------------
   !!------------------------------------------------------------------------------
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
   !!------------------------------------------------------------------------------
   !!UP-X, OCTANT 1----------------------------------------------------------------
   !!------------------------------------------------------------------------------
   reset = nxsd-1
   DO k = 1, nzsd
      i1 = MOD(k+1,2) + 1
      DO j = i1, nysd, 2
         sdio = (k-1)*nysd*nxsd + j*nxsd
         sdii = sdio - reset
         IF (xrank/=npx-1) THEN
           transfer_array(1:neqx,k,j) = psio(ieq:ieq+neqx-1,1,sdio)
           !!CALL MPI_SEND(psio(ieq,1,sdio),neqx,MPI_DOUBLE_PRECISION,xnxt,xnxtag,xcomm,ierr)
         END IF
       !!IF (xrank/=0) CALL MPI_RECV(psii(ieq,1,sdii),neqx,MPI_DOUBLE_PRECISION,xprv,myxtag,xcomm,istat,ierr)
     END DO
   END DO

   IF (xrank/=npx-1) CALL MPI_SEND(transfer_array,neqx * nzsd * nysd,MPI_DOUBLE_PRECISION,xnxt,xnxtag,xcomm,ierr)
   IF (xrank/=0)     CALL MPI_RECV(transfer_array,neqx * nzsd * nysd,MPI_DOUBLE_PRECISION,xprv,myxtag,xcomm,istat,ierr)

   DO k = 1, nzsd
      i1 = MOD(k+1,2) + 1
      DO j = i1, nysd, 2
         sdio = (k-1)*nysd*nxsd + j*nxsd
         sdii = sdio - reset
         IF (xrank/=0) psii(ieq:ieq+neqx-1,1,sdii) = transfer_array(1:neqx,k,j)
       !!IF (xrank/=0) CALL MPI_RECV(psii(ieq,1,sdii),neqx,MPI_DOUBLE_PRECISION,xprv,myxtag,xcomm,istat,ierr)
     END DO
   END DO


   !!------------------------------------------------------------------------------
   !!UP-X, OCTANT 4----------------------------------------------------------------
   !!------------------------------------------------------------------------------
   DO k = 1, nzsd
      i1 = MOD(k+1,2) + 1
      DO j = i1, nysd, 2
         sdio = (k-1)*nysd*nxsd + j*nxsd
         sdii = sdio - reset
         IF (xrank/=npx-1) THEN
           transfer_array(1:neqx,k,j) = psio(ieq:ieq+neqx-1,4,sdio)
           !!CALL MPI_SEND(psio(ieq,1,sdio),neqx,MPI_DOUBLE_PRECISION,xnxt,xnxtag,xcomm,ierr)
         END IF
       !!IF (xrank/=0) CALL MPI_RECV(psii(ieq,1,sdii),neqx,MPI_DOUBLE_PRECISION,xprv,myxtag,xcomm,istat,ierr)
     END DO
   END DO

   IF (xrank/=npx-1) CALL MPI_SEND(transfer_array,neqx * nzsd * nysd,MPI_DOUBLE_PRECISION,xnxt,xnxtag,xcomm,ierr)
   IF (xrank/=0)     CALL MPI_RECV(transfer_array,neqx * nzsd * nysd,MPI_DOUBLE_PRECISION,xprv,myxtag,xcomm,istat,ierr)

   DO k = 1, nzsd
      i1 = MOD(k+1,2) + 1
      DO j = i1, nysd, 2
         sdio = (k-1)*nysd*nxsd + j*nxsd
         sdii = sdio - reset
         IF (xrank/=0) psii(ieq:ieq+neqx-1,4,sdii) = transfer_array(1:neqx,k,j)
       !!IF (xrank/=0) CALL MPI_RECV(psii(ieq,1,sdii),neqx,MPI_DOUBLE_PRECISION,xprv,myxtag,xcomm,istat,ierr)
     END DO
   END DO

   !!------------------------------------------------------------------------------
   !!UP-X, OCTANT 5----------------------------------------------------------------
   !!------------------------------------------------------------------------------
   DO k = 1, nzsd
      i1 = MOD(k+1,2) + 1
      DO j = i1, nysd, 2
         sdio = (k-1)*nysd*nxsd + j*nxsd
         sdii = sdio - reset
         IF (xrank/=npx-1) THEN
           transfer_array(1:neqx,k,j) = psio(ieq:ieq+neqx-1,5,sdio)
           !!CALL MPI_SEND(psio(ieq,1,sdio),neqx,MPI_DOUBLE_PRECISION,xnxt,xnxtag,xcomm,ierr)
         END IF
       !!IF (xrank/=0) CALL MPI_RECV(psii(ieq,1,sdii),neqx,MPI_DOUBLE_PRECISION,xprv,myxtag,xcomm,istat,ierr)
     END DO
   END DO

   IF (xrank/=npx-1) CALL MPI_SEND(transfer_array,neqx * nzsd * nysd,MPI_DOUBLE_PRECISION,xnxt,xnxtag,xcomm,ierr)
   IF (xrank/=0)     CALL MPI_RECV(transfer_array,neqx * nzsd * nysd,MPI_DOUBLE_PRECISION,xprv,myxtag,xcomm,istat,ierr)

   DO k = 1, nzsd
      i1 = MOD(k+1,2) + 1
      DO j = i1, nysd, 2
         sdio = (k-1)*nysd*nxsd + j*nxsd
         sdii = sdio - reset
         IF (xrank/=0) psii(ieq:ieq+neqx-1,5,sdii) = transfer_array(1:neqx,k,j)
       !!IF (xrank/=0) CALL MPI_RECV(psii(ieq,1,sdii),neqx,MPI_DOUBLE_PRECISION,xprv,myxtag,xcomm,istat,ierr)
     END DO
   END DO
   !!------------------------------------------------------------------------------
   !!UP-X, OCTANT 8----------------------------------------------------------------
   !!------------------------------------------------------------------------------
   DO k = 1, nzsd
      i1 = MOD(k+1,2) + 1
      DO j = i1, nysd, 2
         sdio = (k-1)*nysd*nxsd + j*nxsd
         sdii = sdio - reset
         IF (xrank/=npx-1) THEN
           transfer_array(1:neqx,k,j) = psio(ieq:ieq+neqx-1,8,sdio)
           !!CALL MPI_SEND(psio(ieq,1,sdio),neqx,MPI_DOUBLE_PRECISION,xnxt,xnxtag,xcomm,ierr)
         END IF
       !!IF (xrank/=0) CALL MPI_RECV(psii(ieq,1,sdii),neqx,MPI_DOUBLE_PRECISION,xprv,myxtag,xcomm,istat,ierr)
     END DO
   END DO

   IF (xrank/=npx-1) CALL MPI_SEND(transfer_array,neqx * nzsd * nysd,MPI_DOUBLE_PRECISION,xnxt,xnxtag,xcomm,ierr)
   IF (xrank/=0)     CALL MPI_RECV(transfer_array,neqx * nzsd * nysd,MPI_DOUBLE_PRECISION,xprv,myxtag,xcomm,istat,ierr)

   DO k = 1, nzsd
      i1 = MOD(k+1,2) + 1
      DO j = i1, nysd, 2
         sdio = (k-1)*nysd*nxsd + j*nxsd
         sdii = sdio - reset
         IF (xrank/=0) psii(ieq:ieq+neqx-1,8,sdii) = transfer_array(1:neqx,k,j)
       !!IF (xrank/=0) CALL MPI_RECV(psii(ieq,1,sdii),neqx,MPI_DOUBLE_PRECISION,xprv,myxtag,xcomm,istat,ierr)
     END DO
   END DO


   !!------------------------------------------------------------------------------
   !!------------------------------------------------------------------------------
   !!DN-X--------------------------------------------------------------------------
   !!------------------------------------------------------------------------------
   !!------------------------------------------------------------------------------


   !!------------------------------------------------------------------------------
   !!DN-X, INTERNAL----------------------------------------------------------------
   !!------------------------------------------------------------------------------
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
   !!------------------------------------------------------------------------------
   !!DN-X, OCTANT 2----------------------------------------------------------------
   !!------------------------------------------------------------------------------
   DO k = 1, nzsd
      i1 = MOD(k,2) + 1
      DO j = i1, nysd, 2
         sdii = (k-1)*nysd*nxsd + j*nxsd
         sdio = sdii - reset
         IF (xrank/=0) THEN
           transfer_array(1:neqx,k,j) = psio(ieq:ieq+neqx-1,2,sdio)
           !!CALL MPI_SEND(psio(ieq,2,sdio),neqx,MPI_DOUBLE_PRECISION,xprv,xprtag,xcomm,ierr)
         END IF
         !!IF (xrank/=npx-1) CALL MPI_RECV(psii(ieq,2,sdii),neqx,MPI_DOUBLE_PRECISION,xnxt,myxtag,xcomm,istat,ierr)
       END DO
    END DO
    IF (xrank/=0)     CALL MPI_SEND(transfer_array,neqx * nzsd * nysd,MPI_DOUBLE_PRECISION,xprv,xprtag,xcomm,ierr)
    IF (xrank/=npx-1) CALL MPI_RECV(transfer_array,neqx * nzsd * nysd,MPI_DOUBLE_PRECISION,xnxt,myxtag,xcomm,istat,ierr)

    DO k = 1, nzsd
       i1 = MOD(k,2) + 1
       DO j = i1, nysd, 2
          sdii = (k-1)*nysd*nxsd + j*nxsd
          sdio = sdii - reset
          IF (xrank/=npx-1) psii(ieq:ieq+neqx-1,2,sdii) = transfer_array(1:neqx,k,j)
          !!IF (xrank/=npx-1) CALL MPI_RECV(psii(ieq,2,sdii),neqx,MPI_DOUBLE_PRECISION,xnxt,myxtag,xcomm,istat,ierr)
        END DO
     END DO
   !!------------------------------------------------------------------------------
   !!DN-X, OCTANT 3----------------------------------------------------------------
   !!------------------------------------------------------------------------------
   DO k = 1, nzsd
      i1 = MOD(k,2) + 1
      DO j = i1, nysd, 2
         sdii = (k-1)*nysd*nxsd + j*nxsd
         sdio = sdii - reset
         IF (xrank/=0) THEN
           transfer_array(1:neqx,k,j) = psio(ieq:ieq+neqx-1,3,sdio)
           !!CALL MPI_SEND(psio(ieq,2,sdio),neqx,MPI_DOUBLE_PRECISION,xprv,xprtag,xcomm,ierr)
         END IF
         !!IF (xrank/=npx-1) CALL MPI_RECV(psii(ieq,2,sdii),neqx,MPI_DOUBLE_PRECISION,xnxt,myxtag,xcomm,istat,ierr)
       END DO
    END DO
    IF (xrank/=0)     CALL MPI_SEND(transfer_array,neqx * nzsd * nysd,MPI_DOUBLE_PRECISION,xprv,xprtag,xcomm,ierr)
    IF (xrank/=npx-1) CALL MPI_RECV(transfer_array,neqx * nzsd * nysd,MPI_DOUBLE_PRECISION,xnxt,myxtag,xcomm,istat,ierr)

    DO k = 1, nzsd
       i1 = MOD(k,2) + 1
       DO j = i1, nysd, 2
          sdii = (k-1)*nysd*nxsd + j*nxsd
          sdio = sdii - reset
          IF (xrank/=npx-1) psii(ieq:ieq+neqx-1,3,sdii) = transfer_array(1:neqx,k,j)
          !!IF (xrank/=npx-1) CALL MPI_RECV(psii(ieq,2,sdii),neqx,MPI_DOUBLE_PRECISION,xnxt,myxtag,xcomm,istat,ierr)
        END DO
     END DO
   !!------------------------------------------------------------------------------
   !!DN-X, OCTANT 6----------------------------------------------------------------
   !!------------------------------------------------------------------------------
   DO k = 1, nzsd
      i1 = MOD(k,2) + 1
      DO j = i1, nysd, 2
         sdii = (k-1)*nysd*nxsd + j*nxsd
         sdio = sdii - reset
         IF (xrank/=0) THEN
           transfer_array(1:neqx,k,j) = psio(ieq:ieq+neqx-1,6,sdio)
           !!CALL MPI_SEND(psio(ieq,2,sdio),neqx,MPI_DOUBLE_PRECISION,xprv,xprtag,xcomm,ierr)
         END IF
         !!IF (xrank/=npx-1) CALL MPI_RECV(psii(ieq,2,sdii),neqx,MPI_DOUBLE_PRECISION,xnxt,myxtag,xcomm,istat,ierr)
       END DO
    END DO
    IF (xrank/=0)     CALL MPI_SEND(transfer_array,neqx * nzsd * nysd,MPI_DOUBLE_PRECISION,xprv,xprtag,xcomm,ierr)
    IF (xrank/=npx-1) CALL MPI_RECV(transfer_array,neqx * nzsd * nysd,MPI_DOUBLE_PRECISION,xnxt,myxtag,xcomm,istat,ierr)

    DO k = 1, nzsd
       i1 = MOD(k,2) + 1
       DO j = i1, nysd, 2
          sdii = (k-1)*nysd*nxsd + j*nxsd
          sdio = sdii - reset
          IF (xrank/=npx-1) psii(ieq:ieq+neqx-1,6,sdii) = transfer_array(1:neqx,k,j)
          !!IF (xrank/=npx-1) CALL MPI_RECV(psii(ieq,2,sdii),neqx,MPI_DOUBLE_PRECISION,xnxt,myxtag,xcomm,istat,ierr)
        END DO
     END DO
   !!------------------------------------------------------------------------------
   !!DN-X, OCTANT 7----------------------------------------------------------------
   !!------------------------------------------------------------------------------
   DO k = 1, nzsd
      i1 = MOD(k,2) + 1
      DO j = i1, nysd, 2
         sdii = (k-1)*nysd*nxsd + j*nxsd
         sdio = sdii - reset
         IF (xrank/=0) THEN
           transfer_array(1:neqx,k,j) = psio(ieq:ieq+neqx-1,7,sdio)
           !!CALL MPI_SEND(psio(ieq,2,sdio),neqx,MPI_DOUBLE_PRECISION,xprv,xprtag,xcomm,ierr)
         END IF
         !!IF (xrank/=npx-1) CALL MPI_RECV(psii(ieq,2,sdii),neqx,MPI_DOUBLE_PRECISION,xnxt,myxtag,xcomm,istat,ierr)
       END DO
    END DO
    IF (xrank/=0)     CALL MPI_SEND(transfer_array,neqx * nzsd * nysd,MPI_DOUBLE_PRECISION,xprv,xprtag,xcomm,ierr)
    IF (xrank/=npx-1) CALL MPI_RECV(transfer_array,neqx * nzsd * nysd,MPI_DOUBLE_PRECISION,xnxt,myxtag,xcomm,istat,ierr)

    DO k = 1, nzsd
       i1 = MOD(k,2) + 1
       DO j = i1, nysd, 2
          sdii = (k-1)*nysd*nxsd + j*nxsd
          sdio = sdii - reset
          IF (xrank/=npx-1) psii(ieq:ieq+neqx-1,7,sdii) = transfer_array(1:neqx,k,j)
          !!IF (xrank/=npx-1) CALL MPI_RECV(psii(ieq,2,sdii),neqx,MPI_DOUBLE_PRECISION,xnxt,myxtag,xcomm,istat,ierr)
        END DO
     END DO

   DEALLOCATE(transfer_array)

   ! Reflective boundary conditions
   IF (xrank == 0 .AND. bc(5) == 1) THEN               ! Reflective at low x
      DO i = 1, nsdp, nxsd                             ! All x = 1 sub-domains
         psii(ieq:bcs,1,i) = psio(ieq:bcs,2,i)
         psii(ieq:bcs,4,i) = psio(ieq:bcs,3,i)
         psii(ieq:bcs,5,i) = psio(ieq:bcs,6,i)
         psii(ieq:bcs,8,i) = psio(ieq:bcs,7,i)
      END DO
   END IF
   IF (xrank == npx-1 .AND. bc(6) == 1) THEN           ! Reflective at low x
      DO i = nxsd, nsdp, nxsd                          ! All x = 1 sub-domains
         psii(ieq:bcs,2,i) = psio(ieq:bcs,1,i)
         psii(ieq:bcs,3,i) = psio(ieq:bcs,4,i)
         psii(ieq:bcs,6,i) = psio(ieq:bcs,5,i)
         psii(ieq:bcs,7,i) = psio(ieq:bcs,8,i)
      END DO
   END IF
   !--------------------------------------------------------------------------------------

ELSE IF (tmpr < err) THEN
   bcnvf = 1
   tsolve = MPI_WTIME()
   IF (irank == root) THEN
      WRITE(8,'(2X,A,I5,A,ES11.3,/)') "System converged in ", bit, " iterations. Max error = ", tmpr !dfmx
  END IF
ELSE IF (bit == bitmx) THEN
   bcnvf = 0
   tsolve = MPI_WTIME()
   IF (irank == root) THEN
      WRITE (8,'(/,2X,A,I5,A,ES11.3)') "ERROR: System did not converge in ", bit, " iterations . Max error= ", tmpr !dfmx
      WRITE (8,*) "Will adversely affect other results. Aborting program."
   END IF
END IF

111 FORMAT(2X,'It ', I5,' Dfmx ', ES11.3,' phi, res ', 2ES12.3)

RETURN
END SUBROUTINE pgsblk
