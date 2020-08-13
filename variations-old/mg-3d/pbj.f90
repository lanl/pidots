SUBROUTINE pbj(v,spx,spy,spz,bit,itmx,cnvf)

!-------------------------------------------------------------
!
!  Parallel Block Jacobi
!
!  Use the scalar flux solution to check for convergence
!  If not converged, move to next iteration by passing angular
!   flux moments to neighbors. Neighbors store in psi* vectors
!  If not converged and out of iterations, store everything
!   and exit this group
!  If converged, store everything and move to next group
!
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER, INTENT(IN) :: v, spx, spy, spz, bit, itmx
INTEGER, INTENT(OUT) :: cnvf
INTEGER :: myztag, myytag, myxtag, znxt, zprv, ynxt, yprv, xnxt, xprv
INTEGER :: znxtag, zprtag, ynxtag, yprtag, xnxtag, xprtag, neqz, neqy, neqx, ieq, ierr
INTEGER, DIMENSION(3) :: istat
REAL*8 :: df, dfmx

INCLUDE 'mpif.h'

! Set up some tags
myztag = 100 + zrank
myytag = 100 + yrank
myxtag = 100 + xrank
znxt = zrank + spz
znxtag = 100 + znxt
zprv = zrank - spz
zprtag = 100 + zprv
ynxt = yrank + spy
ynxtag = 100 + ynxt
yprv = yrank - spy
yprtag = 100 + yprv
xnxt = xrank + spx
xnxtag = 100 + xnxt
xprv = xrank - spx
xprtag = 100 + xprv

!-----------------------------------------------------------------------------------
! Check convergence of scalar flux if at finest grid
IF (v == 1) THEN
   ! Each process checks its own fluxes, then perform MPI Reduction
   IF (MINVAL(ABS(phiold)) >= tolr) THEN
      df = MAXVAL(ABS((phi(:,v) - phiold)/phiold))
   ELSE
      df = MAXVAL(ABS(phi(:,v) - phiold))
   END IF

   IF (isize == 1) THEN
      IF (bc(1)/=1 .AND. bc(2)/=1 .AND. bc(3)/=1 .AND. bc(4)/=1 .AND. bc(5)/=1 .AND. bc(6)/=1) THEN
         df = 0.0
      END IF
   END IF

   ! Get the overall max difference
   CALL MPI_ALLREDUCE(df,dfmx,1,MPI_DOUBLE_PRECISION,MPI_MAX,allcomm,ierr)

   ! Only print iteration information if requested
   IF (irank == root) THEN
      IF (itp == 1) WRITE(8,111) bit, dfmx
   END IF

   ! Set cnvf to 1 if converged
   IF (dfmx < err) cnvf = 1

   ! Set the previous iterate of the fine grid flux equal to the current
   phiold = phi(:,v)
END IF
111 FORMAT(2X,'It ', I5,' Dfmx ', ES11.3)

!-----------------------------------------------------------------------------------------------
! Get copy of old psii for computing the residual
IF (v == 1 .AND. bit == itmx) THEN   ! Only 1 iteration at v>1 stages, so oldpsi is already zero
!IF (bit == itmx) THEN
   oldpsi = psii(:,:,v)
END IF

!-----------------------------------------------------------------------------------------------
! Send psio to neighbors. Inner blocks have 24 communications (8 octants to
! 3 neigbors each). Then set the incoming data into psii

! Perform the sends 'up z'
neqz = apo*xys
! Octant 1
IF (zrank /= npz-spz) CALL MPI_SEND(psio(1,1,v),neqz,MPI_DOUBLE_PRECISION,znxt,znxtag,zcomm,ierr)
IF (zrank /= 0) CALL MPI_RECV(psii(1,1,v),neqz,MPI_DOUBLE_PRECISION,zprv,myztag,zcomm,istat,ierr)
! Octant 2
IF (zrank /= npz-spz) CALL MPI_SEND(psio(1,2,v),neqz,MPI_DOUBLE_PRECISION,znxt,znxtag,zcomm,ierr)
IF (zrank /= 0) CALL MPI_RECV(psii(1,2,v),neqz,MPI_DOUBLE_PRECISION,zprv,myztag,zcomm,istat,ierr)
! Octant 3
IF (zrank /= npz-spz) CALL MPI_SEND(psio(1,3,v),neqz,MPI_DOUBLE_PRECISION,znxt,znxtag,zcomm,ierr)
IF (zrank /= 0) CALL MPI_RECV(psii(1,3,v),neqz,MPI_DOUBLE_PRECISION,zprv,myztag,zcomm,istat,ierr)
! Octant 4
IF (zrank /= npz-spz) CALL MPI_SEND(psio(1,4,v),neqz,MPI_DOUBLE_PRECISION,znxt,znxtag,zcomm,ierr)
IF (zrank /= 0) CALL MPI_RECV(psii(1,4,v),neqz,MPI_DOUBLE_PRECISION,zprv,myztag,zcomm,istat,ierr)

! Perform the sends 'down z'
! Octant 5
IF (zrank /= 0) CALL MPI_SEND(psio(1,5,v),neqz,MPI_DOUBLE_PRECISION,zprv,zprtag,zcomm,ierr)
IF (zrank /= npz-spz) CALL MPI_RECV(psii(1,5,v),neqz,MPI_DOUBLE_PRECISION,znxt,myztag,zcomm,istat,ierr)
! Octant 6
IF (zrank /= 0) CALL MPI_SEND(psio(1,6,v),neqz,MPI_DOUBLE_PRECISION,zprv,zprtag,zcomm,ierr)
IF (zrank /= npz-spz) CALL MPI_RECV(psii(1,6,v),neqz,MPI_DOUBLE_PRECISION,znxt,myztag,zcomm,istat,ierr)
! Octant 7
IF (zrank /= 0) CALL MPI_SEND(psio(1,7,v),neqz,MPI_DOUBLE_PRECISION,zprv,zprtag,zcomm,ierr)
IF (zrank /= npz-spz) CALL MPI_RECV(psii(1,7,v),neqz,MPI_DOUBLE_PRECISION,znxt,myztag,zcomm,istat,ierr)
! Octant 8
IF (zrank /= 0) CALL MPI_SEND(psio(1,8,v),neqz,MPI_DOUBLE_PRECISION,zprv,zprtag,zcomm,ierr)
IF (zrank /= npz-spz) CALL MPI_RECV(psii(1,8,v),neqz,MPI_DOUBLE_PRECISION,znxt,myztag,zcomm,istat,ierr)

!-------------------------------------------------------------------------------------------------
! Perform the sends 'up y'
ieq = neqz + 1
neqy = apo*xzs
! Quadrant 1
IF (yrank /= npy-spy) CALL MPI_SEND(psio(ieq,1,v),neqy,MPI_DOUBLE_PRECISION,ynxt,ynxtag,ycomm,ierr)
IF (yrank /= 0) CALL MPI_RECV(psii(ieq,1,v),neqy,MPI_DOUBLE_PRECISION,yprv,myytag,ycomm,istat,ierr)
! Quadrant 2
IF (yrank /= npy-spy) CALL MPI_SEND(psio(ieq,2,v),neqy,MPI_DOUBLE_PRECISION,ynxt,ynxtag,ycomm,ierr)
IF (yrank /= 0) CALL MPI_RECV(psii(ieq,2,v),neqy,MPI_DOUBLE_PRECISION,yprv,myytag,ycomm,istat,ierr)
! Quadrant 5
IF (yrank /= npy-spy) CALL MPI_SEND(psio(ieq,5,v),neqy,MPI_DOUBLE_PRECISION,ynxt,ynxtag,ycomm,ierr)
IF (yrank /= 0) CALL MPI_RECV(psii(ieq,5,v),neqy,MPI_DOUBLE_PRECISION,yprv,myytag,ycomm,istat,ierr)
! Quadrant 6
IF (yrank /= npy-spy) CALL MPI_SEND(psio(ieq,6,v),neqy,MPI_DOUBLE_PRECISION,ynxt,ynxtag,ycomm,ierr)
IF (yrank /= 0) CALL MPI_RECV(psii(ieq,6,v),neqy,MPI_DOUBLE_PRECISION,yprv,myytag,ycomm,istat,ierr)

! Perform the sends 'down y'
! Quadrant 3
IF (yrank /= 0) CALL MPI_SEND(psio(ieq,3,v),neqy,MPI_DOUBLE_PRECISION,yprv,yprtag,ycomm,ierr)
IF (yrank /= npy-spy) CALL MPI_RECV(psii(ieq,3,v),neqy,MPI_DOUBLE_PRECISION,ynxt,myytag,ycomm,istat,ierr)
! Quadrant 4
IF (yrank /= 0) CALL MPI_SEND(psio(ieq,4,v),neqy,MPI_DOUBLE_PRECISION,yprv,yprtag,ycomm,ierr)
IF (yrank /= npy-spy) CALL MPI_RECV(psii(ieq,4,v),neqy,MPI_DOUBLE_PRECISION,ynxt,myytag,ycomm,istat,ierr)
! Quadrant 7
IF (yrank /= 0) CALL MPI_SEND(psio(ieq,7,v),neqy,MPI_DOUBLE_PRECISION,yprv,yprtag,ycomm,ierr)
IF (yrank /= npy-spy) CALL MPI_RECV(psii(ieq,7,v),neqy,MPI_DOUBLE_PRECISION,ynxt,myytag,ycomm,istat,ierr)
! Quadrant 8
IF (yrank /= 0) CALL MPI_SEND(psio(ieq,8,v),neqy,MPI_DOUBLE_PRECISION,yprv,yprtag,ycomm,ierr)
IF (yrank /= npy-spy) CALL MPI_RECV(psii(ieq,8,v),neqy,MPI_DOUBLE_PRECISION,ynxt,myytag,ycomm,istat,ierr)

!--------------------------------------------------------------------------------------------------
! Perform the sends 'up x'
ieq = ieq + neqy
neqx = apo*yzs
! Quadrant 1
IF (xrank /= npx-spx) CALL MPI_SEND(psio(ieq,1,v),neqx,MPI_DOUBLE_PRECISION,xnxt,xnxtag,xcomm,ierr)
IF (xrank /= 0) CALL MPI_RECV(psii(ieq,1,v),neqx,MPI_DOUBLE_PRECISION,xprv,myxtag,xcomm,istat,ierr)
! Quadrant 4
IF (xrank /= npx-spx) CALL MPI_SEND(psio(ieq,4,v),neqx,MPI_DOUBLE_PRECISION,xnxt,xnxtag,xcomm,ierr)
IF (xrank /= 0) CALL MPI_RECV(psii(ieq,4,v),neqx,MPI_DOUBLE_PRECISION,xprv,myxtag,xcomm,istat,ierr)
! Quadrant 5
IF (xrank /= npx-spx) CALL MPI_SEND(psio(ieq,5,v),neqx,MPI_DOUBLE_PRECISION,xnxt,xnxtag,xcomm,ierr)
IF (xrank /= 0) CALL MPI_RECV(psii(ieq,5,v),neqx,MPI_DOUBLE_PRECISION,xprv,myxtag,xcomm,istat,ierr)
! Quadrant 8
IF (xrank /= npx-spx) CALL MPI_SEND(psio(ieq,8,v),neqx,MPI_DOUBLE_PRECISION,xnxt,xnxtag,xcomm,ierr)
IF (xrank /= 0) CALL MPI_RECV(psii(ieq,8,v),neqx,MPI_DOUBLE_PRECISION,xprv,myxtag,xcomm,istat,ierr)

! Perform the sends 'down x'
! Quadrant 2
IF (xrank /= 0) CALL MPI_SEND(psio(ieq,2,v),neqx,MPI_DOUBLE_PRECISION,xprv,xprtag,xcomm,ierr)
IF (xrank /= npx-spx) CALL MPI_RECV(psii(ieq,2,v),neqx,MPI_DOUBLE_PRECISION,xnxt,myxtag,xcomm,istat,ierr)
! Quadrant 3
IF (xrank /= 0) CALL MPI_SEND(psio(ieq,3,v),neqx,MPI_DOUBLE_PRECISION,xprv,xprtag,xcomm,ierr)
IF (xrank /= npx-spx) CALL MPI_RECV(psii(ieq,3,v),neqx,MPI_DOUBLE_PRECISION,xnxt,myxtag,xcomm,istat,ierr)
! Quadrant 6
IF (xrank /= 0) CALL MPI_SEND(psio(ieq,6,v),neqx,MPI_DOUBLE_PRECISION,xprv,xprtag,xcomm,ierr)
IF (xrank /= npx-spx) CALL MPI_RECV(psii(ieq,6,v),neqx,MPI_DOUBLE_PRECISION,xnxt,myxtag,xcomm,istat,ierr)
! Quadrant 7
IF (xrank /= 0) CALL MPI_SEND(psio(ieq,7,v),neqx,MPI_DOUBLE_PRECISION,xprv,xprtag,xcomm,ierr)
IF (xrank /= npx-spx) CALL MPI_RECV(psii(ieq,7,v),neqx,MPI_DOUBLE_PRECISION,xnxt,myxtag,xcomm,istat,ierr)

RETURN
END SUBROUTINE pbj
