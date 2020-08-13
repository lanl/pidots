SUBROUTINE pbj(spx,spy,spz,fpi,fpo)

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
IMPLICIT NONE
INTEGER, INTENT(IN) :: spx, spy, spz
INTEGER :: myztag, myytag, myxtag, znxt, zprv, ynxt, yprv, xnxt, xprv
INTEGER :: znxtag, zprtag, ynxtag, yprtag, xnxtag, xprtag, neqz, neqy, neqx, ieq, ierr
INTEGER, DIMENSION(3) :: istat
REAL*8, DIMENSION(bcs,8), INTENT(INOUT) :: fpi
REAL*8, DIMENSION(bcs,8), INTENT(IN) :: fpo

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

!-----------------------------------------------------------------------------------------------
! Get copy of old fpi for computing the residual
! Only 1 iteration at v>1 stages, so foldpsi is already zero
! IF (bit == itmx) THEN
!   foldpsi = fpi
! END IF

!-----------------------------------------------------------------------------------------------
! Send fpo to neighbors. Inner blocks have 24 communications (8 octants to
! 3 neigbors each). Then set the incoming data into fpi

! Perform the sends 'up z'
neqz = apo*xys
! Octant 1
IF (zrank /= npz-spz) CALL MPI_SEND(fpo(1,1),neqz,MPI_DOUBLE_PRECISION,znxt,znxtag,zcomm,ierr)
IF (zrank /= 0) CALL MPI_RECV(fpi(1,1),neqz,MPI_DOUBLE_PRECISION,zprv,myztag,zcomm,istat,ierr)
! Octant 2
IF (zrank /= npz-spz) CALL MPI_SEND(fpo(1,2),neqz,MPI_DOUBLE_PRECISION,znxt,znxtag,zcomm,ierr)
IF (zrank /= 0) CALL MPI_RECV(fpi(1,2),neqz,MPI_DOUBLE_PRECISION,zprv,myztag,zcomm,istat,ierr)
! Octant 3
IF (zrank /= npz-spz) CALL MPI_SEND(fpo(1,3),neqz,MPI_DOUBLE_PRECISION,znxt,znxtag,zcomm,ierr)
IF (zrank /= 0) CALL MPI_RECV(fpi(1,3),neqz,MPI_DOUBLE_PRECISION,zprv,myztag,zcomm,istat,ierr)
! Octant 4
IF (zrank /= npz-spz) CALL MPI_SEND(fpo(1,4),neqz,MPI_DOUBLE_PRECISION,znxt,znxtag,zcomm,ierr)
IF (zrank /= 0) CALL MPI_RECV(fpi(1,4),neqz,MPI_DOUBLE_PRECISION,zprv,myztag,zcomm,istat,ierr)

! Perform the sends 'down z'
! Octant 5
IF (zrank /= 0) CALL MPI_SEND(fpo(1,5),neqz,MPI_DOUBLE_PRECISION,zprv,zprtag,zcomm,ierr)
IF (zrank /= npz-spz) CALL MPI_RECV(fpi(1,5),neqz,MPI_DOUBLE_PRECISION,znxt,myztag,zcomm,istat,ierr)
! Octant 6
IF (zrank /= 0) CALL MPI_SEND(fpo(1,6),neqz,MPI_DOUBLE_PRECISION,zprv,zprtag,zcomm,ierr)
IF (zrank /= npz-spz) CALL MPI_RECV(fpi(1,6),neqz,MPI_DOUBLE_PRECISION,znxt,myztag,zcomm,istat,ierr)
! Octant 7
IF (zrank /= 0) CALL MPI_SEND(fpo(1,7),neqz,MPI_DOUBLE_PRECISION,zprv,zprtag,zcomm,ierr)
IF (zrank /= npz-spz) CALL MPI_RECV(fpi(1,7),neqz,MPI_DOUBLE_PRECISION,znxt,myztag,zcomm,istat,ierr)
! Octant 8
IF (zrank /= 0) CALL MPI_SEND(fpo(1,8),neqz,MPI_DOUBLE_PRECISION,zprv,zprtag,zcomm,ierr)
IF (zrank /= npz-spz) CALL MPI_RECV(fpi(1,8),neqz,MPI_DOUBLE_PRECISION,znxt,myztag,zcomm,istat,ierr)

!-------------------------------------------------------------------------------------------------
! Perform the sends 'up y'
ieq = neqz + 1
neqy = apo*xzs
! Quadrant 1
IF (yrank /= npy-spy) CALL MPI_SEND(fpo(ieq,1),neqy,MPI_DOUBLE_PRECISION,ynxt,ynxtag,ycomm,ierr)
IF (yrank /= 0) CALL MPI_RECV(fpi(ieq,1),neqy,MPI_DOUBLE_PRECISION,yprv,myytag,ycomm,istat,ierr)
! Quadrant 2
IF (yrank /= npy-spy) CALL MPI_SEND(fpo(ieq,2),neqy,MPI_DOUBLE_PRECISION,ynxt,ynxtag,ycomm,ierr)
IF (yrank /= 0) CALL MPI_RECV(fpi(ieq,2),neqy,MPI_DOUBLE_PRECISION,yprv,myytag,ycomm,istat,ierr)
! Quadrant 5
IF (yrank /= npy-spy) CALL MPI_SEND(fpo(ieq,5),neqy,MPI_DOUBLE_PRECISION,ynxt,ynxtag,ycomm,ierr)
IF (yrank /= 0) CALL MPI_RECV(fpi(ieq,5),neqy,MPI_DOUBLE_PRECISION,yprv,myytag,ycomm,istat,ierr)
! Quadrant 6
IF (yrank /= npy-spy) CALL MPI_SEND(fpo(ieq,6),neqy,MPI_DOUBLE_PRECISION,ynxt,ynxtag,ycomm,ierr)
IF (yrank /= 0) CALL MPI_RECV(fpi(ieq,6),neqy,MPI_DOUBLE_PRECISION,yprv,myytag,ycomm,istat,ierr)

! Perform the sends 'down y'
! Quadrant 3
IF (yrank /= 0) CALL MPI_SEND(fpo(ieq,3),neqy,MPI_DOUBLE_PRECISION,yprv,yprtag,ycomm,ierr)
IF (yrank /= npy-spy) CALL MPI_RECV(fpi(ieq,3),neqy,MPI_DOUBLE_PRECISION,ynxt,myytag,ycomm,istat,ierr)
! Quadrant 4
IF (yrank /= 0) CALL MPI_SEND(fpo(ieq,4),neqy,MPI_DOUBLE_PRECISION,yprv,yprtag,ycomm,ierr)
IF (yrank /= npy-spy) CALL MPI_RECV(fpi(ieq,4),neqy,MPI_DOUBLE_PRECISION,ynxt,myytag,ycomm,istat,ierr)
! Quadrant 7
IF (yrank /= 0) CALL MPI_SEND(fpo(ieq,7),neqy,MPI_DOUBLE_PRECISION,yprv,yprtag,ycomm,ierr)
IF (yrank /= npy-spy) CALL MPI_RECV(fpi(ieq,7),neqy,MPI_DOUBLE_PRECISION,ynxt,myytag,ycomm,istat,ierr)
! Quadrant 8
IF (yrank /= 0) CALL MPI_SEND(fpo(ieq,8),neqy,MPI_DOUBLE_PRECISION,yprv,yprtag,ycomm,ierr)
IF (yrank /= npy-spy) CALL MPI_RECV(fpi(ieq,8),neqy,MPI_DOUBLE_PRECISION,ynxt,myytag,ycomm,istat,ierr)

!--------------------------------------------------------------------------------------------------
! Perform the sends 'up x'
ieq = ieq + neqy
neqx = apo*yzs
! Quadrant 1
IF (xrank /= npx-spx) CALL MPI_SEND(fpo(ieq,1),neqx,MPI_DOUBLE_PRECISION,xnxt,xnxtag,xcomm,ierr)
IF (xrank /= 0) CALL MPI_RECV(fpi(ieq,1),neqx,MPI_DOUBLE_PRECISION,xprv,myxtag,xcomm,istat,ierr)
! Quadrant 4
IF (xrank /= npx-spx) CALL MPI_SEND(fpo(ieq,4),neqx,MPI_DOUBLE_PRECISION,xnxt,xnxtag,xcomm,ierr)
IF (xrank /= 0) CALL MPI_RECV(fpi(ieq,4),neqx,MPI_DOUBLE_PRECISION,xprv,myxtag,xcomm,istat,ierr)
! Quadrant 5
IF (xrank /= npx-spx) CALL MPI_SEND(fpo(ieq,5),neqx,MPI_DOUBLE_PRECISION,xnxt,xnxtag,xcomm,ierr)
IF (xrank /= 0) CALL MPI_RECV(fpi(ieq,5),neqx,MPI_DOUBLE_PRECISION,xprv,myxtag,xcomm,istat,ierr)
! Quadrant 8
IF (xrank /= npx-spx) CALL MPI_SEND(fpo(ieq,8),neqx,MPI_DOUBLE_PRECISION,xnxt,xnxtag,xcomm,ierr)
IF (xrank /= 0) CALL MPI_RECV(fpi(ieq,8),neqx,MPI_DOUBLE_PRECISION,xprv,myxtag,xcomm,istat,ierr)

! Perform the sends 'down x'
! Quadrant 2
IF (xrank /= 0) CALL MPI_SEND(fpo(ieq,2),neqx,MPI_DOUBLE_PRECISION,xprv,xprtag,xcomm,ierr)
IF (xrank /= npx-spx) CALL MPI_RECV(fpi(ieq,2),neqx,MPI_DOUBLE_PRECISION,xnxt,myxtag,xcomm,istat,ierr)
! Quadrant 3
IF (xrank /= 0) CALL MPI_SEND(fpo(ieq,3),neqx,MPI_DOUBLE_PRECISION,xprv,xprtag,xcomm,ierr)
IF (xrank /= npx-spx) CALL MPI_RECV(fpi(ieq,3),neqx,MPI_DOUBLE_PRECISION,xnxt,myxtag,xcomm,istat,ierr)
! Quadrant 6
IF (xrank /= 0) CALL MPI_SEND(fpo(ieq,6),neqx,MPI_DOUBLE_PRECISION,xprv,xprtag,xcomm,ierr)
IF (xrank /= npx-spx) CALL MPI_RECV(fpi(ieq,6),neqx,MPI_DOUBLE_PRECISION,xnxt,myxtag,xcomm,istat,ierr)
! Quadrant 7
IF (xrank /= 0) CALL MPI_SEND(fpo(ieq,7),neqx,MPI_DOUBLE_PRECISION,xprv,xprtag,xcomm,ierr)
IF (xrank /= npx-spx) CALL MPI_RECV(fpi(ieq,7),neqx,MPI_DOUBLE_PRECISION,xnxt,myxtag,xcomm,istat,ierr)

RETURN
END SUBROUTINE pbj
