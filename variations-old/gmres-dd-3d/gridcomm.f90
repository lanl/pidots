SUBROUTINE gridcomm(ov,iv)

!-------------------------------------------------------------
!
!  Communicate vectors similar to psio/psii in PBJ code.
!  Communications done by surface, octant to nearest neighbors
!   to facilitate matrix vector operation in GMRES.
!
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER :: myztag, myytag, myxtag, znxt, zprv, ynxt, yprv, xnxt, xprv
INTEGER :: znxtag, zprtag, ynxtag, yprtag, xnxtag, xprtag, neqz, neqy, neqx, ieq, ierr
INTEGER, DIMENSION(3) :: istat
REAL*8, DIMENSION(bcs,8), INTENT(IN) :: ov
REAL*8, DIMENSION(bcs,8), INTENT(OUT) :: iv

INCLUDE 'mpif.h'

! Initialize vector
iv = 0.0

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

! Send ov to neighbors. Inner blocks have 24 communications (8 octants to
! 3 neigbors each). Then set the incoming data into iv

! Perform the sends 'up z'
neqz = apo*xys
! Octant 1
IF (zrank /= npz-1) CALL MPI_SEND(ov(1,1),neqz,MPI_DOUBLE_PRECISION,znxt,znxtag,zcomm,ierr)
IF (zrank /= 0) CALL MPI_RECV(iv(1,1),neqz,MPI_DOUBLE_PRECISION,zprv,myztag,zcomm,istat,ierr)
! Octant 2
IF (zrank /= npz-1) CALL MPI_SEND(ov(1,2),neqz,MPI_DOUBLE_PRECISION,znxt,znxtag,zcomm,ierr)
IF (zrank /= 0) CALL MPI_RECV(iv(1,2),neqz,MPI_DOUBLE_PRECISION,zprv,myztag,zcomm,istat,ierr)
! Octant 3
IF (zrank /= npz-1) CALL MPI_SEND(ov(1,3),neqz,MPI_DOUBLE_PRECISION,znxt,znxtag,zcomm,ierr)
IF (zrank /= 0) CALL MPI_RECV(iv(1,3),neqz,MPI_DOUBLE_PRECISION,zprv,myztag,zcomm,istat,ierr)
! Octant 4
IF (zrank /= npz-1) CALL MPI_SEND(ov(1,4),neqz,MPI_DOUBLE_PRECISION,znxt,znxtag,zcomm,ierr)
IF (zrank /= 0) CALL MPI_RECV(iv(1,4),neqz,MPI_DOUBLE_PRECISION,zprv,myztag,zcomm,istat,ierr)

! Perform the sends 'down z'
! Octant 5
IF (zrank /= 0) CALL MPI_SEND(ov(1,5),neqz,MPI_DOUBLE_PRECISION,zprv,zprtag,zcomm,ierr)
IF (zrank /= npz-1) CALL MPI_RECV(iv(1,5),neqz,MPI_DOUBLE_PRECISION,znxt,myztag,zcomm,istat,ierr)
! Octant 6
IF (zrank /= 0) CALL MPI_SEND(ov(1,6),neqz,MPI_DOUBLE_PRECISION,zprv,zprtag,zcomm,ierr)
IF (zrank /= npz-1) CALL MPI_RECV(iv(1,6),neqz,MPI_DOUBLE_PRECISION,znxt,myztag,zcomm,istat,ierr)
! Octant 7
IF (zrank /= 0) CALL MPI_SEND(ov(1,7),neqz,MPI_DOUBLE_PRECISION,zprv,zprtag,zcomm,ierr)
IF (zrank /= npz-1) CALL MPI_RECV(iv(1,7),neqz,MPI_DOUBLE_PRECISION,znxt,myztag,zcomm,istat,ierr)
! Octant 8
IF (zrank /= 0) CALL MPI_SEND(ov(1,8),neqz,MPI_DOUBLE_PRECISION,zprv,zprtag,zcomm,ierr)
IF (zrank /= npz-1) CALL MPI_RECV(iv(1,8),neqz,MPI_DOUBLE_PRECISION,znxt,myztag,zcomm,istat,ierr)

! Update for reflective boundary conditions at zrank = 0, npz-1
IF (zrank == 0 .AND. bc(1) == 1) THEN
   iv(1:neqz,1) = ov(1:neqz,5)
   iv(1:neqz,2) = ov(1:neqz,6)
   iv(1:neqz,3) = ov(1:neqz,7)
   iv(1:neqz,4) = ov(1:neqz,8)
END IF
IF (zrank == npz-1 .AND. bc(2) == 1) THEN
   iv(1:neqz,5) = ov(1:neqz,1)
   iv(1:neqz,6) = ov(1:neqz,2)
   iv(1:neqz,7) = ov(1:neqz,3)
   iv(1:neqz,8) = ov(1:neqz,4)
END IF

! Perform the sends 'up y'
ieq = neqz + 1
neqy = apo*xzs
! Quadrant 1
IF (yrank /= npy-1) CALL MPI_SEND(ov(ieq,1),neqy,MPI_DOUBLE_PRECISION,ynxt,ynxtag,ycomm,ierr)
IF (yrank /= 0) CALL MPI_RECV(iv(ieq,1),neqy,MPI_DOUBLE_PRECISION,yprv,myytag,ycomm,istat,ierr)
! Quadrant 2
IF (yrank /= npy-1) CALL MPI_SEND(ov(ieq,2),neqy,MPI_DOUBLE_PRECISION,ynxt,ynxtag,ycomm,ierr)
IF (yrank /= 0) CALL MPI_RECV(iv(ieq,2),neqy,MPI_DOUBLE_PRECISION,yprv,myytag,ycomm,istat,ierr)
! Quadrant 5
IF (yrank /= npy-1) CALL MPI_SEND(ov(ieq,5),neqy,MPI_DOUBLE_PRECISION,ynxt,ynxtag,ycomm,ierr)
IF (yrank /= 0) CALL MPI_RECV(iv(ieq,5),neqy,MPI_DOUBLE_PRECISION,yprv,myytag,ycomm,istat,ierr)
! Quadrant 6
IF (yrank /= npy-1) CALL MPI_SEND(ov(ieq,6),neqy,MPI_DOUBLE_PRECISION,ynxt,ynxtag,ycomm,ierr)
IF (yrank /= 0) CALL MPI_RECV(iv(ieq,6),neqy,MPI_DOUBLE_PRECISION,yprv,myytag,ycomm,istat,ierr)

! Perform the sends 'down y'
! Quadrant 3
IF (yrank /= 0) CALL MPI_SEND(ov(ieq,3),neqy,MPI_DOUBLE_PRECISION,yprv,yprtag,ycomm,ierr)
IF (yrank /= npy-1) CALL MPI_RECV(iv(ieq,3),neqy,MPI_DOUBLE_PRECISION,ynxt,myytag,ycomm,istat,ierr)
! Quadrant 4
IF (yrank /= 0) CALL MPI_SEND(ov(ieq,4),neqy,MPI_DOUBLE_PRECISION,yprv,yprtag,ycomm,ierr)
IF (yrank /= npy-1) CALL MPI_RECV(iv(ieq,4),neqy,MPI_DOUBLE_PRECISION,ynxt,myytag,ycomm,istat,ierr)
! Quadrant 7
IF (yrank /= 0) CALL MPI_SEND(ov(ieq,7),neqy,MPI_DOUBLE_PRECISION,yprv,yprtag,ycomm,ierr)
IF (yrank /= npy-1) CALL MPI_RECV(iv(ieq,7),neqy,MPI_DOUBLE_PRECISION,ynxt,myytag,ycomm,istat,ierr)
! Quadrant 8
IF (yrank /= 0) CALL MPI_SEND(ov(ieq,8),neqy,MPI_DOUBLE_PRECISION,yprv,yprtag,ycomm,ierr)
IF (yrank /= npy-1) CALL MPI_RECV(iv(ieq,8),neqy,MPI_DOUBLE_PRECISION,ynxt,myytag,ycomm,istat,ierr)

! Update for reflective boundary conditions at yrank = 0, npy-1
IF (yrank == 0 .AND. bc(3) == 1) THEN
   iv(ieq:(neqz+neqy),1) = ov(ieq:(neqz+neqy),4)
   iv(ieq:(neqz+neqy),2) = ov(ieq:(neqz+neqy),3)
   iv(ieq:(neqz+neqy),5) = ov(ieq:(neqz+neqy),8)
   iv(ieq:(neqz+neqy),6) = ov(ieq:(neqz+neqy),7)
END IF
IF (yrank == npy-1 .AND. bc(4) == 1) THEN
   iv(ieq:(neqz+neqy),4) = ov(ieq:(neqz+neqy),1)
   iv(ieq:(neqz+neqy),3) = ov(ieq:(neqz+neqy),2)
   iv(ieq:(neqz+neqy),8) = ov(ieq:(neqz+neqy),5)
   iv(ieq:(neqz+neqy),7) = ov(ieq:(neqz+neqy),6)
END IF

! Perform the sends 'up x'
ieq = ieq + neqy
neqx = apo*yzs
! Quadrant 1
IF (xrank /= npx-1) CALL MPI_SEND(ov(ieq,1),neqx,MPI_DOUBLE_PRECISION,xnxt,xnxtag,xcomm,ierr)
IF (xrank /= 0) CALL MPI_RECV(iv(ieq,1),neqx,MPI_DOUBLE_PRECISION,xprv,myxtag,xcomm,istat,ierr)
! Quadrant 4
IF (xrank /= npx-1) CALL MPI_SEND(ov(ieq,4),neqx,MPI_DOUBLE_PRECISION,xnxt,xnxtag,xcomm,ierr)
IF (xrank /= 0) CALL MPI_RECV(iv(ieq,4),neqx,MPI_DOUBLE_PRECISION,xprv,myxtag,xcomm,istat,ierr)
! Quadrant 5
IF (xrank /= npx-1) CALL MPI_SEND(ov(ieq,5),neqx,MPI_DOUBLE_PRECISION,xnxt,xnxtag,xcomm,ierr)
IF (xrank /= 0) CALL MPI_RECV(iv(ieq,5),neqx,MPI_DOUBLE_PRECISION,xprv,myxtag,xcomm,istat,ierr)
! Quadrant 8
IF (xrank /= npx-1) CALL MPI_SEND(ov(ieq,8),neqx,MPI_DOUBLE_PRECISION,xnxt,xnxtag,xcomm,ierr)
IF (xrank /= 0) CALL MPI_RECV(iv(ieq,8),neqx,MPI_DOUBLE_PRECISION,xprv,myxtag,xcomm,istat,ierr)

! Perform the sends 'down x'
! Quadrant 2
IF (xrank /= 0) CALL MPI_SEND(ov(ieq,2),neqx,MPI_DOUBLE_PRECISION,xprv,xprtag,xcomm,ierr)
IF (xrank /= npx-1) CALL MPI_RECV(iv(ieq,2),neqx,MPI_DOUBLE_PRECISION,xnxt,myxtag,xcomm,istat,ierr)
! Quadrant 3
IF (xrank /= 0) CALL MPI_SEND(ov(ieq,3),neqx,MPI_DOUBLE_PRECISION,xprv,xprtag,xcomm,ierr)
IF (xrank /= npx-1) CALL MPI_RECV(iv(ieq,3),neqx,MPI_DOUBLE_PRECISION,xnxt,myxtag,xcomm,istat,ierr)
! Quadrant 6
IF (xrank /= 0) CALL MPI_SEND(ov(ieq,6),neqx,MPI_DOUBLE_PRECISION,xprv,xprtag,xcomm,ierr)
IF (xrank /= npx-1) CALL MPI_RECV(iv(ieq,6),neqx,MPI_DOUBLE_PRECISION,xnxt,myxtag,xcomm,istat,ierr)
! Quadrant 7
IF (xrank /= 0) CALL MPI_SEND(ov(ieq,7),neqx,MPI_DOUBLE_PRECISION,xprv,xprtag,xcomm,ierr)
IF (xrank /= npx-1) CALL MPI_RECV(iv(ieq,7),neqx,MPI_DOUBLE_PRECISION,xnxt,myxtag,xcomm,istat,ierr)

! Update for reflective boundary conditions at xrank = 0, npx-1
IF (xrank == 0 .AND. bc(5) == 1) THEN
   iv(ieq:(neqz+neqy+neqx),1) = ov(ieq:(neqz+neqy+neqx),2)
   iv(ieq:(neqz+neqy+neqx),4) = ov(ieq:(neqz+neqy+neqx),3)
   iv(ieq:(neqz+neqy+neqx),5) = ov(ieq:(neqz+neqy+neqx),6)
   iv(ieq:(neqz+neqy+neqx),8) = ov(ieq:(neqz+neqy+neqx),7)
END IF
IF (xrank == npx-1 .AND. bc(6) == 1) THEN
   iv(ieq:(neqz+neqy+neqx),2) = ov(ieq:(neqz+neqy+neqx),1)
   iv(ieq:(neqz+neqy+neqx),3) = ov(ieq:(neqz+neqy+neqx),4)
   iv(ieq:(neqz+neqy+neqx),6) = ov(ieq:(neqz+neqy+neqx),5)
   iv(ieq:(neqz+neqy+neqx),7) = ov(ieq:(neqz+neqy+neqx),8)
END IF

RETURN
END SUBROUTINE gridcomm
