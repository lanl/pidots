SUBROUTINE pbj(bit,its,bcnvf)

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
USE timevar
IMPLICIT NONE
INTEGER, INTENT(IN) :: bit, its
INTEGER, INTENT(OUT) :: bcnvf
INTEGER :: i, j, k, itsm
INTEGER :: myztag, myytag, myxtag, znxt, zprv, ynxt, yprv, xnxt, xprv
INTEGER :: znxtag, zprtag, ynxtag, yprtag, xnxtag, xprtag, neqz, neqy, neqx, ieq, ierr
INTEGER, DIMENSION(3) :: istat
REAL*8 :: df, dfmx

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
   df = MAXVAL(ABS((phi - phiold)/phiold))
ELSE
   df = MAXVAL(ABS(phi - phiold))
END IF

IF (isize == 1) THEN
   IF (bc(1)/=1 .AND. bc(2)/=1 .AND. bc(3)/=1 .AND. bc(4)/=1 .AND. bc(5)/=1 .AND. bc(6)/=1) THEN
      df = 0.0
   END IF
END IF

! Get the overall max difference
CALL MPI_ALLREDUCE(df,dfmx,1,MPI_DOUBLE_PRECISION,MPI_MAX,allcomm,ierr)

! Only get the max number of SI/CG iterations from sub-domains if requested
IF (itp == 1) THEN
   CALL MPI_REDUCE(its,itsm,1,MPI_INTEGER,MPI_MAX,0,allcomm,ierr)
END IF

! Get copy of old psii for computing the residual
!IF (bit == bitmx) THEN
   oldpsi = psii
!END IF

!IF (dfmx > err .AND. bit < bitmx) THEN
   ! Only print iteration information if requested
   IF (irank == root) THEN
      IF (itp == 1) WRITE(8,111) bit, dfmx, itsm
   END IF
   ! Set the previous iterate of the flux equal to the current
   phiold = phi
   ! Send psio to neighbors. Inner blocks have 24 communications (8 octants to
   ! 3 neigbors each). Then set the incoming data into psii

   ! Perform the sends 'up z'
   neqz = apo*xys
   ! Octant 1
   IF (zrank /= npz-1) CALL MPI_SEND(psio(1,1),neqz,MPI_DOUBLE_PRECISION,znxt,znxtag,zcomm,ierr)
   IF (zrank /= 0) CALL MPI_RECV(psii(1,1),neqz,MPI_DOUBLE_PRECISION,zprv,myztag,zcomm,istat,ierr)
   ! Octant 2
   IF (zrank /= npz-1) CALL MPI_SEND(psio(1,2),neqz,MPI_DOUBLE_PRECISION,znxt,znxtag,zcomm,ierr)
   IF (zrank /= 0) CALL MPI_RECV(psii(1,2),neqz,MPI_DOUBLE_PRECISION,zprv,myztag,zcomm,istat,ierr)
   ! Octant 3
   IF (zrank /= npz-1) CALL MPI_SEND(psio(1,3),neqz,MPI_DOUBLE_PRECISION,znxt,znxtag,zcomm,ierr)
   IF (zrank /= 0) CALL MPI_RECV(psii(1,3),neqz,MPI_DOUBLE_PRECISION,zprv,myztag,zcomm,istat,ierr)
   ! Octant 4
   IF (zrank /= npz-1) CALL MPI_SEND(psio(1,4),neqz,MPI_DOUBLE_PRECISION,znxt,znxtag,zcomm,ierr)
   IF (zrank /= 0) CALL MPI_RECV(psii(1,4),neqz,MPI_DOUBLE_PRECISION,zprv,myztag,zcomm,istat,ierr)

   ! Perform the sends 'down z'
   ! Octant 5
   IF (zrank /= 0) CALL MPI_SEND(psio(1,5),neqz,MPI_DOUBLE_PRECISION,zprv,zprtag,zcomm,ierr)
   IF (zrank /= npz-1) CALL MPI_RECV(psii(1,5),neqz,MPI_DOUBLE_PRECISION,znxt,myztag,zcomm,istat,ierr)
   ! Octant 6
   IF (zrank /= 0) CALL MPI_SEND(psio(1,6),neqz,MPI_DOUBLE_PRECISION,zprv,zprtag,zcomm,ierr)
   IF (zrank /= npz-1) CALL MPI_RECV(psii(1,6),neqz,MPI_DOUBLE_PRECISION,znxt,myztag,zcomm,istat,ierr)
   ! Octant 7
   IF (zrank /= 0) CALL MPI_SEND(psio(1,7),neqz,MPI_DOUBLE_PRECISION,zprv,zprtag,zcomm,ierr)
   IF (zrank /= npz-1) CALL MPI_RECV(psii(1,7),neqz,MPI_DOUBLE_PRECISION,znxt,myztag,zcomm,istat,ierr)
   ! Octant 8
   IF (zrank /= 0) CALL MPI_SEND(psio(1,8),neqz,MPI_DOUBLE_PRECISION,zprv,zprtag,zcomm,ierr)
   IF (zrank /= npz-1) CALL MPI_RECV(psii(1,8),neqz,MPI_DOUBLE_PRECISION,znxt,myztag,zcomm,istat,ierr)

   ! Update for reflective boundary conditions at zrank = 0 and npz-1. Only
   ! apply new boundary conditions to those on the z-face for each incoming octant
   IF (zrank == 0 .AND. bc(1) == 1) THEN         ! Reflective at low-z
      psii(1:neqz,1) = psio(1:neqz,5)
      psii(1:neqz,2) = psio(1:neqz,6)
      psii(1:neqz,3) = psio(1:neqz,7)
      psii(1:neqz,4) = psio(1:neqz,8)
   END IF
   IF (zrank == npz-1 .AND. bc(2) == 1) THEN     ! Reflective at high-z
      psii(1:neqz,5) = psio(1:neqz,1)
      psii(1:neqz,6) = psio(1:neqz,2)
      psii(1:neqz,7) = psio(1:neqz,3)
      psii(1:neqz,8) = psio(1:neqz,4)
   END IF

   ! Perform the sends 'up y'
   ieq = neqz + 1
   neqy = apo*xzs
   ! Quadrant 1
   IF (yrank /= npy-1) CALL MPI_SEND(psio(ieq,1),neqy,MPI_DOUBLE_PRECISION,ynxt,ynxtag,ycomm,ierr)
   IF (yrank /= 0) CALL MPI_RECV(psii(ieq,1),neqy,MPI_DOUBLE_PRECISION,yprv,myytag,ycomm,istat,ierr)
   ! Quadrant 2
   IF (yrank /= npy-1) CALL MPI_SEND(psio(ieq,2),neqy,MPI_DOUBLE_PRECISION,ynxt,ynxtag,ycomm,ierr)
   IF (yrank /= 0) CALL MPI_RECV(psii(ieq,2),neqy,MPI_DOUBLE_PRECISION,yprv,myytag,ycomm,istat,ierr)
   ! Quadrant 5
   IF (yrank /= npy-1) CALL MPI_SEND(psio(ieq,5),neqy,MPI_DOUBLE_PRECISION,ynxt,ynxtag,ycomm,ierr)
   IF (yrank /= 0) CALL MPI_RECV(psii(ieq,5),neqy,MPI_DOUBLE_PRECISION,yprv,myytag,ycomm,istat,ierr)
   ! Quadrant 6
   IF (yrank /= npy-1) CALL MPI_SEND(psio(ieq,6),neqy,MPI_DOUBLE_PRECISION,ynxt,ynxtag,ycomm,ierr)
   IF (yrank /= 0) CALL MPI_RECV(psii(ieq,6),neqy,MPI_DOUBLE_PRECISION,yprv,myytag,ycomm,istat,ierr)

   ! Perform the sends 'down y'
   ! Quadrant 3
   IF (yrank /= 0) CALL MPI_SEND(psio(ieq,3),neqy,MPI_DOUBLE_PRECISION,yprv,yprtag,ycomm,ierr)
   IF (yrank /= npy-1) CALL MPI_RECV(psii(ieq,3),neqy,MPI_DOUBLE_PRECISION,ynxt,myytag,ycomm,istat,ierr)
   ! Quadrant 4
   IF (yrank /= 0) CALL MPI_SEND(psio(ieq,4),neqy,MPI_DOUBLE_PRECISION,yprv,yprtag,ycomm,ierr)
   IF (yrank /= npy-1) CALL MPI_RECV(psii(ieq,4),neqy,MPI_DOUBLE_PRECISION,ynxt,myytag,ycomm,istat,ierr)
   ! Quadrant 7
   IF (yrank /= 0) CALL MPI_SEND(psio(ieq,7),neqy,MPI_DOUBLE_PRECISION,yprv,yprtag,ycomm,ierr)
   IF (yrank /= npy-1) CALL MPI_RECV(psii(ieq,7),neqy,MPI_DOUBLE_PRECISION,ynxt,myytag,ycomm,istat,ierr)
   ! Quadrant 8
   IF (yrank /= 0) CALL MPI_SEND(psio(ieq,8),neqy,MPI_DOUBLE_PRECISION,yprv,yprtag,ycomm,ierr)
   IF (yrank /= npy-1) CALL MPI_RECV(psii(ieq,8),neqy,MPI_DOUBLE_PRECISION,ynxt,myytag,ycomm,istat,ierr)

   ! Update for reflective boundary conditions at yrank = 0 and npy-1. Only
   ! apply new boundary conditions to those on the y-face for each incoming octant
   IF (yrank == 0 .AND. bc(3) == 1) THEN         ! Reflective at low-y
      psii(ieq:(neqz+neqy),1) = psio(ieq:(neqz+neqy),4)
      psii(ieq:(neqz+neqy),2) = psio(ieq:(neqz+neqy),3)
      psii(ieq:(neqz+neqy),5) = psio(ieq:(neqz+neqy),8)
      psii(ieq:(neqz+neqy),6) = psio(ieq:(neqz+neqy),7)
   END IF
   IF (yrank == npy-1 .AND. bc(4) == 1) THEN     ! Reflective at high-y
      psii(ieq:(neqz+neqy),4) = psio(ieq:(neqz+neqy),1)
      psii(ieq:(neqz+neqy),3) = psio(ieq:(neqz+neqy),2)
      psii(ieq:(neqz+neqy),8) = psio(ieq:(neqz+neqy),5)
      psii(ieq:(neqz+neqy),7) = psio(ieq:(neqz+neqy),6)
   END IF

   ! Perform the sends 'up x'
   ieq = ieq + neqy
   neqx = apo*yzs
   ! Quadrant 1
   IF (xrank /= npx-1) CALL MPI_SEND(psio(ieq,1),neqx,MPI_DOUBLE_PRECISION,xnxt,xnxtag,xcomm,ierr)
   IF (xrank /= 0) CALL MPI_RECV(psii(ieq,1),neqx,MPI_DOUBLE_PRECISION,xprv,myxtag,xcomm,istat,ierr)
   ! Quadrant 4
   IF (xrank /= npx-1) CALL MPI_SEND(psio(ieq,4),neqx,MPI_DOUBLE_PRECISION,xnxt,xnxtag,xcomm,ierr)
   IF (xrank /= 0) CALL MPI_RECV(psii(ieq,4),neqx,MPI_DOUBLE_PRECISION,xprv,myxtag,xcomm,istat,ierr)
   ! Quadrant 5
   IF (xrank /= npx-1) CALL MPI_SEND(psio(ieq,5),neqx,MPI_DOUBLE_PRECISION,xnxt,xnxtag,xcomm,ierr)
   IF (xrank /= 0) CALL MPI_RECV(psii(ieq,5),neqx,MPI_DOUBLE_PRECISION,xprv,myxtag,xcomm,istat,ierr)
   ! Quadrant 8
   IF (xrank /= npx-1) CALL MPI_SEND(psio(ieq,8),neqx,MPI_DOUBLE_PRECISION,xnxt,xnxtag,xcomm,ierr)
   IF (xrank /= 0) CALL MPI_RECV(psii(ieq,8),neqx,MPI_DOUBLE_PRECISION,xprv,myxtag,xcomm,istat,ierr)

   ! Perform the sends 'down x'
   ! Quadrant 2
   IF (xrank /= 0) CALL MPI_SEND(psio(ieq,2),neqx,MPI_DOUBLE_PRECISION,xprv,xprtag,xcomm,ierr)
   IF (xrank /= npx-1) CALL MPI_RECV(psii(ieq,2),neqx,MPI_DOUBLE_PRECISION,xnxt,myxtag,xcomm,istat,ierr)
   ! Quadrant 3
   IF (xrank /= 0) CALL MPI_SEND(psio(ieq,3),neqx,MPI_DOUBLE_PRECISION,xprv,xprtag,xcomm,ierr)
   IF (xrank /= npx-1) CALL MPI_RECV(psii(ieq,3),neqx,MPI_DOUBLE_PRECISION,xnxt,myxtag,xcomm,istat,ierr)
   ! Quadrant 6
   IF (xrank /= 0) CALL MPI_SEND(psio(ieq,6),neqx,MPI_DOUBLE_PRECISION,xprv,xprtag,xcomm,ierr)
   IF (xrank /= npx-1) CALL MPI_RECV(psii(ieq,6),neqx,MPI_DOUBLE_PRECISION,xnxt,myxtag,xcomm,istat,ierr)
   ! Quadrant 7
   IF (xrank /= 0) CALL MPI_SEND(psio(ieq,7),neqx,MPI_DOUBLE_PRECISION,xprv,xprtag,xcomm,ierr)
   IF (xrank /= npx-1) CALL MPI_RECV(psii(ieq,7),neqx,MPI_DOUBLE_PRECISION,xnxt,myxtag,xcomm,istat,ierr)

   ! Update for reflective boundary conditions at xrank = 0 and npx-1. Only
   ! apply new boundary conditions to those on the x-face for each incoming octant
   IF (xrank == 0 .AND. bc(5) == 1) THEN         ! Reflective at low-x
      psii(ieq:(neqz+neqy+neqx),1) = psio(ieq:(neqz+neqy+neqx),2)
      psii(ieq:(neqz+neqy+neqx),4) = psio(ieq:(neqz+neqy+neqx),3)
      psii(ieq:(neqz+neqy+neqx),5) = psio(ieq:(neqz+neqy+neqx),6)
      psii(ieq:(neqz+neqy+neqx),8) = psio(ieq:(neqz+neqy+neqx),7)
   END IF
   IF (xrank == npx-1 .AND. bc(6) == 1) THEN     ! Reflective at high-x
      psii(ieq:(neqz+neqy+neqx),2) = psio(ieq:(neqz+neqy+neqx),1)
      psii(ieq:(neqz+neqy+neqx),3) = psio(ieq:(neqz+neqy+neqx),4)
      psii(ieq:(neqz+neqy+neqx),6) = psio(ieq:(neqz+neqy+neqx),5)
      psii(ieq:(neqz+neqy+neqx),7) = psio(ieq:(neqz+neqy+neqx),8)
   END IF

!ELSE IF (dfmx < err) THEN
IF (dfmx < err) THEN
   bcnvf = 1
!   CALL CPU_TIME(tsolve)
!   IF (irank == root) THEN
!      WRITE(8,'(2X,A,I5,A,ES11.3,/)') "System converged in ", bit, " iterations. Max error = ", dfmx
!      IF (itp == 1 .AND. isize == 1) THEN
!         IF (meth == 1 .AND. idos == 1) WRITE(8,'(2X,A,I5,A)') "Serial ITM-CG run used ", its, " iterations"
!         IF (meth == 1 .AND. idos == 2) WRITE(8,'(2X,A,I5,A)') "Serial ITM-SOR run used ", its, " iterations"
!      END IF
!   END IF
ELSE IF (bit == bitmx) THEN
   bcnvf = 0
!   CALL CPU_TIME(tsolve)
!   IF (irank == root) THEN
!      WRITE (8,'(/,2X,A,ES11.3)') "ERROR: Did not converge. Max error= ", dfmx
!      WRITE (8,*) "Will adversely affect other results. Aborting program."
!   END IF
END IF

111 FORMAT(2X,'It ', I5,' Dfmx ', ES11.3,' CG Its ', I5)

RETURN
END SUBROUTINE pbj
