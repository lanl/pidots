SUBROUTINE pbj(g,bit,root)

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
INTEGER, INTENT(IN) :: g, bit, root
INTEGER :: myytag, myxtag, ynxt, yprv, xnxt, xprv
INTEGER :: ynxtag, yprtag, xnxtag, xprtag, neqy, neqx, ieq, ierr
INTEGER, DIMENSION(3) :: istat
REAL*8 :: df, dfmx

INCLUDE 'mpif.h'

! Set up some tags
myytag = 100 + yrank
myxtag = 100 + xrank
ynxt = yrank + 1
ynxtag = 100 + ynxt
yprv = yrank - 1
yprtag = 100 + yprv
xnxt = xrank + 1
xnxtag = 100 + xnxt
xprv = xrank - 1
xprtag = 100 + xprv

! Check convergence of scalar flux depending on the method
! Each process checks its own fluxes, then perform MPI Reduction
IF (meth == 0) THEN
   IF (MINVAL(ABS(fold)) >= tolr) THEN
      df = MAXVAL(ABS((f(:,:,:,:,g) - fold)/fold))
   ELSE
      df = MAXVAL(ABS(f(:,:,:,:,g) - fold))
   END IF
ELSE
   IF (MINVAL(ABS(phiold)) >= tolr) THEN
      df = MAXVAL(ABS((phi(:,g) - phiold)/phiold))
   ELSE
      df = MAXVAL(ABS(phi(:,g) - phiold))
   END IF
END IF

IF (isize == 1) df = 0.0

! Get the overall max difference
CALL MPI_ALLREDUCE(df,dfmx,1,MPI_DOUBLE_PRECISION,MPI_MAX,allcomm,ierr)

IF (dfmx > err .AND. bit < bitmx) THEN
   IF (irank == root) WRITE(8,111) g, bit, dfmx
   ! Set the previous iterate of the flux equal to the current
   IF (meth == 0) THEN
      fold = f(:,:,:,:,g)
   ELSE
      phiold = phi(:,g)
   END IF
   ! Send psio to neighbors. Inner blocks have 8 communications (4 quadrants to
   ! 2 neigbors each). Then set the incoming data into psi#

   ! Perform the sends 'up y'
   neqy = apo*nx*order
   ! Quadrant 1
   IF (yrank /= npy-1) CALL MPI_SEND(psio(1,1,g),neqy,MPI_DOUBLE_PRECISION,ynxt,ynxtag,ycomm,ierr)
   IF (yrank /= 0) CALL MPI_RECV(psi1(1),neqy,MPI_DOUBLE_PRECISION,yprv,myytag,ycomm,istat,ierr)
   ! Quadrant 2
   IF (yrank /= npy-1) CALL MPI_SEND(psio(1,2,g),neqy,MPI_DOUBLE_PRECISION,ynxt,ynxtag,ycomm,ierr)
   IF (yrank /= 0) CALL MPI_RECV(psi2(1),neqy,MPI_DOUBLE_PRECISION,yprv,myytag,ycomm,istat,ierr)

   ! Perform the sends 'down y'
   ! Quadrant 3
   IF (yrank /= 0) CALL MPI_SEND(psio(1,3,g),neqy,MPI_DOUBLE_PRECISION,yprv,yprtag,ycomm,ierr)
   IF (yrank /= npy-1) CALL MPI_RECV(psi3(1),neqy,MPI_DOUBLE_PRECISION,ynxt,myytag,ycomm,istat,ierr)
   ! Quadrant 4
   IF (yrank /= 0) CALL MPI_SEND(psio(1,4,g),neqy,MPI_DOUBLE_PRECISION,yprv,yprtag,ycomm,ierr)
   IF (yrank /= npy-1) CALL MPI_RECV(psi4(1),neqy,MPI_DOUBLE_PRECISION,ynxt,myytag,ycomm,istat,ierr)

   ieq = apo*nx*order + 1
   neqx = apo*ny*order
   ! Perform the sends 'up x'
   ! Quadrant 1
   IF (xrank /= npx-1) CALL MPI_SEND(psio(ieq,1,g),neqx,MPI_DOUBLE_PRECISION,xnxt,xnxtag,xcomm,ierr)
   IF (xrank /= 0) CALL MPI_RECV(psi1(ieq),neqx,MPI_DOUBLE_PRECISION,xprv,myxtag,xcomm,istat,ierr)
   ! Quadrant 4
   IF (xrank /= npx-1) CALL MPI_SEND(psio(ieq,4,g),neqx,MPI_DOUBLE_PRECISION,xnxt,xnxtag,xcomm,ierr)
   IF (xrank /= 0) CALL MPI_RECV(psi4(ieq),neqx,MPI_DOUBLE_PRECISION,xprv,myxtag,xcomm,istat,ierr)

   ! Perform the sends 'down x'
   ! Quadrant 2
   IF (xrank /= 0) CALL MPI_SEND(psio(ieq,2,g),neqx,MPI_DOUBLE_PRECISION,xprv,xprtag,xcomm,ierr)
   IF (xrank /= npx-1) CALL MPI_RECV(psi2(ieq),neqx,MPI_DOUBLE_PRECISION,xnxt,myxtag,xcomm,istat,ierr)
   ! Quadrant 3
   IF (xrank /= 0) CALL MPI_SEND(psio(ieq,3,g),neqx,MPI_DOUBLE_PRECISION,xprv,xprtag,xcomm,ierr)
   IF (xrank /= npx-1) CALL MPI_RECV(psi3(ieq),neqx,MPI_DOUBLE_PRECISION,xnxt,myxtag,xcomm,istat,ierr)

ELSE IF (dfmx < err) THEN
   bcnvf(g) = 1
   IF (irank == root) THEN
      WRITE(8,'(2X,A,I3,A,I5,A,ES11.3,/)') "Group ", g, " converged in ", bit, " iterations. Max error = ", dfmx
   END IF
ELSE IF (bit == bitmx) THEN
   bcnvf(g) = 0
   IF (irank == root) THEN
      WRITE (8,'(/,2X,A,I2,A,ES11.3)') "ERROR: Group ", g, " did not converge. Max error= ", dfmx
      WRITE (8,*) "Will adversely affect other results. Aborting program."
   END IF
END IF

111 FORMAT(2X,'Gr',I3,' It ', I5,' Dfmx ', ES11.3)

RETURN
END SUBROUTINE pbj
