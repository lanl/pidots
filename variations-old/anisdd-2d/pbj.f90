SUBROUTINE pbj(g,neq,bit,its)

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
INTEGER, INTENT(IN) :: g, neq, bit, its
INTEGER :: i, j, k, itsm
INTEGER :: myztag, myytag, myxtag, znxt, zprv, ynxt, yprv, xnxt, xprv
INTEGER :: znxtag, zprtag, ynxtag, yprtag, xnxtag, xprtag, neqz, neqy, neqx, ieq, ierr
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
dfmx = -1.0
IF (meth == 0) THEN
   IF (MINVAL(ABS(fold(1,:,:))) >= tolr) THEN
      dfmx = MAXVAL(ABS((f(1,:,:,g)-fold(1,:,:))/fold(1,:,:)))
   ELSE
      dfmx = MAXVAL(ABS(f(1,:,:,g) - fold(1,:,:)))
   END IF
ELSE
   DO i = 1, neq, nmom
      IF (ABS(phiold(i)) >= tolr) THEN
         df = ABS((phi(i,g) - phiold(i))/phiold(i))
      ELSE
         df = ABS(phi(i,g) - phiold(i))
      END IF
   END DO
   IF (df > dfmx) dfmx = df
END IF

IF (isize == 1) THEN
   IF (bc(1)/=1 .AND. bc(2)/=1 .AND. bc(3)/=1 .AND. bc(4)/=1) THEN
      dfmx = 0.0
   END IF
END IF

df = dfmx

! Get the overall max difference
CALL MPI_ALLREDUCE(df,dfmx,1,MPI_DOUBLE_PRECISION,MPI_MAX,allcomm,ierr)

! Only get the max number of SI/CG iterations from sub-domains if requested
IF (itp == 1) THEN
   CALL MPI_REDUCE(its,itsm,1,MPI_INTEGER,MPI_MAX,0,allcomm,ierr)
END IF

IF (dfmx > err .AND. bit < bitmx) THEN
   ! Only print iteration information if requested
   IF (irank == root) THEN
      IF (itp == 1) WRITE(8,111) g, bit, dfmx, itsm
   END IF
   ! Set the previous iterate of the flux equal to the current
   IF (meth == 0) THEN
      fold = f(:,:,:,g)
   ELSE
      phiold = phi(:,g)
   END IF
   ! Send psio to neighbors. Inner blocks have 12 communications (4 octants to
   ! 2 neigbors each). Then set the incoming data into psii

   ! Perform the sends 'up y'
   neqy = apo*nx
   ! Quadrant 1
   IF (yrank /= npy-1) CALL MPI_SEND(psio(1,1,g),neqy,MPI_DOUBLE_PRECISION,ynxt,ynxtag,ycomm,ierr)
   IF (yrank /= 0) CALL MPI_RECV(psii(1,1),neqy,MPI_DOUBLE_PRECISION,yprv,myytag,ycomm,istat,ierr)
   ! Quadrant 2
   IF (yrank /= npy-1) CALL MPI_SEND(psio(1,2,g),neqy,MPI_DOUBLE_PRECISION,ynxt,ynxtag,ycomm,ierr)
   IF (yrank /= 0) CALL MPI_RECV(psii(1,2),neqy,MPI_DOUBLE_PRECISION,yprv,myytag,ycomm,istat,ierr)

   ! Perform the sends 'down y'
   ! Quadrant 3
   IF (yrank /= 0) CALL MPI_SEND(psio(1,3,g),neqy,MPI_DOUBLE_PRECISION,yprv,yprtag,ycomm,ierr)
   IF (yrank /= npy-1) CALL MPI_RECV(psii(1,3),neqy,MPI_DOUBLE_PRECISION,ynxt,myytag,ycomm,istat,ierr)
   ! Quadrant 4
   IF (yrank /= 0) CALL MPI_SEND(psio(1,4,g),neqy,MPI_DOUBLE_PRECISION,yprv,yprtag,ycomm,ierr)
   IF (yrank /= npy-1) CALL MPI_RECV(psii(1,4),neqy,MPI_DOUBLE_PRECISION,ynxt,myytag,ycomm,istat,ierr)

   ! Update for reflective boundary conditions at yrank = 0 and npy-1. Only
   ! apply new boundary conditions to those on the y-face for each incoming octant
   IF (yrank == 0 .AND. bc(1) == 1) THEN         ! Reflective at low-y
      psii(1:neqy,1) = psio(1:neqy,4,g)
      psii(1:neqy,2) = psio(1:neqy,3,g)
   END IF
   IF (yrank == npy-1 .AND. bc(2) == 1) THEN     ! Reflective at high-y
      psii(1:neqy,4) = psio(1:neqy,1,g)
      psii(1:neqy,3) = psio(1:neqy,2,g)
   END IF

   ! Perform the sends 'up x'
   ieq = 1 + neqy
   neqx = apo*ny
   ! Quadrant 1
   IF (xrank /= npx-1) CALL MPI_SEND(psio(ieq,1,g),neqx,MPI_DOUBLE_PRECISION,xnxt,xnxtag,xcomm,ierr)
   IF (xrank /= 0) CALL MPI_RECV(psii(ieq,1),neqx,MPI_DOUBLE_PRECISION,xprv,myxtag,xcomm,istat,ierr)
   ! Quadrant 4
   IF (xrank /= npx-1) CALL MPI_SEND(psio(ieq,4,g),neqx,MPI_DOUBLE_PRECISION,xnxt,xnxtag,xcomm,ierr)
   IF (xrank /= 0) CALL MPI_RECV(psii(ieq,4),neqx,MPI_DOUBLE_PRECISION,xprv,myxtag,xcomm,istat,ierr)

   ! Perform the sends 'down x'
   ! Quadrant 2
   IF (xrank /= 0) CALL MPI_SEND(psio(ieq,2,g),neqx,MPI_DOUBLE_PRECISION,xprv,xprtag,xcomm,ierr)
   IF (xrank /= npx-1) CALL MPI_RECV(psii(ieq,2),neqx,MPI_DOUBLE_PRECISION,xnxt,myxtag,xcomm,istat,ierr)
   ! Quadrant 3
   IF (xrank /= 0) CALL MPI_SEND(psio(ieq,3,g),neqx,MPI_DOUBLE_PRECISION,xprv,xprtag,xcomm,ierr)
   IF (xrank /= npx-1) CALL MPI_RECV(psii(ieq,3),neqx,MPI_DOUBLE_PRECISION,xnxt,myxtag,xcomm,istat,ierr)

   ! Update for reflective boundary conditions at xrank = 0 and npx-1. Only
   ! apply new boundary conditions to those on the x-face for each incoming octant
   IF (xrank == 0 .AND. bc(3) == 1) THEN         ! Reflective at low-x
      psii(ieq:(neqy+neqx),1) = psio(ieq:(neqy+neqx),2,g)
      psii(ieq:(neqy+neqx),4) = psio(ieq:(neqy+neqx),3,g)
   END IF
   IF (xrank == npx-1 .AND. bc(4) == 1) THEN     ! Reflective at high-x
      psii(ieq:(neqy+neqx),2) = psio(ieq:(neqy+neqx),1,g)
      psii(ieq:(neqy+neqx),3) = psio(ieq:(neqy+neqx),4,g)
   END IF

ELSE IF (dfmx < err) THEN
   bcnvf(g) = 1
   CALL CPU_TIME(tsolve)
   IF (irank == root) THEN
      WRITE(8,'(2X,A,I3,A,I5,A,ES11.3,/)') "Group ", g, " converged in ", bit, " iterations. Max error = ", dfmx
      IF (itp == 1 .AND. isize == 1) THEN
         IF (meth == 0) WRITE(8,'(2X,A,I5,A)') "Serial SI run used ", its, " iterations"
         IF (meth == 1 .AND. idos == 1) WRITE(8,'(2X,A,I5,A)') "Serial ITM-CG run used ", its, " iterations"
      END IF
  END IF
ELSE IF (bit == bitmx) THEN
   bcnvf(g) = 0
   CALL CPU_TIME(tsolve)
   IF (irank == root) THEN
      WRITE (8,'(/,2X,A,I2,A,ES11.3)') "ERROR: Group ", g, " did not converge. Max error= ", dfmx
      WRITE (8,*) "Will adversely affect other results. Aborting program."
   END IF
END IF

111 FORMAT(2X,'Gr',I3,' It ', I5,' Dfmx ', ES11.3,' SI/CG Its ', I5)

RETURN
END SUBROUTINE pbj
