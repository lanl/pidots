SUBROUTINE pbj(g,nsm,bit,its)

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
INTEGER, INTENT(IN) :: g, nsm, bit, its
INTEGER :: i, itsm, tag, next, prv, ntag, ptag, ierr
INTEGER, DIMENSION(3) :: istat
REAL*8 :: df, dfmx

INCLUDE 'mpif.h'

! Set up tags and ranks
tag = 100 + irank
next = irank + 1
prv = irank - 1
ntag = 100 + next
ptag = 100 + prv

! Check convergence of scalar flux depending on the method
! Each process checks its own fluxes, then perform MPI Reduction
dfmx = -1.0
DO i = 1, nsm, sord
   IF (ABS(fold(i)) >= tolr) THEN
      df = ABS((f(i,g)-fold(i))/fold(i))
   ELSE
      df = ABS(f(i,g) - fold(i))
   END IF
   IF (df > dfmx) dfmx = df
END DO

IF (isize == 1) THEN
   IF (bc(1)/=1 .AND. bc(2)/=1) THEN
      df = 0.0
   END IF
END IF

! Get the overall max difference
CALL MPI_ALLREDUCE(df,dfmx,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)

! Only get the max number of SI/CG iterations from sub-domains if requested
IF (itp == 1) THEN
   CALL MPI_REDUCE(its,itsm,1,MPI_INTEGER,MPI_MAX,root,MPI_COMM_WORLD,ierr)
END IF

IF (dfmx > err .AND. bit < bitmx) THEN
   ! Only print iteration information if requested
   IF (irank == root) THEN
      IF (itp == 1) WRITE(8,111) g, bit, dfmx, itsm
   END IF
   ! Set the previous iterate of the flux equal to the current
   fold = f(:,g)

   ! Send psio to neighbors. Then set incoming data into psii

   ! Perform the sends 'up x'
   ! Direction 1
   IF (irank /= npx-1) CALL MPI_SEND(psio(1,1,g),apo,MPI_DOUBLE_PRECISION,next,ntag,MPI_COMM_WORLD,ierr)
   IF (irank /= 0) CALL MPI_RECV(psii(1,1),apo,MPI_DOUBLE_PRECISION,prv,tag,MPI_COMM_WORLD,istat,ierr)

   ! Perform the sends 'down x'
   ! Direction 2
   IF (irank /= 0) CALL MPI_SEND(psio(1,2,g),apo,MPI_DOUBLE_PRECISION,prv,ptag,MPI_COMM_WORLD,ierr)
   IF (irank /= npx-1) CALL MPI_RECV(psii(1,2),apo,MPI_DOUBLE_PRECISION,next,tag,MPI_COMM_WORLD,istat,ierr)

   ! Update for reflective boundary conditions at irank = 0 and npx-1.
   IF (irank == 0 .AND. bc(1) == 1) THEN         ! Reflective at low-x
      psii(:,1) = psio(:,2,g)
   END IF
   IF (irank == npx-1 .AND. bc(2) == 1) THEN     ! Reflective at high-x
      psii(:,2) = psio(:,1,g)
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
