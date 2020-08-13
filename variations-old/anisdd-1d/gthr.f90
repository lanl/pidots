SUBROUTINE gthr

!-------------------------------------------------------------
!
! Gathers all the f vectors and puts solution into flux 
! 
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER :: tag, indx, pi, eqs, i, g, ierr
INTEGER, DIMENSION(3) ::  istat
REAL*8, DIMENSION(:,:), ALLOCATABLE :: tmpf

INCLUDE 'mpif.h'

! Bring all the blocks' solutions for all the groups back to root for output
IF (irank /= root) THEN
   tag = 100 + irank
   CALL MPI_SEND(f,nx*sord*ng,MPI_DOUBLE_PRECISION,root,tag,MPI_COMM_WORLD,ierr)
ELSE
   indx = 0
   DO pi = 0, npx-1
      eqs = nxvec(pi)*sord*ng
      IF (pi /= 0) THEN
         ALLOCATE(tmpf(nxvec(pi)*sord,ng))
         tmpf = 0.0
         tag = 100 + pi
         CALL MPI_RECV(tmpf,eqs,MPI_DOUBLE_PRECISION,pi,tag,MPI_COMM_WORLD,istat,ierr)
         DO g = 1, ng
            DO i = 1, nxvec(pi)*sord
               flux(i+indx,g) = tmpf(i,g)
            END DO
         END DO
         DEALLOCATE(tmpf)
      ELSE IF (pi == 0) THEN
         ! Root copies its own solution into flux
         DO g = 1, ng
            DO i = 1, nx*sord
               flux(i+indx,g) = f(i,g)
            END DO
         END DO
      END IF
   indx = indx + nxvec(pi)*sord
   END DO
END IF

IF (irank == root) DEALLOCATE(nxvec)

RETURN
END SUBROUTINE gthr
