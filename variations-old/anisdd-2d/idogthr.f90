SUBROUTINE idogthr(neq)

!-------------------------------------------------------------
!
! Gathers all the phi vectors and puts solution into flux 
! 
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER, INTENT(IN) :: neq
INTEGER :: root2, tag, jndx, indx, pj, pi, eqs, trank, ieq
INTEGER :: i, j, g, ierr
INTEGER, DIMENSION(2) :: coord
INTEGER, DIMENSION(3) :: istat
REAL*8, DIMENSION(:,:), ALLOCATABLE :: tmphi

INCLUDE 'mpif.h'

! Setup for parallel
IF (irank == root) root2 = nrank
CALL MPI_BCAST(root2,1,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)

! Bring all the blocks' solutions for all the groups back to root for output
IF (nrank /= root2) THEN
   tag = 100 + nrank
   CALL MPI_SEND(phi,neq*ng,MPI_DOUBLE_PRECISION,root2,tag,allcomm,ierr)
ELSE
   jndx = 0
   DO pj = 0, npy-1
      coord(2) = pj
      indx = 0
      DO pi = 0, npx-1
         coord(1) = pi
         eqs = nxvec(pi)*nyvec(pj)*nmom
         CALL MPI_CART_RANK(allcomm,coord,trank,ierr)
         IF (trank /= nrank) THEN
            ALLOCATE(tmphi(eqs,ng))
            tmphi = 0.0
            tag = 100 + trank
            CALL MPI_RECV(tmphi,eqs*ng,MPI_DOUBLE_PRECISION,trank,tag,allcomm,istat,ierr)
            DO g = 1, ng
               DO j = 1, nyvec(pj)
                  DO i = 1, nxvec(pi)
                     ieq = ((i-1) + (j-1)*nxvec(pi))*nmom + 1
                     flux(i+indx,j+jndx,g) = tmphi(ieq,g)
                  END DO
               END DO
            END DO
            DEALLOCATE(tmphi)
         ELSE IF (trank == nrank) THEN
            ! Root copies its own solution into flux
            DO g = 1, ng
               DO j = 1, ny
                  DO i = 1, nx
                     ieq = ((i-1) + (j-1)*nx)*nmom + 1
                     flux(i+indx,j+jndx,g) = phi(ieq,g)
                  END DO
               END DO
            END DO
         END IF
         indx = indx + nxvec(pi)
      END DO
      jndx = jndx + nyvec(pj)
   END DO
END IF

RETURN
END SUBROUTINE idogthr
