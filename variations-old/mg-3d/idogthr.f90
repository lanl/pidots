SUBROUTINE idogthr

!-------------------------------------------------------------
!
! Gathers all the phi vectors and puts solution into flux 
! 
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER :: root2, tag, kndx, jndx, indx, pk, pj, pi, trank, ieq
INTEGER :: i, j, k, ierr
INTEGER, DIMENSION(3) :: coord, istat
REAL*8, DIMENSION(neq) :: tmphi

INCLUDE 'mpif.h'

! Setup for parallel
IF (irank == root) root2 = nrank
CALL MPI_BCAST(root2,1,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)

! Bring all the blocks' solutions for all the groups back to root for output
IF (nrank /= root2) THEN
   tag = 100 + nrank
   CALL MPI_SEND(phi(:,1),neq,MPI_DOUBLE_PRECISION,root2,tag,allcomm,ierr)
ELSE
   DO pk = 0, npz-1
      coord(1) = pk
      kndx = pk*nz
      DO pj = 0, npy-1
         coord(2) = pj
         jndx = pj*ny
         DO pi = 0, npx-1
            coord(3) = pi
            indx = pi*nx
            CALL MPI_CART_RANK(allcomm,coord,trank,ierr)
            IF (trank /= nrank) THEN
               tmphi = 0.0
               tag = 100 + trank
               CALL MPI_RECV(tmphi,neq,MPI_DOUBLE_PRECISION,trank,tag,allcomm,istat,ierr)
               DO k = 1, nz
                  DO j = 1, ny
                     DO i = 1, nx
                        ieq = i + (j-1)*nx + (k-1)*xys
                        flux(i+indx,j+jndx,k+kndx) = tmphi(ieq)
                     END DO
                  END DO
               END DO
            ELSE IF (trank == nrank) THEN
               ! Root copies its own solution into flux
               DO k = 1, nz
                  DO j = 1, ny
                     DO i = 1, nx
                        ieq = i + (j-1)*nx + (k-1)*xys
                        flux(i+indx,j+jndx,k+kndx) = phi(ieq,1)
                     END DO
                  END DO
               END DO
            END IF
         END DO
      END DO
   END DO
END IF

RETURN
END SUBROUTINE idogthr
