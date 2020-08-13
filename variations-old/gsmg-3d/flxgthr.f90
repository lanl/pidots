SUBROUTINE flxgthr

!-------------------------------------------------------------
!
! Gathers all the flux arrays and puts solution into fluxt 
! 
!-------------------------------------------------------------

USE invar
USE totvar
IMPLICIT NONE
INTEGER :: root2, tag, kndx, jndx, indx, pk, pj, pi, eqs, trank
INTEGER :: i, j, k, ierr
INTEGER, DIMENSION(3) :: coord, istat
REAL*8, DIMENSION(nx,ny,nz) :: tmpf

INCLUDE 'mpif.h'

! Setup for parallel
IF (irank == root) root2 = nrank
CALL MPI_BCAST(root2,1,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)

eqs = nx*ny*nz

! Bring all the blocks' solutions for all the groups back to root for output
IF (nrank /= root2) THEN
   tag = 100 + nrank
   CALL MPI_SEND(flux,eqs,MPI_DOUBLE_PRECISION,root2,tag,allcomm,ierr)
ELSE
   kndx = 0
   DO pk = 0, npz-1
      coord(1) = pk
      jndx = 0
      DO pj = 0, npy-1
         coord(2) = pj
         indx = 0
         DO pi = 0, npx-1
            coord(3) = pi
            CALL MPI_CART_RANK(allcomm,coord,trank,ierr)
            IF (trank /= nrank) THEN
               tmpf = 0.0
               tag = 100 + trank
               CALL MPI_RECV(tmpf,eqs,MPI_DOUBLE_PRECISION,trank,tag,allcomm,istat,ierr)
               DO k = 1, nz
                  DO j = 1, ny
                     DO i = 1, nx
                        fluxt(i+indx,j+jndx,k+kndx) = tmpf(i,j,k)
                     END DO
                  END DO
               END DO
            ELSE IF (trank == nrank) THEN
               ! Root copies its own solution into flux
               DO k = 1, nz
                  DO j = 1, ny
                     DO i = 1, nx
                        fluxt(i+indx,j+jndx,k+kndx) = flux(i,j,k)
                     END DO
                  END DO
               END DO
            END IF
            indx = indx + nx
         END DO
         jndx = jndx + ny
      END DO
      kndx = kndx + nz
   END DO
END IF

RETURN
END SUBROUTINE flxgthr
