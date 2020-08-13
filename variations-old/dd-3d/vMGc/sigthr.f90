SUBROUTINE sigthr

!-------------------------------------------------------------
!
! Gathers all the f vectors and puts solution into flux 
! 
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER :: root2, tag, kndx, jndx, indx, pk, pj, pi, eqs, trank
INTEGER :: i, j, k, g, ierr
INTEGER, DIMENSION(3) :: coord, istat
REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE :: tmpf

INCLUDE 'mpif.h'

! Setup for parallel
IF (irank == root) root2 = nrank
CALL MPI_BCAST(root2,1,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)

! Bring all the blocks' solutions for all the groups back to root for output
IF (nrank /= root2) THEN
   tag = 100 + nrank
   CALL MPI_SEND(f,neq*ng,MPI_DOUBLE_PRECISION,root2,tag,allcomm,ierr)
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
            eqs = nxvec(pi)*nyvec(pj)*nzvec(pk)*ng
            CALL MPI_CART_RANK(allcomm,coord,trank,ierr)
            IF (trank /= nrank) THEN
               ALLOCATE(tmpf(nxvec(pi),nyvec(pj),nzvec(pk),ng))
               tmpf = 0.0
               tag = 100 + trank
               CALL MPI_RECV(tmpf,eqs,MPI_DOUBLE_PRECISION,trank,tag,allcomm,istat,ierr)
               DO g = 1, ng
                  DO k = 1, nzvec(pk)
                     DO j = 1, nyvec(pj)
                        DO i = 1, nxvec(pi)
                           flux(i+indx,j+jndx,k+kndx,g) = tmpf(i,j,k,g)
                        END DO
                     END DO
                  END DO
               END DO
               DEALLOCATE(tmpf)
            ELSE IF (trank == nrank) THEN
               ! Root copies its own solution into flux
               DO g = 1, ng
                  DO k = 1, nz
                     DO j = 1, ny
                        DO i = 1, nx
                           flux(i+indx,j+jndx,k+kndx,g) = f(i,j,k,g)
                        END DO
                     END DO
                  END DO
               END DO
            END IF
            indx = indx + nxvec(pi)
         END DO
         jndx = jndx + nyvec(pj)
      END DO
      kndx = kndx + nzvec(pk)
   END DO
END IF

RETURN
END SUBROUTINE sigthr
