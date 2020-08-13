SUBROUTINE psigthr

!-------------------------------------------------------------
!
! Gathers all the rpave vectors and puts them into rpsi 
! 
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER :: root2, tag, pk, pj, pi, trank, ieq, n, ierr, i, j, k
INTEGER, DIMENSION(3) :: coord, istat
REAL*8, DIMENSION(bcs,8) :: tmpsi
REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE :: afz, afy, afx

INCLUDE 'mpif.h'

IF (irank == root) THEN
   ALLOCATE(afz(nxt,nyt,apo,8),afy(nxt,nzt,apo,8),afx(nyt,nzt,apo,8))
END IF

! Setup for parallel
IF (irank == root) root2 = nrank
CALL MPI_BCAST(root2,1,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)

! Bring all the blocks' solutions for all the groups back to root for output
IF (nrank /= root2) THEN
   tag = 100 + nrank
   CALL MPI_SEND(psio(:,:),bcs*8,MPI_DOUBLE_PRECISION,root2,tag,allcomm,ierr)
ELSE
   DO pk = 0, npz-1
      coord(1) = pk
      DO pj = 0, npy-1
         coord(2) = pj
         DO pi = 0, npx-1
            coord(3) = pi
            CALL MPI_CART_RANK(allcomm,coord,trank,ierr)
            IF (trank /= nrank) THEN
               tmpsi = 0.0
               tag = 100 + trank
               CALL MPI_RECV(tmpsi,bcs*8,MPI_DOUBLE_PRECISION,trank,tag,allcomm,istat,ierr)
               IF (pk == 0 .AND. bc(1) == 1) THEN
                  DO n = 1, apo
                     DO j = 1, ny
                        DO i = 1, nx
                           ieq = (n-1)*xys + (j-1)*nx + i
                           afz(pi*nx+i,pj*ny+j,n,1) = tmpsi(ieq,5)
                           afz(pi*nx+i,pj*ny+j,n,2) = tmpsi(ieq,6)
                           afz(pi*nx+i,pj*ny+j,n,3) = tmpsi(ieq,7)
                           afz(pi*nx+i,pj*ny+j,n,4) = tmpsi(ieq,8)
                        END DO
                     END DO
                  END DO
               END IF
               IF (pk == npz-1 .AND. bc(2) == 1) THEN
                  DO n = 1, apo
                     DO j = 1, ny
                        DO i = 1, nx
                           ieq = (n-1)*xys + (j-1)*nx + i
                           afz(pi*nx+i,pj*ny+j,n,5) = tmpsi(ieq,1)
                           afz(pi*nx+i,pj*ny+j,n,6) = tmpsi(ieq,2)
                           afz(pi*nx+i,pj*ny+j,n,7) = tmpsi(ieq,3)
                           afz(pi*nx+i,pj*ny+j,n,8) = tmpsi(ieq,4)
                        END DO
                     END DO
                  END DO
               END IF

               IF (pj == 0 .AND. bc(3) == 1) THEN
                  DO n = 1, apo
                     DO k = 1, nz
                        DO i = 1, nx
                           ieq = apo*xys + (n-1)*xzs + (k-1)*nx + i
                           afy(pi*nx+i,pk*nz+k,n,1) = tmpsi(ieq,4)
                           afy(pi*nx+i,pk*nz+k,n,2) = tmpsi(ieq,3)
                           afy(pi*nx+i,pk*nz+k,n,5) = tmpsi(ieq,8)
                           afy(pi*nx+i,pk*nz+k,n,6) = tmpsi(ieq,7)
                        END DO
                     END DO
                  END DO
               END IF
               IF (pj == npy-1 .AND. bc(4) == 1) THEN
                  DO n = 1, apo
                     DO k = 1, nz
                        DO i = 1, nx
                           ieq = apo*xys + (n-1)*xzs + (k-1)*nx + i
                           afy(pi*nx+i,pk*nz+k,n,4) = tmpsi(ieq,1)
                           afy(pi*nx+i,pk*nz+k,n,3) = tmpsi(ieq,2)
                           afy(pi*nx+i,pk*nz+k,n,8) = tmpsi(ieq,5)
                           afy(pi*nx+i,pk*nz+k,n,7) = tmpsi(ieq,6)
                        END DO
                     END DO
                  END DO
               END IF

               IF (pi == 0 .AND. bc(5) == 1) THEN
                  DO n = 1, apo
                     DO k = 1, nz
                        DO j = 1, ny
                           ieq = apo*(xys+xzs) + (n-1)*yzs + (k-1)*ny + i
                           afx(pj*ny+j,pk*nz+k,n,1) = tmpsi(ieq,2)
                           afx(pj*ny+j,pk*nz+k,n,4) = tmpsi(ieq,3)
                           afx(pj*ny+j,pk*nz+k,n,5) = tmpsi(ieq,6)
                           afx(pj*ny+j,pk*nz+k,n,8) = tmpsi(ieq,7)
                        END DO
                     END DO
                  END DO
               END IF
               IF (pi == npx-1 .AND. bc(6) == 1) THEN
                  DO n = 1, apo
                     DO k = 1, nz
                        DO j = 1, ny
                           ieq = apo*(xys+xzs) + (n-1)*yzs + (k-1)*ny + i
                           afx(pj*ny+j,pk*nz+k,n,2) = tmpsi(ieq,1)
                           afx(pj*ny+j,pk*nz+k,n,3) = tmpsi(ieq,4)
                           afx(pj*ny+j,pk*nz+k,n,6) = tmpsi(ieq,5)
                           afx(pj*ny+j,pk*nz+k,n,7) = tmpsi(ieq,8)
                        END DO
                     END DO
                  END DO
               END IF

            ELSE IF (trank == nrank) THEN
               ! Know that root is 0,0,0
               IF (bc(1) == 1) THEN
                  DO n = 1, apo
                     DO j = 1, ny
                        DO i = 1, nx
                           ieq = (n-1)*xys + (j-1)*nx + i
                           afz(i,j,n,1) = psio(ieq,5)
                           afz(i,j,n,2) = psio(ieq,6)
                           afz(i,j,n,3) = psio(ieq,7)
                           afz(i,j,n,4) = psio(ieq,8)
                        END DO
                     END DO
                  END DO
               END IF
               IF (bc(3) == 1) THEN
                  DO n = 1, apo
                     DO k = 1, nz
                        DO i = 1, nx
                           ieq = apo*xys + (n-1)*xzs + (k-1)*nx + i
                           afy(i,k,n,1) = psio(ieq,4)
                           afy(i,k,n,2) = psio(ieq,3)
                           afy(i,k,n,5) = psio(ieq,8)
                           afy(i,k,n,6) = psio(ieq,7)
                        END DO
                     END DO
                  END DO
               END IF
               IF (bc(5) == 1) THEN
                  DO n = 1, apo
                     DO k = 1, nz
                        DO j = 1, ny
                           ieq = apo*(xys+xzs) + (n-1)*yzs + (k-1)*ny + i
                           afx(j,k,n,1) = psio(ieq,2)
                           afx(j,k,n,4) = psio(ieq,3)
                           afx(j,k,n,5) = psio(ieq,6)
                           afx(j,k,n,8) = psio(ieq,7)
                        END DO
                     END DO
                  END DO
               END IF
            END IF
         END DO
      END DO
   END DO
END IF

RETURN
END SUBROUTINE psigthr
