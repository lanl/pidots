SUBROUTINE rpsigthr(rpsi,rpave)

!-------------------------------------------------------------
!
! Gathers all the rpave vectors and puts them into rpsi 
! 
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER :: root2, tag, pk, pj, pi, trank, ieq, n, ierr
INTEGER, DIMENSION(3) :: coord, istat
REAL*8, DIMENSION(3*apo,8) :: tmpsi
REAL*8, DIMENSION(bcs,8), INTENT(INOUT) :: rpsi
REAL*8, DIMENSION(3*apo,8), INTENT(IN) :: rpave

INCLUDE 'mpif.h'

! Setup for parallel
IF (irank == root) root2 = nrank
CALL MPI_BCAST(root2,1,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)

! Bring all the blocks' solutions for all the groups back to root for output
IF (nrank /= root2) THEN
   tag = 100 + nrank
   CALL MPI_SEND(rpave,3*apo*8,MPI_DOUBLE_PRECISION,root2,tag,allcomm,ierr)
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
               CALL MPI_RECV(tmpsi,3*apo*8,MPI_DOUBLE_PRECISION,trank,tag,allcomm,istat,ierr)
               IF (pk == 0) THEN
                  DO n = 1, apo
                     ieq = (n-1)*xys + pj*nx + pi + 1
                     rpsi(ieq,5) = tmpsi(n,5)
                     rpsi(ieq,6) = tmpsi(n,6)
                     rpsi(ieq,7) = tmpsi(n,7)
                     rpsi(ieq,8) = tmpsi(n,8)
                  END DO
               ELSE
                  DO n = 1, apo
                     ieq = (n-1)*xys + pj*nx + pi + 1
                     rpsi(ieq,1) = tmpsi(n,1)
                     rpsi(ieq,2) = tmpsi(n,2)
                     rpsi(ieq,3) = tmpsi(n,3)
                     rpsi(ieq,4) = tmpsi(n,4)
                  END DO
               END IF
               IF (pj == 0) THEN
                  DO n = 1, apo
                     ieq = apo*xys + (n-1)*xzs + pk*nx + pi + 1
                     rpsi(ieq,3) = tmpsi(apo+n,3)
                     rpsi(ieq,4) = tmpsi(apo+n,4)
                     rpsi(ieq,7) = tmpsi(apo+n,7)
                     rpsi(ieq,8) = tmpsi(apo+n,8)
                  END DO
               ELSE
                  DO n = 1, apo
                     ieq = apo*xys + (n-1)*xzs + pk*nx + pi + 1
                     rpsi(ieq,1) = tmpsi(apo+n,1)
                     rpsi(ieq,2) = tmpsi(apo+n,2)
                     rpsi(ieq,5) = tmpsi(apo+n,5)
                     rpsi(ieq,6) = tmpsi(apo+n,6)
                  END DO
               END IF
               IF (pi == 0) THEN
                  DO n = 1, apo
                     ieq = apo*(xys+xzs) + (n-1)*yzs + pk*ny + pj + 1
                     rpsi(ieq,2) = tmpsi(2*apo+n,2)
                     rpsi(ieq,3) = tmpsi(2*apo+n,3)
                     rpsi(ieq,6) = tmpsi(2*apo+n,6)
                     rpsi(ieq,7) = tmpsi(2*apo+n,7)
                  END DO
               ELSE
                  DO n = 1, apo
                     ieq = apo*(xys+xzs) + (n-1)*yzs + pk*ny + pj + 1
                     rpsi(ieq,1) = tmpsi(2*apo+n,1)
                     rpsi(ieq,4) = tmpsi(2*apo+n,4)
                     rpsi(ieq,5) = tmpsi(2*apo+n,5)
                     rpsi(ieq,8) = tmpsi(2*apo+n,8)
                  END DO
               END IF
            ELSE IF (trank == nrank) THEN
               ! Know that root is 0,0,0
               DO n = 1, apo
                  ieq = (n-1)*xys + 1
                  rpsi(ieq,5) = rpave(n,5)
                  rpsi(ieq,6) = rpave(n,6)
                  rpsi(ieq,7) = rpave(n,7)
                  rpsi(ieq,8) = rpave(n,8)

                  ieq = apo*xys + (n-1)*xzs + 1
                  rpsi(ieq,3) = rpave(apo+n,3)
                  rpsi(ieq,4) = rpave(apo+n,4)
                  rpsi(ieq,7) = rpave(apo+n,7)
                  rpsi(ieq,8) = rpave(apo+n,8)

                  ieq = apo*(xys+xzs) + (n-1)*yzs + 1
                  rpsi(ieq,2) = rpave(2*apo+n,2)
                  rpsi(ieq,3) = rpave(2*apo+n,3)
                  rpsi(ieq,6) = rpave(2*apo+n,6)
                  rpsi(ieq,7) = rpave(2*apo+n,7)
               END DO
            END IF
         END DO
      END DO
   END DO
END IF

RETURN
END SUBROUTINE rpsigthr
