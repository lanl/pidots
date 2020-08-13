SUBROUTINE psds

!-------------------------------------------------------------
!
! Sets up the parallel subdomains
! Calls for the 3D topology creation, then distributes problem
!  data according to their coordinates in topology
!
!-------------------------------------------------------------

USE invar
USE totvar
IMPLICIT NONE
INTEGER :: i, j, n, pi, pj, pk, recrnk, sndrnk
INTEGER :: xys, xzs, yzs, bcs, kndx, jndx, indx, indx1, ieq, ierr
INTEGER :: tag, tagi, tagj, tagk, tagm, tags, tagp
INTEGER, DIMENSION(3) :: coord, istat
INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: tempmat
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: temps
REAL*8, DIMENSION(:,:), ALLOCATABLE :: tempsi

INCLUDE 'mpif.h'

! Call for the processgrid to be created
CALL procgrid

! Now all the processes have their ranks and locations in the grid
! Can use info to distribute data: place topology on top of physical mesh
xys = nx*ny
xzs = nx*nz
yzs = ny*nz
bcs = apo*(xys+xzs+yzs)
ALLOCATE(dx(nx), dy(ny), dz(nz), mat(nx,ny,nz), s(nx,ny,nz))
ALLOCATE(psii(bcs,8,ns))
dx = 0.0
dy = 0.0
dz = 0.0
mat = 0
s = 0.0
psii = 0.0

! Root (from MPI_COMM_WORLD) has full data
IF (irank == root) sndrnk = nrank
CALL MPI_BCAST(sndrnk,1,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)

! Loop over the grid and send data to the processes: dx, dy, dz, mat, s, psiin
IF (irank == root) THEN
   ! Allocate temp matrices
   ALLOCATE(tempmat(nx,ny,nz), temps(nx,ny,nz))
   ALLOCATE(tempsi(bcs,8))

   ! Organize the data to be sent, get the rank on the receiving end, send
   ! Use dummy index to know how much to send and where to start the buffer
   DO pk = 0, npz-1
      coord(1) = pk
      kndx = pk*nz
      DO pj = 0, npy-1
         coord(2) = pj
         jndx = pj*ny
         DO pi = 0, npx-1
            coord(3) = pi
            indx = pi*nx

            ! Set up temporary matrices
            tempmat = matt((indx+1):(indx+nx),(jndx+1):(jndx+ny),(kndx+1):(kndx+nz))
            temps = st((indx+1):(indx+nx),(jndx+1):(jndx+ny),(kndx+1):(kndx+nz))

            ! Initialize the angular flux BCs, then write over with actual data
            tempsi = 0.0
            IF (pk == 0 .AND. bc(1) == 2) THEN
               DO n = 1, apo
                  indx1 = (n-1)*xys
                  DO j = 1, ny
                     DO i = 1, nx
                        ieq = indx1 + (j-1)*nx + i
                        tempsi(ieq,1) = posinz(indx+i,jndx+j,n,1)
                        tempsi(ieq,2) = posinz(indx+i,jndx+j,n,2)
                        tempsi(ieq,3) = posinz(indx+i,jndx+j,n,3)
                        tempsi(ieq,4) = posinz(indx+i,jndx+j,n,4)
                     END DO
                  END DO
               END DO
            END IF
            IF (pk == npz-1 .AND. bc(2) == 2) THEN
               DO n = 1, apo
                  indx1 = (n-1)*xys
                  DO j = 1, ny
                     DO i = 1, nx
                        ieq = indx1 + (j-1)*nx + i
                        tempsi(ieq,5) = neginz(indx+i,jndx+j,n,1)
                        tempsi(ieq,6) = neginz(indx+i,jndx+j,n,2)
                        tempsi(ieq,7) = neginz(indx+i,jndx+j,n,3)
                        tempsi(ieq,8) = neginz(indx+i,jndx+j,n,4)
                     END DO
                  END DO
               END DO
            END IF
            IF (pj == 0 .AND. bc(3) == 2) THEN
               DO n = 1, apo
                  indx1 = apo*xys + (n-1)*xzs
                  DO j = 1, nz
                     DO i = 1, nx
                        ieq = indx1 + (j-1)*nx + i
                        tempsi(ieq,1) = posiny(indx+i,kndx+j,n,1)
                        tempsi(ieq,2) = posiny(indx+i,kndx+j,n,2)
                        tempsi(ieq,5) = posiny(indx+i,kndx+j,n,3)
                        tempsi(ieq,6) = posiny(indx+i,kndx+j,n,4)
                     END DO
                  END DO
               END DO
            END IF
            IF (pj == npy-1 .AND. bc(4) == 2) THEN
               DO n = 1, apo
                  indx1 = apo*xys + (n-1)*xzs
                  DO j = 1, nz
                     DO i = 1, nx
                        ieq = indx1 + (j-1)*nx + i
                        tempsi(ieq,3) = neginy(indx+i,kndx+j,n,1)
                        tempsi(ieq,4) = neginy(indx+i,kndx+j,n,2)
                        tempsi(ieq,7) = neginy(indx+i,kndx+j,n,3)
                        tempsi(ieq,8) = neginy(indx+i,kndx+j,n,4)
                     END DO
                  END DO
               END DO
            END IF
            IF (pi == 0 .AND. bc(5) == 2) THEN
               DO n = 1, apo
                  indx1 = apo*(xys+xzs) + (n-1)*yzs
                  DO j = 1, nz
                     DO i = 1, ny
                        ieq = indx1 + (j-1)*ny + i
                        tempsi(ieq,1) = posinx(jndx+i,kndx+j,n,1)
                        tempsi(ieq,4) = posinx(jndx+i,kndx+j,n,2)
                        tempsi(ieq,5) = posinx(jndx+i,kndx+j,n,3)
                        tempsi(ieq,8) = posinx(jndx+i,kndx+j,n,4)
                     END DO
                  END DO
               END DO
            END IF
            IF (pi == npx-1 .AND. bc(6) == 2) THEN
               DO n = 1, apo
                  indx1 = apo*(xys+xzs) + (n-1)*yzs
                  DO j = 1, nz
                     DO i = 1, ny
                        ieq = indx1 + (j-1)*ny + i
                        tempsi(ieq,2) = neginx(jndx+i,kndx+j,n,1)
                        tempsi(ieq,3) = neginx(jndx+i,kndx+j,n,2)
                        tempsi(ieq,6) = neginx(jndx+i,kndx+j,n,3)
                        tempsi(ieq,7) = neginx(jndx+i,kndx+j,n,4)
                     END DO
                  END DO
               END DO
            END IF

            ! Get the rank to send data to and set up the tags
            CALL MPI_CART_RANK(allcomm,coord,recrnk,ierr)
            tagi = 200 + recrnk
            tagj = 300 + recrnk
            tagk = 400 + recrnk
            tagm = 500 + recrnk
            tags = 600 + recrnk
            tagp = 700 + recrnk
            ! Send if not on the root's process, just copy if it is
            IF (recrnk /= sndrnk) THEN
               CALL MPI_SEND(dxt(indx+1),nx,MPI_DOUBLE_PRECISION,recrnk,tagi,allcomm,ierr)
               CALL MPI_SEND(dyt(jndx+1),ny,MPI_DOUBLE_PRECISION,recrnk,tagj,allcomm,ierr)
               CALL MPI_SEND(dzt(kndx+1),nz,MPI_DOUBLE_PRECISION,recrnk,tagk,allcomm,ierr)
               CALL MPI_SEND(tempmat,nx*ny*nz,MPI_INTEGER,recrnk,tagm,allcomm,ierr)
               CALL MPI_SEND(temps,nx*ny*nz,MPI_DOUBLE_PRECISION,recrnk,tags,allcomm,ierr)
               CALL MPI_SEND(tempsi,bcs*8,MPI_DOUBLE_PRECISION,recrnk,tagp,allcomm,ierr)
            ELSE
               dx = dxt((indx+1):(indx+nx))
               dy = dyt((jndx+1):(jndx+ny))
               dz = dzt((kndx+1):(kndx+nz))
               mat = tempmat
               s = temps
               psii(:,:,1) = tempsi
            END IF
         END DO
      END DO
   END DO
   ! Deallocate
   DEALLOCATE(tempmat,temps,tempsi)
ELSE
   ! Other processes receive their message and copy into their variables accordingly
   tagi = 200 + nrank
   tagj = 300 + nrank
   tagk = 400 + nrank
   tagm = 500 + nrank
   tags = 600 + nrank
   tagp = 700 + nrank
   CALL MPI_RECV(dx,nx,MPI_DOUBLE_PRECISION,sndrnk,tagi,allcomm,istat,ierr)
   CALL MPI_RECV(dy,ny,MPI_DOUBLE_PRECISION,sndrnk,tagj,allcomm,istat,ierr)
   CALL MPI_RECV(dz,nz,MPI_DOUBLE_PRECISION,sndrnk,tagk,allcomm,istat,ierr)
   CALL MPI_RECV(mat,nx*ny*nz,MPI_INTEGER,sndrnk,tagm,allcomm,istat,ierr)
   CALL MPI_RECV(s,nx*ny*nz,MPI_DOUBLE_PRECISION,sndrnk,tags,allcomm,istat,ierr)
   CALL MPI_RECV(psii(:,:,1),bcs*8,MPI_DOUBLE_PRECISION,sndrnk,tagp,allcomm,istat,ierr)
END IF

RETURN
END SUBROUTINE psds 
