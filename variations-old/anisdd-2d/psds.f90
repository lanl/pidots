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
INTEGER :: bnc, i, n, pi, pj, tnx, tny, recrnk, sndrnk
INTEGER :: bcs, jndx, indx, indx1, ieq, ierr
INTEGER :: tag, tagi, tagj, tagm, tags, tagp
INTEGER, DIMENSION(2) :: ipak, coord
INTEGER, DIMENSION(3) :: istat
INTEGER, DIMENSION(:,:), ALLOCATABLE :: tempmat
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: temps
REAL*8, DIMENSION(:,:), ALLOCATABLE :: tempsi

INCLUDE 'mpif.h'

! Call for the processgrid to be created
CALL procgrid

! Now all the processes have their ranks and locations in the grid
! Can use info to distribute data: place topology on top of physical mesh

! Start with getting the appropriate number of cells for each process
IF (irank == root) THEN
   ! Allocate arrays to store
   ALLOCATE(nxvec(0:(npx-1)), nyvec(0:(npy-1)))
   DO i = 0, npx-1
      nxvec(i) = bnc(nxt,npx,i)
   END DO
   DO i = 0, npy-1
      nyvec(i) = bnc(nyt,npy,i)
   END DO
END IF

! Initialize ipak
ipak = 0

! Send to each process their size (nx,ny) so that they can allocate space
IF (irank == root) THEN
   ! Send the root's nrank to everyone so they know the message source in allcomm
   sndrnk = nrank
   CALL MPI_BCAST(sndrnk,1,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
   DO pj = 0, npy-1
      coord(2) = pj
      ipak(2) = nyvec(pj)
      tny = ipak(2)
      DO pi = 0, npx-1
         coord(1) = pi
         ipak(1) = nxvec(pi)
         tnx = ipak(1)
         ! Get the rank according to alcomm and set up tags
         CALL MPI_CART_RANK(allcomm,coord,recrnk,ierr)
         tag = 100 + recrnk
         IF (recrnk /= nrank) THEN
            ! Send the dimensions to the appropriate process in allcomm
            CALL MPI_SEND(ipak,2,MPI_INTEGER,recrnk,tag,allcomm,ierr)
         ELSE
            ! Set dimensions and allocate arays for later use and initialize
            nx = ipak(1)
            ny = ipak(2)
            bcs = apo*(nx+ny)
            ALLOCATE(dx(nx), dy(ny), mat(nx,ny), s(nx,ny,ng))
            ALLOCATE(psii(bcs,4))
            dx = 0.0
            dy = 0.0
            mat = 0
            s = 0.0
            psii = 0.0
         END IF
      END DO
   END DO
ELSE
   ! Get the root's sending rank for allcomm
   CALL MPI_BCAST(sndrnk,1,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
   ! Ranks that are not root need to receive data and allocate arrays
   tag = 100 + nrank
   CALL MPI_RECV(ipak,2,MPI_INTEGER,sndrnk,tag,allcomm,istat,ierr)
   nx = ipak(1)
   ny = ipak(2)
   bcs = apo*(nx+ny)
   ALLOCATE(dx(nx), dy(ny), mat(nx,ny), s(nx,ny,ng))
   ALLOCATE(psii(bcs,4))
   dx = 0.0
   dy = 0.0
   mat = 0
   s = 0.0
   psii = 0.0
END IF

! Root (from MPI_COMM_WORLD) has full data
! Loop over the grid and send data to the processes: dx, dy, dz, mat, s, psiin
IF (irank == root) THEN
   ! Organize the data to be sent, get the rank on the receiving end, send
   ! Use dummy index to know how much to send and where to start the buffer
   jndx = 0
   DO pj = 0, npy-1
      coord(2) = pj
      tny = nyvec(pj)
      indx = 0
      DO pi = 0, npx-1
         coord(1) = pi
         tnx = nxvec(pi)

         ! Allocate the matrices
         bcs = apo*(tnx+tny)
         ALLOCATE(tempmat(tnx,tny), temps(tnx,tny,ng))
         ALLOCATE(tempsi(bcs,4))
         tempmat(1:tnx,1:tny) = matt((indx+1):(indx+tnx),(jndx+1):(jndx+tny))
         temps(1:tnx,1:tny,:) = st((indx+1):(indx+tnx),(jndx+1):(jndx+tny),:)

         ! Initialize the angular flux BCs, then write over with actual data
         tempsi = 0.0
         IF (pj == 0 .AND. bc(1) == 2) THEN
            DO n = 1, apo
              indx1 = (n-1)*(tnx)
               DO i = 1, tnx
                  ieq = indx1 + i
                  tempsi(ieq,1) = posiny(indx+i,n,1)
                  tempsi(ieq,2) = posiny(indx+i,n,2)
               END DO
            END DO
         END IF
         IF (pj == npy-1 .AND. bc(2) == 2) THEN
            DO n = 1, apo
               indx1 = (n-1)*(tnx)
               DO i = 1, tnx
                 ieq = indx1 + i
                  tempsi(ieq,3) = neginy(indx+i,n,1)
                  tempsi(ieq,4) = neginy(indx+i,n,2)
               END DO
            END DO
         END IF
         IF (pi == 0 .AND. bc(3) == 2) THEN
            DO n = 1, apo
               indx1 = apo*(tnx) + (n-1)*(tny)
               DO i = 1, tny
                  ieq = indx1 + i
                  tempsi(ieq,1) = posinx(jndx+i,n,1)
                  tempsi(ieq,4) = posinx(jndx+i,n,2)
               END DO
            END DO
         END IF
         IF (pi == npx-1 .AND. bc(4) == 2) THEN
            DO n = 1, apo
               indx1 = apo*(tnx) + (n-1)*(tny)
               DO i = 1, tny
                  ieq = indx1 + i
                  tempsi(ieq,2) = neginx(jndx+i,n,1)
                  tempsi(ieq,3) = neginx(jndx+i,n,2)
               END DO
            END DO
         END IF

         ! Get the rank to send data to and set up the tags
         CALL MPI_CART_RANK(allcomm,coord,recrnk,ierr)
         tagi = 200 + recrnk
         tagj = 300 + recrnk
         tagm = 400 + recrnk
         tags = 500 + recrnk
         tagp = 600 + recrnk
         ! Send if not on the root's process, just copy if it is
         IF (recrnk /= sndrnk) THEN
            CALL MPI_SEND(dxt(indx+1),tnx,MPI_DOUBLE_PRECISION,recrnk,tagi,allcomm,ierr)
            CALL MPI_SEND(dyt(jndx+1),tny,MPI_DOUBLE_PRECISION,recrnk,tagj,allcomm,ierr)
            CALL MPI_SEND(tempmat,tnx*tny,MPI_INTEGER,recrnk,tagm,allcomm,ierr)
            CALL MPI_SEND(temps,tnx*tny*ng,MPI_DOUBLE_PRECISION,recrnk,tags,allcomm,ierr)
            CALL MPI_SEND(tempsi,bcs*4,MPI_DOUBLE_PRECISION,recrnk,tagp,allcomm,ierr)
         ELSE
            dx(1:nx) = dxt((indx+1):(indx+tnx))
            dy(1:ny) = dyt((jndx+1):(jndx+tny))
            mat(1:nx,1:ny) = tempmat(1:tnx,1:tny)
            s(1:nx,1:ny,:) = temps(1:tnx,1:tny,:)
            psii(1:(apo*(nx+ny)),1:4) = tempsi(1:bcs,1:4)
         END IF

         ! Update indx before cycling pi
         indx = indx + tnx

         ! Deallocate the matrices so they can change size afterwards
         DEALLOCATE(tempmat, temps, tempsi)
      END DO

      ! Update jndx before cycling pj
      jndx = jndx + tny
   END DO

ELSE
   ! Other processes receive their message and copy into their variables accordingly
   tagi = 200 + nrank
   tagj = 300 + nrank
   tagm = 400 + nrank
   tags = 500 + nrank
   tagp = 600 + nrank
   CALL MPI_RECV(dx,nx,MPI_DOUBLE_PRECISION,sndrnk,tagi,allcomm,istat,ierr)
   CALL MPI_RECV(dy,ny,MPI_DOUBLE_PRECISION,sndrnk,tagj,allcomm,istat,ierr)
   CALL MPI_RECV(mat,nx*ny,MPI_INTEGER,sndrnk,tagm,allcomm,istat,ierr)
   CALL MPI_RECV(s,nx*ny*ng,MPI_DOUBLE_PRECISION,sndrnk,tags,allcomm,istat,ierr)
   CALL MPI_RECV(psii,bcs*4,MPI_DOUBLE_PRECISION,sndrnk,tagp,allcomm,istat,ierr)
END IF

RETURN
END SUBROUTINE psds 
