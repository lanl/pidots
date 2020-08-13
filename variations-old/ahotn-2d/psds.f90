SUBROUTINE psds(root)

!-------------------------------------------------------------
!
! Sets up the parallel subdomains
! Calls for the 2D topology creation, then distributes problem
!  data according to their coordinates in topology
!
!-------------------------------------------------------------

USE invar
USE totvar
IMPLICIT NONE
INTEGER, INTENT(IN) :: root
INTEGER :: bnc, i, ierr, pi, pj, sndrnk, recrnk, tnx, tny, n, indx, jndx, indx1, indx2
INTEGER :: tag, tagi, tagj, tagm, tags, tagp1, tagp2, tagp3, tagp4
INTEGER, DIMENSION(2) :: ipak, coord
INTEGER, DIMENSION(3) :: istat
INTEGER, DIMENSION(:,:), ALLOCATABLE :: tempmat
REAL*8, DIMENSION(:,:,:,:,:), ALLOCATABLE :: temps
REAL*8, DIMENSION(:), ALLOCATABLE :: tmpsi1, tmpsi2, tmpsi3, tmpsi4

INCLUDE 'mpif.h'

! Call for the process grid to be created
CALL procgrid

! Now all processes have their ranks and locations in the grid
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
   ! Send the root's nrank to everyone so they know the message source
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
         ! Get the rank according to allcomm and set up tags
         CALL MPI_CART_RANK(allcomm,coord,recrnk,ierr)
         tag = 100 + recrnk
         IF (recrnk /= nrank) THEN
            ! Send the dimensions to the appropriate process in allcomm
            CALL MPI_SEND(ipak,2,MPI_INTEGER,recrnk,tag,allcomm,ierr)
         ELSE
            ! Set dimensions and allocate arrays for later use and initialize
            nx = ipak(1)
            ny = ipak(2)
            ALLOCATE(dx(nx), dy(ny), mat(nx,ny), s(nx,ny,0:lambda,0:lambda,ng))
            ALLOCATE(psi1(apo*(nx+ny)*order), psi2(apo*(nx+ny)*order))
            ALLOCATE(psi3(apo*(nx+ny)*order), psi4(apo*(nx+ny)*order))
            dx = 0.0
            dy = 0.0
            mat = 0
            s = 0.0
            psi1 = 0.0
            psi2 = 0.0
            psi3 = 0.0
            psi4 = 0.0
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
   ALLOCATE(dx(nx), dy(ny), mat(nx,ny), s(nx,ny,0:lambda,0:lambda,ng))
   ALLOCATE(psi1(apo*(nx+ny)*order), psi2(apo*(nx+ny)*order))
   ALLOCATE(psi3(apo*(nx+ny)*order), psi4(apo*(nx+ny)*order))
   dx = 0.0
   dy = 0.0
   mat = 0
   s = 0.0
   psi1 = 0.0
   psi2 = 0.0
   psi3 = 0.0
   psi4 = 0.0
END IF

! Root (from MPI_COMM_WORLD) has full data
! Loop over the grid and send data to the processes: dx, dy, mat, s, psi#
IF (irank == root) THEN
   sndrnk = nrank
   jndx = 0
   ! Organize the data to be sent, get the rank on the receiving end, send
   DO pj = 0, npy-1
      coord(2) = pj
      tny = nyvec(pj)
      indx = 0
      DO pi = 0, npx-1
         coord(1) = pi
         tnx = nxvec(pi)

         ! Allocate the matrices
         i = apo*(tnx+tny)*order
         ALLOCATE(tempmat(tnx,tny), temps(tnx,tny,0:lambda,0:lambda,ng))
         ALLOCATE(tmpsi1(i), tmpsi2(i), tmpsi3(i), tmpsi4(i))

         tempmat(1:tnx,1:tny) = matt((indx+1):(indx+tnx),(jndx+1):(jndx+tny))
         temps(1:tnx,1:tny,:,:,:) = st((indx+1):(indx+tnx),(jndx+1):(jndx+tny),:,:,:)

         ! Initialize the angular flux BCs, then write over with actual data
         tmpsi1 = 0.0
         tmpsi2 = 0.0
         tmpsi3 = 0.0
         tmpsi4 = 0.0
         ! Prepare the ybc/xbc array
         IF (pj == 0) THEN
            DO n = 1, apo
               indx1 = (n-1)*nxt*order + indx*order
               indx2 = (n-1)*tnx*order
               tmpsi1((indx2+1):(indx2+tnx*order)) = psi1t((indx1+1):(indx1+tnx*order))
               tmpsi2((indx2+1):(indx2+tnx*order)) = psi2t((indx1+1):(indx1+tnx*order))
            END DO
         END IF
         IF (pj == npy-1) THEN
            DO n = 1, apo
               indx1 = (n-1)*nxt*order + indx*order
               indx2 = (n-1)*tnx*order
               tmpsi3((indx2+1):(indx2+tnx*order)) = psi3t((indx1+1):(indx1+tnx*order))
               tmpsi4((indx2+1):(indx2+tnx*order)) = psi4t((indx1+1):(indx1+tnx*order))
            END DO
         END IF
         IF (pi == 0) THEN
            DO n = 1, apo
               indx1 = apo*nxt*order + (n-1)*nyt*order + jndx*order
               indx2 = apo*tnx*order + (n-1)*tny*order
               tmpsi1((indx2+1):(indx2+tny*order)) = psi1t((indx1+1):(indx1+tny*order))
               tmpsi4((indx2+1):(indx2+tny*order)) = psi4t((indx1+1):(indx1+tny*order))
            END DO
         END IF
         IF (pi == npx-1) THEN
            DO n = 1, apo
               indx1 = apo*nxt*order + (n-1)*nyt*order + jndx*order
               indx2 = apo*tnx*order + (n-1)*tny*order
               tmpsi2((indx2+1):(indx2+tny*order)) = psi2t((indx1+1):(indx1+tny*order))
               tmpsi3((indx2+1):(indx2+tny*order)) = psi3t((indx1+1):(indx1+tny*order))
            END DO
         END IF

         ! Get the rank to send data to and set up the tags
         CALL MPI_CART_RANK(allcomm,coord,recrnk,ierr)
         tagi = 200 + recrnk
         tagj = 300 + recrnk
         tagm = 400 + recrnk
         tags = 500 + recrnk
         tagp1 = 600 + recrnk
         tagp2 = 700 + recrnk
         tagp3 = 800 + recrnk
         tagp4 = 900 + recrnk
         ! Send if not on the root's process, just copy if it is
         IF (recrnk /= sndrnk) THEN
            CALL MPI_SEND(dxt(indx+1),tnx,MPI_DOUBLE_PRECISION,recrnk,tagi,allcomm,ierr)
            CALL MPI_SEND(dyt(jndx+1),tny,MPI_DOUBLE_PRECISION,recrnk,tagj,allcomm,ierr)
            CALL MPI_SEND(tempmat,tnx*tny,MPI_INTEGER,recrnk,tagm,allcomm,ierr)
            CALL MPI_SEND(temps,tnx*tny*ordsq*ng,MPI_DOUBLE_PRECISION,recrnk,tags,allcomm,ierr)
            CALL MPI_SEND(tmpsi1,apo*(tnx+tny)*order,MPI_DOUBLE_PRECISION,recrnk,tagp1,allcomm,ierr)
            CALL MPI_SEND(tmpsi2,apo*(tnx+tny)*order,MPI_DOUBLE_PRECISION,recrnk,tagp2,allcomm,ierr)
            CALL MPI_SEND(tmpsi3,apo*(tnx+tny)*order,MPI_DOUBLE_PRECISION,recrnk,tagp3,allcomm,ierr)
            CALL MPI_SEND(tmpsi4,apo*(tnx+tny)*order,MPI_DOUBLE_PRECISION,recrnk,tagp4,allcomm,ierr)
         ELSE
            ! For the process with all the info, cut down to what is needed, ignore the rest
            dx(1:nx) = dxt((indx+1):(indx+tnx))
            dy(1:ny) = dyt((jndx+1):(jndx+tny))
            mat(1:nx,1:ny) = tempmat(1:tnx,1:tny)
            s(1:nx,1:ny,:,:,:) = temps(1:tnx,1:tny,:,:,:)
            psi1(1:(apo*(nx+ny)*order)) = tmpsi1(1:(apo*(tnx+tny)*order))
            psi2(1:(apo*(nx+ny)*order)) = tmpsi2(1:(apo*(tnx+tny)*order))
            psi3(1:(apo*(nx+ny)*order)) = tmpsi3(1:(apo*(tnx+tny)*order))
            psi4(1:(apo*(nx+ny)*order)) = tmpsi4(1:(apo*(tnx+tny)*order))
         END IF
         
         ! Update indx before cycling pi
         indx = indx + tnx

         ! Deallocate the matrices so they can change size afterwards
         DEALLOCATE(tempmat, temps, tmpsi1, tmpsi2, tmpsi3, tmpsi4)
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
   tagp1 = 600 + nrank
   tagp2 = 700 + nrank
   tagp3 = 800 + nrank
   tagp4 = 900 + nrank
   CALL MPI_RECV(dx,nx,MPI_DOUBLE_PRECISION,sndrnk,tagi,allcomm,istat,ierr)
   CALL MPI_RECV(dy,ny,MPI_DOUBLE_PRECISION,sndrnk,tagj,allcomm,istat,ierr)
   CALL MPI_RECV(mat,nx*ny,MPI_INTEGER,sndrnk,tagm,allcomm,istat,ierr)
   CALL MPI_RECV(s,nx*ny*ordsq*ng,MPI_DOUBLE_PRECISION,sndrnk,tags,allcomm,istat,ierr)
   CALL MPI_RECV(psi1,apo*(nx+ny)*order,MPI_DOUBLE_PRECISION,sndrnk,tagp1,allcomm,istat,ierr)
   CALL MPI_RECV(psi2,apo*(nx+ny)*order,MPI_DOUBLE_PRECISION,sndrnk,tagp2,allcomm,istat,ierr)
   CALL MPI_RECV(psi3,apo*(nx+ny)*order,MPI_DOUBLE_PRECISION,sndrnk,tagp3,allcomm,istat,ierr)
   CALL MPI_RECV(psi4,apo*(nx+ny)*order,MPI_DOUBLE_PRECISION,sndrnk,tagp4,allcomm,istat,ierr)
END IF

RETURN
END SUBROUTINE psds 
