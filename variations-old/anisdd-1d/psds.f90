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
INTEGER :: bnc, i, pi, tnx, indx, ierr, tag, tagi, tagm, tags, tagp
INTEGER, DIMENSION(3) :: istat
INTEGER, DIMENSION(:), ALLOCATABLE :: tempmat
REAL*8, DIMENSION(:,:), ALLOCATABLE :: temps
REAL*8, DIMENSION(apo,2) :: tempsi

INCLUDE 'mpif.h'

! Call for the processgrid to be created

! Now all the processes have their ranks and locations in the grid
! Can use info to distribute data: place topology on top of physical mesh

! Start with getting the appropriate number of cells for each process
IF (irank == root) THEN
   ! Allocate arrays to store
   ALLOCATE(nxvec(0:(npx-1)))
   DO i = 0, npx-1
      nxvec(i) = bnc(nxt,npx,i)
   END DO
END IF

! Send to each process their size (nx) so that they can allocate space
IF (irank == root) THEN
   DO pi = 0, npx-1
      tnx = nxvec(pi)
      tag = 100 + pi
      IF (pi /= irank) THEN
         ! Send the dimensions to the appropriate process in allcomm
         CALL MPI_SEND(tnx,1,MPI_INTEGER,pi,tag,MPI_COMM_WORLD,ierr)
      ELSE
         ! Set dimensions and allocate arays for later use and initialize
         nx = tnx
         ALLOCATE(dx(nx), mat(nx), s(nx,ng))
         ALLOCATE(psii(apo,2))
         dx = 0.0
         mat = 0
         s = 0.0
         psii = 0.0
      END IF
   END DO
ELSE
   tag = 100 + irank
   CALL MPI_RECV(tnx,1,MPI_INTEGER,root,tag,MPI_COMM_WORLD,istat,ierr)
   nx = tnx
   ALLOCATE(dx(nx), mat(nx), s(nx,ng))
   ALLOCATE(psii(apo,2))
   dx = 0.0
   mat = 0
   s = 0.0
   psii = 0.0
END IF

! Root (from MPI_COMM_WORLD) has full data
! Loop over the grid and send data to the processes: dx, mat, s, psii
IF (irank == root) THEN
   ! Use dummy index to know how much to send and where to start the buffer
   indx = 0
   DO pi = 0, npx-1
      tnx = nxvec(pi)

      ! Allocate the matrices
      ALLOCATE(tempmat(tnx), temps(tnx,ng))
      tempmat(1:tnx) = matt((indx+1):(indx+tnx))
      temps(1:tnx,:) = st((indx+1):(indx+tnx),:)

      ! Initialize the angular flux BCs, then write over with actual data
      tempsi = 0.0
      IF (pi == 0 .AND. bc(1) == 2) THEN
         tempsi(:,1) = posinx
      END IF
      IF (pi == npx-1 .AND. bc(2) == 2) THEN
         tempsi(:,2) = neginx
      END IF

      tagi = 200 + pi
      tagm = 300 + pi
      tags = 400 + pi
      tagp = 500 + pi
      ! Send if not on the root's process, just copy if it is
      IF (pi /= irank) THEN
         CALL MPI_SEND(dxt(indx+1),tnx,MPI_DOUBLE_PRECISION,pi,tagi,MPI_COMM_WORLD,ierr)
         CALL MPI_SEND(tempmat,tnx,MPI_INTEGER,pi,tagm,MPI_COMM_WORLD,ierr)
         CALL MPI_SEND(temps,tnx*ng,MPI_DOUBLE_PRECISION,pi,tags,MPI_COMM_WORLD,ierr)
         CALL MPI_SEND(tempsi,apo*2,MPI_DOUBLE_PRECISION,pi,tagp,MPI_COMM_WORLD,ierr)
      ELSE
         dx(1:nx) = dxt(1:tnx)
         mat(1:nx) = tempmat(1:tnx)
         s(1:nx,:) = temps(1:tnx,:)
         psii = tempsi
      END IF

      ! Update indx before cycling pi
      indx = indx + tnx

      ! Deallocate the matrices so they can change size afterwards
      DEALLOCATE(tempmat, temps)
   END DO

ELSE
   ! Other processes receive their message and copy into their variables accordingly
   tagi = 200 + irank
   tagm = 300 + irank
   tags = 400 + irank
   tagp = 500 + irank
   CALL MPI_RECV(dx,nx,MPI_DOUBLE_PRECISION,root,tagi,MPI_COMM_WORLD,istat,ierr)
   CALL MPI_RECV(mat,nx,MPI_INTEGER,root,tagm,MPI_COMM_WORLD,istat,ierr)
   CALL MPI_RECV(s,nx*ng,MPI_DOUBLE_PRECISION,root,tags,MPI_COMM_WORLD,istat,ierr)
   CALL MPI_RECV(psii,apo*2,MPI_DOUBLE_PRECISION,root,tagp,MPI_COMM_WORLD,istat,ierr)
END IF

RETURN
END SUBROUTINE psds 
