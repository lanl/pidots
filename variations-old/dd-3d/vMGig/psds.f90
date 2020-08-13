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
INTEGER :: bnc, i, j, n, pi, pj, pk, tnx, tny, tnz, recrnk, sndrnk
INTEGER :: bcs, kndx, jndx, indx, indx1, ieq, ierr
INTEGER :: tag, tagi, tagj, tagk, tagm, tags, tagp, tagf
INTEGER, DIMENSION(3) :: ipak, coord, istat
INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: tempmat
REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE :: temps
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: tempf
REAL*8, DIMENSION(:,:), ALLOCATABLE :: tempsi

INCLUDE 'mpif.h'

! Call for the processgrid to be created
CALL procgrid

! Now all the processes have their ranks and locations in the grid
! Can use info to distribute data: place topology on top of physical mesh

! Start with getting the appropriate number of cells for each process
IF (irank == root) THEN
   ! Allocate arrays to store
   ALLOCATE(nxvec(0:(npx-1)), nyvec(0:(npy-1)), nzvec(0:(npz-1)))
   DO i = 0, npx-1
      nxvec(i) = bnc(nxt,npx,i)
   END DO
   DO i = 0, npy-1
      nyvec(i) = bnc(nyt,npy,i)
   END DO
   DO i = 0, npz-1
      nzvec(i) = bnc(nzt,npz,i)
   END DO
END IF

! Initialize ipak
ipak = 0

! Send to each process their size (nx,ny,nz) so that they can allocate space
IF (irank == root) THEN
   ! Send the root's nrank to everyone so they know the message source in allcomm
   sndrnk = nrank
   CALL MPI_BCAST(sndrnk,1,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
   DO pk = 0, npz-1
      coord(1) = pk
      ipak(3) = nzvec(pk)
      tnz = ipak(3)
      DO pj = 0, npy-1
         coord(2) = pj
         ipak(2) = nyvec(pj)
         tny = ipak(2)
         DO pi = 0, npx-1
            coord(3) = pi
            ipak(1) = nxvec(pi)
            tnx = ipak(1)
            ! Get the rank according to alcomm and set up tags
            CALL MPI_CART_RANK(allcomm,coord,recrnk,ierr)
            tag = 100 + recrnk
            IF (recrnk /= nrank) THEN
               ! Send the dimensions to the appropriate process in allcomm
               CALL MPI_SEND(ipak,3,MPI_INTEGER,recrnk,tag,allcomm,ierr)
            ELSE
               ! Set dimensions and allocate arays for later use and initialize
               nx = ipak(1)
               ny = ipak(2)
               nz = ipak(3)
               bcs = apo*(nx*ny+nx*nz+ny*nz)
               ALLOCATE(dx(nx), dy(ny), dz(nz), mat(nx,ny,nz), s(nx,ny,nz,ng), fig(nx,ny,nz))
               ALLOCATE(psii(bcs,8))
               dx = 0.0
               dy = 0.0
               dz = 0.0
               mat = 0
               s = 0.0
               fig = 0.0
               psii = 0.0
            END IF
         END DO
      END DO
   END DO
ELSE
   ! Get the root's sending rank for allcomm
   CALL MPI_BCAST(sndrnk,1,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
   ! Ranks that are not root need to receive data and allocate arrays
   tag = 100 + nrank
   CALL MPI_RECV(ipak,3,MPI_INTEGER,sndrnk,tag,allcomm,istat,ierr)
   nx = ipak(1)
   ny = ipak(2)
   nz = ipak(3)
   bcs = apo*(nx*ny+nx*nz+ny*nz)
   ALLOCATE(dx(nx), dy(ny), dz(nz), mat(nx,ny,nz), s(nx,ny,nz,ng), fig(nx,ny,nz))
   ALLOCATE(psii(bcs,8))
   dx = 0.0
   dy = 0.0
   dz = 0.0
   mat = 0
   s = 0.0
   fig = 0.0
   psii = 0.0
END IF

! Root (from MPI_COMM_WORLD) has full data
! Loop over the grid and send data to the processes: dx, dy, dz, mat, s, psiin
IF (irank == root) THEN
   ! Organize the data to be sent, get the rank on the receiving end, send
   ! Use dummy index to know how much to send and where to start the buffer
   kndx = 0
   DO pk = 0, npz-1
      coord(1) = pk
      tnz = nzvec(pk)
      jndx = 0
      DO pj = 0, npy-1
         coord(2) = pj
         tny = nyvec(pj)
         indx = 0
         DO pi = 0, npx-1
            coord(3) = pi
            tnx = nxvec(pi)

            ! Allocate the matrices
            bcs = apo*(tnx*tny+tnx*tnz+tny*tnz)
            ALLOCATE(tempmat(tnx,tny,tnz), temps(tnx,tny,tnz,ng), tempf(tnx,tny,tnz))
            ALLOCATE(tempsi(bcs,8))
            tempmat(1:tnx,1:tny,1:tnz) = matt((indx+1):(indx+tnx),(jndx+1):(jndx+tny),(kndx+1):(kndx+tnz))
            temps(1:tnx,1:tny,1:tnz,:) = st((indx+1):(indx+tnx),(jndx+1):(jndx+tny),(kndx+1):(kndx+tnz),:)
            tempf(1:tnx,1:tny,1:tnz) = flux((indx+1):(indx+tnx),(jndx+1):(jndx+tny),(kndx+1):(kndx+tnz),1)
            
            ! Include interpolated error (copy)
            tempf(1:tnx,1:tny,1:tnz) = tempf(1:tnx,1:tny,1:tnz) + cerr(pi+1,pj+1,pk+1)

            ! Initialize the angular flux BCs, then write over with actual data
            tempsi = 0.0
            IF (pk == 0) THEN
               DO n = 1, apo
                  indx1 = (n-1)*(tnx*tny)
                  DO j = 1, tny
                     DO i = 1, tnx
                        ieq = indx1 + (j-1)*tnx + i
                        tempsi(ieq,1) = afz(indx+i,jndx+j,n,1) + zerr(pi+1,pj+1,n,1)
                        tempsi(ieq,2) = afz(indx+i,jndx+j,n,2) + zerr(pi+1,pj+1,n,2)
                        tempsi(ieq,3) = afz(indx+i,jndx+j,n,3) + zerr(pi+1,pj+1,n,3)
                        tempsi(ieq,4) = afz(indx+i,jndx+j,n,4) + zerr(pi+1,pj+1,n,4)
                     END DO
                  END DO
               END DO
            END IF
            IF (pk == npz-1) THEN
               DO n = 1, apo
                  indx1 = (n-1)*(tnx*tny)
                  DO j = 1, tny
                     DO i = 1, tnx
                        ieq = indx1 + (j-1)*tnx + i
                        tempsi(ieq,5) = afz(indx+i,jndx+j,n,5) + zerr(pi+1,pj+1,n,5)
                        tempsi(ieq,6) = afz(indx+i,jndx+j,n,6) + zerr(pi+1,pj+1,n,6)
                        tempsi(ieq,7) = afz(indx+i,jndx+j,n,7) + zerr(pi+1,pj+1,n,7)
                        tempsi(ieq,8) = afz(indx+i,jndx+j,n,8) + zerr(pi+1,pj+1,n,8)
                     END DO
                  END DO
               END DO
            END IF
            IF (pj == 0) THEN
               DO n = 1, apo
                  indx1 = apo*(tnx*tny) + (n-1)*(tnx*tnz)
                  DO j = 1, tnz
                     DO i = 1, tnx
                        ieq = indx1 + (j-1)*tnx + i
                        tempsi(ieq,1) = afy(indx+i,kndx+j,n,1) + yerr(pi+1,pk+1,n,1)
                        tempsi(ieq,2) = afy(indx+i,kndx+j,n,2) + yerr(pi+1,pk+1,n,2)
                        tempsi(ieq,5) = afy(indx+i,kndx+j,n,5) + yerr(pi+1,pk+1,n,5)
                        tempsi(ieq,6) = afy(indx+i,kndx+j,n,6) + yerr(pi+1,pk+1,n,6)
                     END DO
                  END DO
               END DO
            END IF
            IF (pj == npy-1) THEN
               DO n = 1, apo
                  indx1 = apo*(tnx*tny) + (n-1)*(tnx*tnz)
                  DO j = 1, tnz
                     DO i = 1, tnx
                        ieq = indx1 + (j-1)*tnx + i
                        tempsi(ieq,3) = afy(indx+i,kndx+j,n,3) + yerr(pi+1,pk+1,n,3)
                        tempsi(ieq,4) = afy(indx+i,kndx+j,n,4) + yerr(pi+1,pk+1,n,4)
                        tempsi(ieq,7) = afy(indx+i,kndx+j,n,7) + yerr(pi+1,pk+1,n,7)
                        tempsi(ieq,8) = afy(indx+i,kndx+j,n,8) + yerr(pi+1,pk+1,n,8)
                     END DO
                  END DO
               END DO
            END IF
            IF (pi == 0) THEN
               DO n = 1, apo
                  indx1 = apo*(tnx*tny+tnx*tnz) + (n-1)*(tny*tnz)
                  DO j = 1, tnz
                     DO i = 1, tny
                        ieq = indx1 + (j-1)*tny + i
                        tempsi(ieq,1) = afx(jndx+i,kndx+j,n,1) + xerr(pj+1,pk+1,n,1)
                        tempsi(ieq,4) = afx(jndx+i,kndx+j,n,4) + xerr(pj+1,pk+1,n,4)
                        tempsi(ieq,5) = afx(jndx+i,kndx+j,n,5) + xerr(pj+1,pk+1,n,5)
                        tempsi(ieq,8) = afx(jndx+i,kndx+j,n,8) + xerr(pj+1,pk+1,n,8)
                     END DO
                  END DO
               END DO
            END IF
            IF (pi == npx-1) THEN
               DO n = 1, apo
                  indx1 = apo*(tnx*tny+tnx*tnz) + (n-1)*(tny*tnz)
                  DO j = 1, tnz
                     DO i = 1, tny
                        ieq = indx1 + (j-1)*tny + i
                        tempsi(ieq,2) = afx(jndx+i,kndx+j,n,2) + xerr(pj+1,pk+1,n,2)
                        tempsi(ieq,3) = afx(jndx+i,kndx+j,n,3) + xerr(pj+1,pk+1,n,3)
                        tempsi(ieq,6) = afx(jndx+i,kndx+j,n,6) + xerr(pj+1,pk+1,n,6)
                        tempsi(ieq,7) = afx(jndx+i,kndx+j,n,7) + xerr(pj+1,pk+1,n,7)
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
            tagf = 800 + recrnk
            ! Send if not on the root's process, just copy if it is
            IF (recrnk /= sndrnk) THEN
               CALL MPI_SEND(dxt(indx+1),tnx,MPI_DOUBLE_PRECISION,recrnk,tagi,allcomm,ierr)
               CALL MPI_SEND(dyt(jndx+1),tny,MPI_DOUBLE_PRECISION,recrnk,tagj,allcomm,ierr)
               CALL MPI_SEND(dzt(kndx+1),tnz,MPI_DOUBLE_PRECISION,recrnk,tagk,allcomm,ierr)
               CALL MPI_SEND(tempmat,tnx*tny*tnz,MPI_INTEGER,recrnk,tagm,allcomm,ierr)
               CALL MPI_SEND(temps,tnx*tny*tnz*ng,MPI_DOUBLE_PRECISION,recrnk,tags,allcomm,ierr)
               CALL MPI_SEND(tempf,tnx*tny*tnz,MPI_DOUBLE_PRECISION,recrnk,tagf,allcomm,ierr)
               CALL MPI_SEND(tempsi,bcs*8,MPI_DOUBLE_PRECISION,recrnk,tagp,allcomm,ierr)
            ELSE
               dx(1:nx) = dxt((indx+1):(indx+tnx))
               dy(1:ny) = dyt((jndx+1):(jndx+tny))
               dz(1:nz) = dzt((kndx+1):(kndx+tnz))
               mat(1:nx,1:ny,1:nz) = tempmat(1:tnx,1:tny,1:tnz)
               s(1:nx,1:ny,1:nz,:) = temps(1:tnx,1:tny,1:tnz,:)
               fig(1:nx,1:ny,1:nz) = tempf(1:tnx,1:tny,1:tnz)
               psii(1:(apo*(nx*ny+nx*nz+ny*nz)),1:8) = tempsi(1:bcs,1:8)
            END IF

            ! Update indx before cycling pi
            indx = indx + tnx

            ! Deallocate the matrices so they can change size afterwards
            DEALLOCATE(tempmat, temps, tempf, tempsi)
         END DO

         ! Update jndx before cycling pj
         jndx = jndx + tny
      END DO

      ! Update kndx before cycling
      kndx = kndx + tnz
   END DO
ELSE
   ! Other processes receive their message and copy into their variables accordingly
   tagi = 200 + nrank
   tagj = 300 + nrank
   tagk = 400 + nrank
   tagm = 500 + nrank
   tags = 600 + nrank
   tagp = 700 + nrank
   tagf = 800 + nrank
   CALL MPI_RECV(dx,nx,MPI_DOUBLE_PRECISION,sndrnk,tagi,allcomm,istat,ierr)
   CALL MPI_RECV(dy,ny,MPI_DOUBLE_PRECISION,sndrnk,tagj,allcomm,istat,ierr)
   CALL MPI_RECV(dz,nz,MPI_DOUBLE_PRECISION,sndrnk,tagk,allcomm,istat,ierr)
   CALL MPI_RECV(mat,nx*ny*nz,MPI_INTEGER,sndrnk,tagm,allcomm,istat,ierr)
   CALL MPI_RECV(s,nx*ny*nz*ng,MPI_DOUBLE_PRECISION,sndrnk,tags,allcomm,istat,ierr)
   CALL MPI_RECV(fig,nx*ny*nz,MPI_DOUBLE_PRECISION,sndrnk,tagf,allcomm,istat,ierr)
   CALL MPI_RECV(psii,bcs*8,MPI_DOUBLE_PRECISION,sndrnk,tagp,allcomm,istat,ierr)
END IF

RETURN
END SUBROUTINE psds 
