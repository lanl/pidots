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
INTEGER :: i, j, k, n, pi, pj, pk, recrnk, sndrnk, fneq
INTEGER :: ii, jj, kk, isd, jsd, ksd, sdi
INTEGER :: kndx, jndx, indx, indx1, indx2, ieq, ierr
INTEGER :: tag, tagi, tagj, tagk, tagm, tags, tagp
INTEGER, DIMENSION(3) :: coord, istat
INTEGER, DIMENSION(nx,ny,nz) :: tempmat
REAL*8, DIMENSION(nx,ny,nz) :: temps
!REAL*8, DIMENSION(bcs,8,nsdp) :: tempsi

INCLUDE 'mpif.h'

! Call for the processgrid to be created
CALL procgrid

fneq = nx*ny*nz

! Allocate array sizes from number of cells on each processor and initialize
ALLOCATE(dx(nx), dy(ny), dz(nz), mat(nx,ny,nz), s(nx,ny,nz))
ALLOCATE(psii(bcs,8,nsdp))
dx = 0.0d0
dy = 0.0d0
dz = 0.0d0
mat = 0
s = 0.0d0
psii = 0.0d0

! Send the root's nrank to everyone so they know the message source in allcomm
IF (irank == root) sndrnk = nrank
CALL MPI_BCAST(sndrnk,1,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)

! Root (from MPI_COMM_WORLD) has full data
! Loop over the grid and send data to the processes: dx, dy, dz, mat, s, psiin
IF (irank == root) THEN
   ! Organize the data to be sent, get the rank on the receiving end, send
   DO pk = 0, npz-1
      coord(1) = pk
      kndx = pk*nz
      DO pj = 0, npy-1
         coord(2) = pj
         jndx = pj*ny
         DO pi = 0, npx-1
            coord(3) = pi
            indx = pi*nx

            tempmat = matt((indx+1):(indx+nx),(jndx+1):(jndx+ny),(kndx+1):(kndx+nz))
            temps = st((indx+1):(indx+nx),(jndx+1):(jndx+ny),(kndx+1):(kndx+nz))

!            ! Initialize the angular flux BCs, then write over with actual data
!            tempsi = 0.0
!            IF (pk == 0 .AND. bc(1) == 2) THEN
!               DO n = 1, apo
!                  indx1 = (n-1)*xys
!                  DO j = 1, ny
!                     jsd = (j-1)/sdny + 1
!                     jj = MOD((j-1),sdny) + 1
!                     DO i = 1, nx
!                        isd = (i-1)/sdnx + 1
!                        ii = MOD((i-1),sdnx) + 1
!                        ieq = indx1 + (jj-1)*sdnx + ii
!                        sdi = (jsd-1)*nxsd + isd
!                        tempsi(ieq,1,sdi) = posinz(indx+i,jndx+j,n,1)
!                        tempsi(ieq,2,sdi) = posinz(indx+i,jndx+j,n,2)
!                        tempsi(ieq,3,sdi) = posinz(indx+i,jndx+j,n,3)
!                        tempsi(ieq,4,sdi) = posinz(indx+i,jndx+j,n,4)
!                     END DO
!                  END DO
!               END DO
!            END IF
!            IF (pk == npz-1 .AND. bc(2) == 2) THEN
!               indx2 = (nzsd-1)*nysd*nxsd
!               DO n = 1, apo
!                  indx1 = (n-1)*xys
!                  DO j = 1, ny
!                     jsd = (j-1)/sdny + 1
!                     jj = MOD((j-1),sdny) + 1
!                     DO i = 1, nx
!                        isd = (i-1)/sdnx + 1
!                        ii = MOD((i-1),sdnx) + 1
!                        ieq = indx1 + (jj-1)*sdnx + ii
!                        sdi = indx2 + (jsd-1)*nxsd + isd
!                        tempsi(ieq,5,sdi) = neginz(indx+i,jndx+j,n,1)
!                        tempsi(ieq,6,sdi) = neginz(indx+i,jndx+j,n,2)
!                        tempsi(ieq,7,sdi) = neginz(indx+i,jndx+j,n,3)
!                        tempsi(ieq,8,sdi) = neginz(indx+i,jndx+j,n,4)
!                     END DO
!                  END DO
!               END DO
!            END IF
!            IF (pj == 0 .AND. bc(3) == 2) THEN
!               DO n = 1, apo
!                  indx1 = apo*xys + (n-1)*xzs
!                  DO k = 1, nz
!                     ksd = (k-1)/sdnz + 1
!                     kk = MOD((k-1),sdnz) + 1
!                     DO i = 1, nx
!                        isd = (i-1)/sdnx + 1
!                        ii = MOD((i-1),sdnx) + 1
!                        ieq = indx1 + (kk-1)*sdnx + ii
!                        sdi = (ksd-1)*nysd*nxsd + isd
!                        tempsi(ieq,1,sdi) = posiny(indx+i,kndx+k,n,1)
!                        tempsi(ieq,2,sdi) = posiny(indx+i,kndx+k,n,2)
!                        tempsi(ieq,5,sdi) = posiny(indx+i,kndx+k,n,3)
!                        tempsi(ieq,6,sdi) = posiny(indx+i,kndx+k,n,4)
!                     END DO
!                  END DO
!               END DO
!            END IF
!            IF (pj == npy-1 .AND. bc(4) == 2) THEN
!               indx2 = (nysd-1)*nxsd
!               DO n = 1, apo
!                  indx1 = apo*xys + (n-1)*xzs
!                  DO k = 1, nz
!                     ksd = (k-1)/sdnz + 1
!                     kk = MOD((k-1),sdnz) + 1
!                     DO i = 1, nx
!                        isd = (i-1)/sdnx + 1
!                        ii = MOD((i-1),sdnx) + 1
!                        ieq = indx1 + (kk-1)*sdnx + ii
!                        sdi = indx2 + (ksd-1)*nysd*nxsd + isd
!                        tempsi(ieq,3,sdi) = neginy(indx+i,kndx+k,n,1)
!                        tempsi(ieq,4,sdi) = neginy(indx+i,kndx+k,n,2)
!                        tempsi(ieq,7,sdi) = neginy(indx+i,kndx+k,n,3)
!                        tempsi(ieq,8,sdi) = neginy(indx+i,kndx+k,n,4)
!                     END DO
!                  END DO
!               END DO
!            END IF
!            IF (pi == 0 .AND. bc(5) == 2) THEN
!               DO n = 1, apo
!                  indx1 = apo*(xys+xzs) + (n-1)*yzs
!                  DO k = 1, nz
!                     ksd = (k-1)/sdnz + 1
!                     kk = MOD((k-1),sdnz) + 1
!                     DO j = 1, ny
!                        jsd = (j-1)/sdny + 1
!                        jj = MOD((j-1),sdny) + 1
!                        ieq = indx1 + (kk-1)*sdny + jj
!                        sdi = (ksd-1)*nysd*nxsd + (jsd-1)*nxsd + 1
!                        tempsi(ieq,1,sdi) = posinx(jndx+j,kndx+k,n,1)
!                        tempsi(ieq,4,sdi) = posinx(jndx+j,kndx+k,n,2)
!                        tempsi(ieq,5,sdi) = posinx(jndx+j,kndx+k,n,3)
!                        tempsi(ieq,8,sdi) = posinx(jndx+j,kndx+k,n,4)
!                     END DO
!                  END DO
!               END DO
!            END IF
!            IF (pi == npx-1 .AND. bc(6) == 2) THEN
!               DO n = 1, apo
!                  indx1 = apo*(xys+xzs) + (n-1)*yzs
!                  DO k = 1, nz
!                     ksd = (k-1)/sdnz + 1
!                     kk = MOD((k-1),sdnz) + 1
!                     DO j = 1, ny
!                        jsd = (j-1)/sdny + 1
!                        jj = MOD((j-1),sdny) + 1
!                        ieq = indx1 + (kk-1)*sdny + jj
!                        sdi = (ksd-1)*nysd*nxsd + (jsd-1)*nxsd + nxsd
!                        tempsi(ieq,2,sdi) = neginx(jndx+j,kndx+k,n,1)
!                        tempsi(ieq,3,sdi) = neginx(jndx+j,kndx+k,n,2)
!                        tempsi(ieq,6,sdi) = neginx(jndx+j,kndx+k,n,3)
!                        tempsi(ieq,7,sdi) = neginx(jndx+j,kndx+k,n,4)
!                     END DO
!                  END DO
!               END DO
!            END IF

            ! Get the rank to send data to and set up the tags
            CALL MPI_CART_RANK(allcomm,coord,recrnk,ierr)
            tagi = 200 + recrnk
            tagj = 300 + recrnk
            tagk = 400 + recrnk
            tagm = 500 + recrnk
            tags = 600 + recrnk
!            tagp = 700 + recrnk
            ! Send if not on the root's process, just copy if it is
            IF (recrnk /= sndrnk) THEN
               CALL MPI_SEND(dxt(indx+1),nx,MPI_DOUBLE_PRECISION,recrnk,tagi,allcomm,ierr)
               CALL MPI_SEND(dyt(jndx+1),ny,MPI_DOUBLE_PRECISION,recrnk,tagj,allcomm,ierr)
               CALL MPI_SEND(dzt(kndx+1),nz,MPI_DOUBLE_PRECISION,recrnk,tagk,allcomm,ierr)
               CALL MPI_SEND(tempmat,fneq,MPI_INTEGER,recrnk,tagm,allcomm,ierr)
               CALL MPI_SEND(temps,fneq,MPI_DOUBLE_PRECISION,recrnk,tags,allcomm,ierr)
!               CALL MPI_SEND(tempsi,bcs*8*nsdp,MPI_DOUBLE_PRECISION,recrnk,tagp,allcomm,ierr)
            ELSE
               dx = dxt((indx+1):(indx+nx))
               dy = dyt((jndx+1):(jndx+ny))
               dz = dzt((kndx+1):(kndx+nz))
               mat = tempmat
               s = temps
!               psii = tempsi
            END IF
         END DO
      END DO
   END DO
ELSE
   ! Other processes receive their message and copy into their variables accordingly
   tagi = 200 + nrank
   tagj = 300 + nrank
   tagk = 400 + nrank
   tagm = 500 + nrank
   tags = 600 + nrank
!   tagp = 700 + nrank
   CALL MPI_RECV(dx,nx,MPI_DOUBLE_PRECISION,sndrnk,tagi,allcomm,istat,ierr)
   CALL MPI_RECV(dy,ny,MPI_DOUBLE_PRECISION,sndrnk,tagj,allcomm,istat,ierr)
   CALL MPI_RECV(dz,nz,MPI_DOUBLE_PRECISION,sndrnk,tagk,allcomm,istat,ierr)
   CALL MPI_RECV(mat,fneq,MPI_INTEGER,sndrnk,tagm,allcomm,istat,ierr)
   CALL MPI_RECV(s,fneq,MPI_DOUBLE_PRECISION,sndrnk,tags,allcomm,istat,ierr)
!   CALL MPI_RECV(psii,bcs*8*nsdp,MPI_DOUBLE_PRECISION,sndrnk,tagp,allcomm,istat,ierr)
END IF

RETURN
END SUBROUTINE psds 
