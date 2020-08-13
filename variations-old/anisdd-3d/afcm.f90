SUBROUTINE afcm(i,j,k,xs,ys,zs,xe,ye,ze,incx,incy,incz,n,oct,mu,omega)

!-------------------------------------------------------------
!
!  Angular Flux Coefficient Matrix
!  Directs the construction of the matrices for working with
!   the incoming angular flux. Constructs kmat's and kpsi's
!
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER, INTENT(IN) :: i, j, k, xs, ys, zs, xe, ye, ze, incx, incy, incz, n, oct
INTEGER :: ii, jj, kk, ieq, jeq, indx
REAL*8, INTENT(IN) :: mu, omega
REAL*8, DIMENSION(:,:), ALLOCATABLE :: ktmpz, ktmpy, ktmpx

! Save the BC matrices for updates
IF (i /= xs) THEN
   zbcyzo = zbcyz
   ybcyzo = ybcyz
   xbcyzo = xbcyz
END IF
IF (j /= ys) THEN
   zbcxzo = zbcxz(:,:,i)
   ybcxzo = ybcxz(:,:,i)
   xbcxzo = xbcxz(:,:,i)
END IF
IF (k /= zs) THEN
   zbcxyo = zbcxy(:,:,i,j)
   ybcxyo = ybcxy(:,:,i,j)
   xbcxyo = xbcxy(:,:,i,j)
END IF

! Earlier versions that did not give ang. flux out did not solve at xe, ye, ze
! Now solve at xe, ye, ze to use for computing angular flux out

! Update BC matrices with xbcyzo (BC data in x-direction into cell)
IF (i /= xs) THEN
   DO kk = zs, k, incz
      DO jj = ys, j, incy
         xbcyz(jj,kk) = gyzyz*xbcyzo(jj,kk)
         xbcxz(jj,kk,i) = gxzyz*xbcyzo(jj,kk)
         xbcxy(jj,kk,i,j) = gxyyz*xbcyzo(jj,kk)
      END DO
   END DO
END IF

! Update BC matrices with ybcyzo
IF (i /= xs) THEN
   DO kk = zs, k, incz
      DO ii = xs, i-incx, incx
         ybcyz(ii,kk) = gyzyz*ybcyzo(ii,kk)
         ybcxz(ii,kk,i) = gxzyz*ybcyzo(ii,kk)
         ybcxy(ii,kk,i,j) = gxyyz*ybcyzo(ii,kk)
      END DO
   END DO
END IF

! Update BC matrices with zbcyzo
IF (i /= xs) THEN
   DO jj = ys, j, incy
      DO ii = xs, i-incx, incx
         zbcyz(ii,jj) = gyzyz*zbcyzo(ii,jj)
         zbcxz(ii,jj,i) = gxzyz*zbcyzo(ii,jj)
         zbcxy(ii,jj,i,j) = gxyyz*zbcyzo(ii,jj)
      END DO
   END DO
END IF

! Update BC matrices with xbcxzo (BC data in y-dir into cell)
IF (j /= ys) THEN
   DO kk = zs, k, incz
      DO jj = ys, j-incy, incy
         IF (i /= xs) THEN
            xbcyz(jj,kk) = xbcyz(jj,kk) + gyzxz*xbcxzo(jj,kk)
            xbcxz(jj,kk,i) = xbcxz(jj,kk,i) + gxzxz*xbcxzo(jj,kk)
            xbcxy(jj,kk,i,j) = xbcxy(jj,kk,i,j) + gxyxz*xbcxzo(jj,kk)
         ELSE
            xbcyz(jj,kk) = gyzxz*xbcxzo(jj,kk)
            xbcxz(jj,kk,i) = gxzxz*xbcxzo(jj,kk)
            xbcxy(jj,kk,i,j) = gxyxz*xbcxzo(jj,kk)
         END IF
      END DO
   END DO
END IF

! Update BC matrices with ybcxzo
IF (j /= ys) THEN
   DO kk = zs, k, incz
      DO ii = xs, i, incx
         IF (ii /= i) THEN
            ybcyz(ii,kk) = ybcyz(ii,kk) + gyzxz*ybcxzo(ii,kk)
            ybcxz(ii,kk,i) = ybcxz(ii,kk,i) + gxzxz*ybcxzo(ii,kk)
            ybcxy(ii,kk,i,j) = ybcxy(ii,kk,i,j) + gxyxz*ybcxzo(ii,kk)
         ELSE
            ybcyz(ii,kk) = gyzxz*ybcxzo(ii,kk)
            ybcxz(ii,kk,i) = gxzxz*ybcxzo(ii,kk)
            ybcxy(ii,kk,i,j) = gxyxz*ybcxzo(ii,kk)
         END IF
      END DO
   END DO
END IF

! Update BC matrices with zbcxzo
IF (j /= ys) THEN
   DO jj = ys, j-incy, incy
      DO ii = xs, i, incx
         IF (ii /= i) THEN
            zbcyz(ii,jj) = zbcyz(ii,jj) + gyzxz*zbcxzo(ii,jj)
            zbcxz(ii,jj,i) = zbcxz(ii,jj,i) + gxzxz*zbcxzo(ii,jj)
            zbcxy(ii,jj,i,j) = zbcxy(ii,jj,i,j) + gxyxz*zbcxzo(ii,jj)
         ELSE
            zbcyz(ii,jj) = gyzxz*zbcxzo(ii,jj)
            zbcxz(ii,jj,i) = gxzxz*zbcxzo(ii,jj)
            zbcxy(ii,jj,i,j) = gxyxz*zbcxzo(ii,jj)
         END IF
      END DO
   END DO
END IF

! Update BC matrices xbcxyo (BC data in z-dir into cell)
IF (k /= zs) THEN
   DO kk = zs, k-incz, incz
      DO jj = ys, j, incy
         IF (jj /= j .OR. i /= xs) THEN
            xbcyz(jj,kk) = xbcyz(jj,kk) + gyzxy*xbcxyo(jj,kk)
            xbcxz(jj,kk,i) = xbcxz(jj,kk,i) + gxzxy*xbcxyo(jj,kk)
            xbcxy(jj,kk,i,j) = xbcxy(jj,kk,i,j) + gxyxy*xbcxyo(jj,kk)
         ELSE IF (jj == j .AND. i == xs) THEN
            xbcyz(jj,kk) = gyzxy*xbcxyo(jj,kk)
            xbcxz(jj,kk,i) = gxzxy*xbcxyo(jj,kk)
            xbcxy(jj,kk,i,j) = gxyxy*xbcxyo(jj,kk)
         END IF
      END DO
   END DO
END IF

! Update BC matrices ybcxyo
IF (k /= zs) THEN
   DO kk = zs, k-incz, incz
      DO ii = xs, i, incx
         IF (ii /= i .OR. j /= ys) THEN
            ybcyz(ii,kk) = ybcyz(ii,kk) + gyzxy*ybcxyo(ii,kk)
            ybcxz(ii,kk,i) = ybcxz(ii,kk,i) + gxzxy*ybcxyo(ii,kk)
            ybcxy(ii,kk,i,j) = ybcxy(ii,kk,i,j) + gxyxy*ybcxyo(ii,kk)
         ELSE IF (ii == i .AND. j == ys) THEN
            ybcyz(ii,kk) = gyzxy*ybcxyo(ii,kk)
            ybcxz(ii,kk,i) = gxzxy*ybcxyo(ii,kk)
            ybcxy(ii,kk,i,j) = gxyxy*ybcxyo(ii,kk)
         END IF
      END DO
   END DO
END IF

! Update BC matrices zbcxyo
IF (k /= zs) THEN
   DO jj = ys, j, incy
      DO ii = xs, i, incx
         IF (ii /= i .OR. jj /= j) THEN
            zbcyz(ii,jj) = zbcyz(ii,jj) + gyzxy*zbcxyo(ii,jj)
            zbcxz(ii,jj,i) = zbcxz(ii,jj,i) + gxzxy*zbcxyo(ii,jj)
            zbcxy(ii,jj,i,j) = zbcxy(ii,jj,i,j) + gxyxy*zbcxyo(ii,jj)
         ELSE IF (ii == i .AND. jj == j) THEN
            zbcyz(ii,jj) = gyzxy*zbcxyo(ii,jj)
            zbcxz(ii,jj,i) = gxzxy*zbcxyo(ii,jj)
            zbcxy(ii,jj,i,j) = gxyxy*zbcxyo(ii,jj)
         END IF
      END DO
   END DO
END IF

! Adjust the values for BCs to their adjacent cells
IF (i == xs) THEN
   xbcyz(j,k) = gyzyz
   xbcxz(j,k,i) = gxzyz
   xbcxy(j,k,i,j) = gxyyz
END IF
IF (j == ys) THEN
   ybcyz(i,k) = gyzxz
   ybcxz(i,k,i) = gxzxz
   ybcxy(i,k,i,j) = gxyxz
END IF
IF (k == zs) THEN
   zbcyz(i,j) = gyzxy
   zbcxz(i,j,i) = gxzxy
   zbcxy(i,j,i,j) = gxyxy
END IF

!-----------------------------------------------------------------------------
! Update the appropriate kmat given the octant of interest
! Call conk to construct the row of kmat, then set it in kmat#
! Set according to 'tpose' flag
IF (tpose == 1) THEN
   ALLOCATE(ktmpz(xys,nmom), ktmpy(xzs,nmom), ktmpx(yzs,nmom))
   jeq = ((i-1) + (j-1)*nx + (k-1)*nx*ny)*nmom
   CALL conkt(n,i,j,k,incx,incy,incz,xs,ys,zs,mu,omega,ktmpz,ktmpy,ktmpx)
   indx = (n-1)*nx*ny
   kmat(indx+1:indx+xys,jeq+1:jeq+nmom,oct) = kmat(indx+1:indx+xys,jeq+1:jeq+nmom,oct) + ktmpz
   indx = apo*nx*ny + (n-1)*nx*nz
   kmat(indx+1:indx+xzs,jeq+1:jeq+nmom,oct) = kmat(indx+1:indx+xzs,jeq+1:jeq+nmom,oct) + ktmpy
   indx = apo*(nx*ny+nx*nz) + (n-1)*ny*nz
   kmat(indx+1:indx+yzs,jeq+1:jeq+nmom,oct) = kmat(indx+1:indx+yzs,jeq+1:jeq+nmom,oct) + ktmpx

   ! Use the data in the xbc, ybc, zbc varialbes at xe ye ze to compute kpsi
   ! Start with z-direction coefficients
   IF (k == ze) THEN
      ! Reset equation index
      jeq = i + (j-1)*nx
      DO jj = ys, j, incy
         DO ii = xs, i, incx
            ieq = ii + (jj-1)*nx
            kpsi(ieq,jeq,n,oct) = zbcxy(ii,jj,i,j)
         END DO
      END DO
      DO kk = zs, k, incz
         DO ii = xs, i, incx
            ieq = (nx*ny) + ii + (kk-1)*nx
            kpsi(ieq,jeq,n,oct) = ybcxy(ii,kk,i,j)
         END DO
      END DO
      DO kk = zs, k, incz
         DO jj = ys, j, incy
            ieq = (nx*ny+nx*nz) + jj + (kk-1)*ny
            kpsi(ieq,jeq,n,oct) = xbcxy(jj,kk,i,j)
         END DO
      END DO
   END IF
   ! Next get the y-direction coefficients
   IF (j == ye) THEN
      ! Reset equation index
      jeq = (nx*ny) + i + (k-1)*nx
      DO jj = ys, j, incy
         DO ii = xs, i, incx
            ieq = ii + (jj-1)*nx
            kpsi(ieq,jeq,n,oct) = zbcxz(ii,jj,i)
         END DO
      END DO
      DO kk = zs, k, incz
         DO ii = xs, i, incx
            ieq = (nx*ny) + ii + (kk-1)*nx
            kpsi(ieq,jeq,n,oct) = ybcxz(ii,kk,i)
         END DO
      END DO
      DO kk = zs, k, incz
         DO jj = ys, j, incy
            ieq = (nx*ny+nx*nz) + jj + (kk-1)*ny
            kpsi(ieq,jeq,n,oct) = xbcxz(jj,kk,i)
         END DO
      END DO
   END IF
   ! Next get the x-direction coefficients
   IF (i == xe) THEN
      ! Reset equation index
      jeq = (nx*ny+nx*nz) + j + (k-1)*ny
      DO jj = ys, j, incy
         DO ii = xs, i, incx
            ieq = ii + (jj-1)*nx
            kpsi(ieq,jeq,n,oct) = zbcyz(ii,jj)
         END DO
      END DO
      DO kk = zs, k, incz
         DO ii = xs, i, incx
            ieq = (nx*ny) + ii + (kk-1)*nx
            kpsi(ieq,jeq,n,oct) = ybcyz(ii,kk)
         END DO
      END DO
      DO kk = zs, k, incz
         DO jj = ys, j, incy
            ieq = (nx*ny+nx*nz) + jj + (kk-1)*ny
            kpsi(ieq,jeq,n,oct) = xbcyz(jj,kk)
         END DO
      END DO
   END IF
!-------------------------------------------------------------------------------
ELSE
   ALLOCATE(ktmpz(nmom,xys), ktmpy(nmom,xzs), ktmpx(nmom,yzs))
   ieq = ((i-1) + (j-1)*nx + (k-1)*nx*ny)*nmom
   CALL conk(n,i,j,k,incx,incy,incz,xs,ys,zs,mu,omega,ktmpz,ktmpy,ktmpx)
   indx = (n-1)*nx*ny
   kmat(ieq+1:ieq+nmom,indx+1:indx+xys,oct) = kmat(ieq+1:ieq+nmom,indx+1:indx+xys,oct) + ktmpz
   indx = apo*nx*ny + (n-1)*nx*nz
   kmat(ieq+1:ieq+nmom,indx+1:indx+xzs,oct) = kmat(ieq+1:ieq+nmom,indx+1:indx+xzs,oct) + ktmpy
   indx = apo*(nx*ny+nx*nz) + (n-1)*ny*nz
   kmat(ieq+1:ieq+nmom,indx+1:indx+yzs,oct) = kmat(ieq+1:ieq+nmom,indx+1:indx+yzs,oct) + ktmpx

   ! Use the data in the xbc, ybc, zbc varialbes at xe ye ze to compute kpsi
   ! Start with z-direction coefficients
   IF (k == ze) THEN
      ! Reset equation index
      ieq = i + (j-1)*nx
      DO jj = ys, j, incy
         DO ii = xs, i, incx
            jeq = ii + (jj-1)*nx
            kpsi(ieq,jeq,n,oct) = zbcxy(ii,jj,i,j)
         END DO
      END DO
      DO kk = zs, k, incz
         DO ii = xs, i, incx
            jeq = (nx*ny) + ii + (kk-1)*nx
            kpsi(ieq,jeq,n,oct) = ybcxy(ii,kk,i,j)
         END DO
      END DO
      DO kk = zs, k, incz
         DO jj = ys, j, incy
            jeq = (nx*ny+nx*nz) + jj + (kk-1)*ny
            kpsi(ieq,jeq,n,oct) = xbcxy(jj,kk,i,j)
         END DO
      END DO
   END IF
   ! Next get the y-direction coefficients
   IF (j == ye) THEN
      ! Reset equation index
      ieq = (nx*ny) + i + (k-1)*nx
      DO jj = ys, j, incy
         DO ii = xs, i, incx
            jeq = ii + (jj-1)*nx
            kpsi(ieq,jeq,n,oct) = zbcxz(ii,jj,i)
         END DO
      END DO
      DO kk = zs, k, incz
         DO ii = xs, i, incx
            jeq = (nx*ny) + ii + (kk-1)*nx
            kpsi(ieq,jeq,n,oct) = ybcxz(ii,kk,i)
         END DO
      END DO
      DO kk = zs, k, incz
         DO jj = ys, j, incy
            jeq = (nx*ny+nx*nz) + jj + (kk-1)*ny
            kpsi(ieq,jeq,n,oct) = xbcxz(jj,kk,i)
         END DO
      END DO
   END IF
   ! Next get the x-direction coefficients
   IF (i == xe) THEN
      ! Reset equation index
      ieq = (nx*ny+nx*nz) + j + (k-1)*ny
      DO jj = ys, j, incy
         DO ii = xs, i, incx
            jeq = ii + (jj-1)*nx
            kpsi(ieq,jeq,n,oct) = zbcyz(ii,jj)
         END DO
      END DO
      DO kk = zs, k, incz
         DO ii = xs, i, incx
            jeq = (nx*ny) + ii + (kk-1)*nx
            kpsi(ieq,jeq,n,oct) = ybcyz(ii,kk)
         END DO
      END DO
      DO kk = zs, k, incz
         DO jj = ys, j, incy
            jeq = (nx*ny+nx*nz) + jj + (kk-1)*ny
            kpsi(ieq,jeq,n,oct) = xbcyz(jj,kk)
         END DO
      END DO
   END IF
END IF

DEALLOCATE(ktmpx, ktmpy, ktmpz)

RETURN
END SUBROUTINE afcm
