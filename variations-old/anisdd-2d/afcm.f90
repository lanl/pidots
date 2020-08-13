SUBROUTINE afcm(i,j,xs,ys,xe,ye,incx,incy,n,oct,mu,omega)

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
INTEGER, INTENT(IN) :: i, j, xs, ys, xe, ye, incx, incy, n, oct
INTEGER :: ii, jj, ieq, jeq, indx
REAL*8, INTENT(IN) :: mu, omega
REAL*8, DIMENSION(:,:), ALLOCATABLE :: ktmpy
REAL*8, DIMENSION(:,:), ALLOCATABLE :: ktmpx

! Save the BC matrices for updates
IF (i /= xs) THEN
   ybcyo = ybcy
   xbcyo = xbcy
END IF
IF (j /= ys) THEN
   ybcxo = ybcx(:,i)
   xbcxo = xbcx(:,i)
END IF

! Update BC matrices with xbcyo (BC data in x-direction into cell)
IF (i /= xs) THEN
   DO jj = ys, j, incy
      xbcy(jj) = gyy*xbcyo(jj)
      xbcx(jj,i) = gxy*xbcyo(jj)
   END DO
END IF

! Update BC matrices with ybcyo
IF (i /= xs) THEN
   DO ii = xs, i-incx, incx
      ybcy(ii) = gyy*ybcyo(ii)
      ybcx(ii,i) = gxy*ybcyo(ii)
   END DO
END IF

! Update BC matrices with xbcxo (BC data in y-dir into cell)
IF (j /= ys) THEN
   DO jj = ys, j-incy, incy
      IF (i /= xs) THEN
         xbcy(jj) = xbcy(jj) + gyx*xbcxo(jj)
         xbcx(jj,i) = xbcx(jj,i) + gxx*xbcxo(jj)
      ELSE
         xbcy(jj) = gyx*xbcxo(jj)
         xbcx(jj,i) = gxx*xbcxo(jj)
      END IF
   END DO
END IF

! Update BC matrices with ybcxo
IF (j /= ys) THEN
   DO ii = xs, i, incx
      IF (ii /= i) THEN
         ybcy(ii) = ybcy(ii) + gyx*ybcxo(ii)
         ybcx(ii,i) = ybcx(ii,i) + gxx*ybcxo(ii)
      ELSE
         ybcy(ii) = gyx*ybcxo(ii)
         ybcx(ii,i) = gxx*ybcxo(ii)
      END IF
   END DO
END IF

! Adjust the values for BCs to their adjacent cells
IF (i == xs) THEN
   xbcy(j) = gyy
   xbcx(j,i) = gxy
END IF
IF (j == ys) THEN
   ybcy(i) = gyx
   ybcx(i,i) = gxx
END IF

!------------------------------------------------------------------------------
! Update the appropriate kmat given the octant of interest
! Call conk to construct the row of kmat, then set it in kmat#
! Set according to 'tpose' flag
IF (tpose == 1) THEN
   ALLOCATE(ktmpy(nx,nmom), ktmpx(ny,nmom))
   jeq = ((i-1) + (j-1)*nx)*nmom
   CALL conkt(n,i,j,incx,incy,xs,ys,mu,omega,ktmpy,ktmpx)
   indx = (n-1)*nx
   kmat(indx+1:indx+nx,jeq+1:jeq+nmom,oct) = kmat(indx+1:indx+nx,jeq+1:jeq+nmom,oct) + ktmpy
   indx = apo*nx + (n-1)*ny
   kmat(indx+1:indx+ny,jeq+1:jeq+nmom,oct) = kmat(indx+1:indx+ny,jeq+1:jeq+nmom,oct) + ktmpx

   ! Use the data in the xbc, ybc varialbes at xe ye to compute kpsi
   ! Get the y-direction coefficients
   IF (j == ye) THEN
      ! Reset equation index
      jeq = i
      DO ii = xs, i, incx
         ieq = ii
         kpsi(ieq,jeq,n,oct) = ybcx(ii,i)
      END DO
      DO jj = ys, j, incy
         ieq = nx + jj
         kpsi(ieq,jeq,n,oct) = xbcx(jj,i)
      END DO
   END IF
   ! Next get the x-direction coefficients
   IF (i == xe) THEN
      ! Reset equation index
      jeq = nx + j
      DO ii = xs, i, incx
         ieq = ii
         kpsi(ieq,jeq,n,oct) = ybcy(ii)
      END DO
      DO jj = ys, j, incy
         ieq = nx + jj
         kpsi(ieq,jeq,n,oct) = xbcy(jj)
      END DO
   END IF
!-------------------------------------------------------------------------------
ELSE
   ALLOCATE(ktmpy(nmom,nx), ktmpx(nmom,ny))
   ieq = ((i-1) + (j-1)*nx)*nmom
   CALL conk(n,i,j,incx,incy,xs,ys,mu,omega,ktmpy,ktmpx)
   indx = (n-1)*nx
   kmat(ieq+1:ieq+nmom,indx+1:indx+nx,oct) = kmat(ieq+1:ieq+nmom,indx+1:indx+nx,oct) + ktmpy
   indx = apo*nx + (n-1)*ny
   kmat(ieq+1:ieq+nmom,indx+1:indx+ny,oct) = kmat(ieq+1:ieq+nmom,indx+1:indx+ny,oct) + ktmpx

   ! Use the data in the xbc, ybc variables at xe ye to compute kpsi
   ! Get the y-direction coefficients
   IF (j == ye) THEN
      ! Reset equation index
      ieq = i
      DO ii = xs, i, incx
         jeq = ii
         kpsi(ieq,jeq,n,oct) = ybcx(ii,i)
      END DO
      DO jj = ys, j, incy
         jeq = nx + jj
         kpsi(ieq,jeq,n,oct) = xbcx(jj,i)
      END DO
   END IF
   ! Next get the x-direction coefficients
   IF (i == xe) THEN
      ! Reset equation index
      ieq = nx + j
      DO ii = xs, i, incx
         jeq = ii
         kpsi(ieq,jeq,n,oct) = ybcy(ii)
      END DO
      DO jj = ys, j, incy
         jeq = nx + jj
         kpsi(ieq,jeq,n,oct) = xbcy(jj)
      END DO
   END IF
END IF

DEALLOCATE(ktmpy, ktmpx)

RETURN
END SUBROUTINE afcm
