SUBROUTINE afcm(i,j,xs,ys,xe,ye,incx,incy,n,quad,bcs)

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
INTEGER, INTENT(IN) :: i, j, xs, ys, xe, ye, incx, incy, n, quad, bcs
INTEGER :: ii, jj, ieq, jeq
REAL*8, DIMENSION(ordsq,bcs) :: ktmp

! Save the BC matrices for updates
IF (i /= xs) THEN
   xbcyo = xbcy
   ybcyo = ybcy
END IF
IF (j /= ys) THEN
   xbcxo = xbcx(:,:,:,i)
   ybcxo = ybcx(:,:,:,i)
END IF

! Earlier versions that did not give ang. flux out did not solve at xe, ye
! Now solve at xe, ye to use for computing angular flux out

! Update BC matrices with xbcyo
IF (i /= xs) THEN
   DO jj = ys, j, incy
      xbcy(:,:,jj) = MATMUL(gyy,xbcyo(:,:,jj))
      xbcx(:,:,jj,i) = MATMUL(gxy,xbcyo(:,:,jj))
   END DO
END IF

! Update BC matrices with ybcyo
IF (i /= xs) THEN
   DO ii = xs, i-incx, incx
      ybcy(:,:,ii) = MATMUL(gyy,ybcyo(:,:,ii))
      ybcx(:,:,ii,i) = MATMUL(gxy,ybcyo(:,:,ii))
   END DO
END IF

! Update BC matrices with xbcxo
IF (j /= ys) THEN
   DO jj = ys, j-incy, incy
      IF (i /= xs) THEN
         xbcy(:,:,jj) = xbcy(:,:,jj) + MATMUL(gyx,xbcxo(:,:,jj))
         xbcx(:,:,jj,i) = xbcx(:,:,jj,i) + MATMUL(gxx,xbcxo(:,:,jj))
      ELSE IF (i == xs) THEN
         xbcy(:,:,jj) = MATMUL(gyx,xbcxo(:,:,jj))
         xbcx(:,:,jj,i) = MATMUL(gxx,xbcxo(:,:,jj))
      END IF
   END DO
END IF

! Update BC matrices with ybcxo
IF (j /= ys) THEN
   DO ii = xs, i, incx
      IF (ii /= i) THEN
         ybcy(:,:,ii) = ybcy(:,:,ii) + MATMUL(gyx,ybcxo(:,:,ii))
         ybcx(:,:,ii,i) = ybcx(:,:,ii,i) + MATMUL(gxx,ybcxo(:,:,ii))
      ELSE IF (ii == i) THEN
         ybcy(:,:,ii) = MATMUL(gyx,ybcxo(:,:,ii))
         ybcx(:,:,ii,i) = MATMUL(gxx,ybcxo(:,:,ii))
      END IF
   END DO
END IF

! Adjust the values for boundary conditions to their adjacent cells
IF (i == xs) THEN
   xbcy(:,:,j) = gyy
   xbcx(:,:,j,i) = gxy
END IF
IF (j == ys) THEN
   ybcy(:,:,i) = gyx
   ybcx(:,:,i,i) = gxx
END IF

! Update the appropriate kmat given the quadrant of interest
! Call conk to construct the row of ktmp, then set it in kmat
ieq = ((i + (j-1)*nx) - 1)*ordsq
CALL conk(n,i,j,incx,incy,xs,ys,bcs,ktmp)
IF (quad == 1) THEN
   kmat((ieq+1):(ieq+ordsq),:,1) = kmat((ieq+1):(ieq+ordsq),:,1) + ktmp
ELSE IF (quad == 2) THEN
   kmat((ieq+1):(ieq+ordsq),:,2) = kmat((ieq+1):(ieq+ordsq),:,2) + ktmp
ELSE IF (quad == 3) THEN
   kmat((ieq+1):(ieq+ordsq),:,3) = kmat((ieq+1):(ieq+ordsq),:,3) + ktmp
ELSE IF (quad == 4) THEN
   kmat((ieq+1):(ieq+ordsq),:,4) = kmat((ieq+1):(ieq+ordsq),:,4) + ktmp
END IF

! Use the data in the xbc ybc variables at xe ye to compute kpsi
! Start with y-direction coefficients
IF (j == ye) THEN
   ! Reset equation index
   ieq = (i-1)*order
   DO ii = xs, i, incx
      jeq = (ii-1)*order
      kpsi((ieq+1):(ieq+order),(jeq+1):(jeq+order),n,quad) = ybcx(:,:,ii,i)
   END DO
   DO jj = ys, j, incy
      jeq = nx*order + (jj-1)*order
      kpsi((ieq+1):(ieq+order),(jeq+1):(jeq+order),n,quad) = xbcx(:,:,jj,i)
   END DO
END IF
! Next get the x-direction coefficients
IF (i == xe) THEN
   ! Reset equation index
   ieq = nx*order + (j-1)*order
   DO ii = xs, i, incx
      jeq = (ii-1)*order
      kpsi((ieq+1):(ieq+order),(jeq+1):(jeq+order),n,quad) = ybcy(:,:,ii)
   END DO
   DO jj = ys, j, incy
      jeq = nx*order + (jj-1)*order
      kpsi((ieq+1):(ieq+order),(jeq+1):(jeq+order),n,quad) = xbcy(:,:,jj)
   END DO
END IF
 
RETURN
END SUBROUTINE afcm
