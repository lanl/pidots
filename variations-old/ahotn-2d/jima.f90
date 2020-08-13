SUBROUTINE jima(i,j,xs,ys,xe,ye,incx,incy,n,c,quad)

!-------------------------------------------------------------
!
!  Jacobian Iteration MAtrix
!  Directs the equations to set up the matrices that work
!   with the phi and source vectors. Takes in parameters from
!   matsweep and solves for the updates to jmat and jpsi
! 
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER, INTENT(IN) :: i, j, xs, ys, xe, ye, incx, incy, n, quad
REAL*8, INTENT(IN) :: c
INTEGER :: ii, jj, ieq, jeq

! Save the X and Y matrices for updates
IF (i /= xs) THEN
   xold = xmat
END IF
IF (j /= ys) THEN
   yold(:,:,:,:) = ymat(i,:,:,:,:)
END IF

! Older version without angular flux out did not compute anything at xe, ye
! Now compute that info and use it to construct jpsi

! Update X and Y matrices with the old X matrix
IF (i /= xs) THEN
   DO jj = ys, j, incy
      DO ii = xs, i, incx
          xmat(:,:,ii,jj) = MATMUL(gyy,xold(:,:,ii,jj))
          ymat(i,:,:,ii,jj) = MATMUL(gxy,xold(:,:,ii,jj))
      END DO
   END DO
END IF
         
! Update X and Y matrices with the old Y matrix
IF (j /= ys) THEN
   DO jj = ys, j, incy
      DO ii = xs, i, incx
         IF (i /= xs) THEN
            xmat(:,:,ii,jj) = xmat(:,:,ii,jj) + MATMUL(gyx,yold(:,:,ii,jj))
            ymat(i,:,:,ii,jj) = ymat(i,:,:,ii,jj) + MATMUL(gxx,yold(:,:,ii,jj))
         ELSE IF (i == xs) THEN
            xmat(:,:,ii,jj) = MATMUL(gyx,yold(:,:,ii,jj))
            ymat(i,:,:,ii,jj) = MATMUL(gxx,yold(:,:,ii,jj))
         END IF
      END DO
   END DO
END IF

! Append X and Y matrices
ymat(i,:,:,i,j) = c*gxa
xmat(:,:,i,j) = c*gya

! Start to formulate Jacobian matrix, diagonal blocks
ieq = ((i + (j-1)*nx) - 1)*ordsq
jmat((ieq+1):(ieq+ordsq),(ieq+1):(ieq+ordsq)) = jmat((ieq+1):(ieq+ordsq),(ieq+1):(ieq+ordsq)) + w(n)*c*gaa

! Off-diagonal contribution from X
IF (i /= xs) THEN
   DO jj = ys, j, incy
      DO ii = xs, i, incx
         jeq = ((ii + (jj-1)*nx) - 1)*ordsq
         jmat((ieq+1):(ieq+ordsq),(jeq+1):(jeq+ordsq)) = jmat((ieq+1):(ieq+ordsq),(jeq+1):(jeq+ordsq)) &
                                                          + w(n)*MATMUL(gay,xold(:,:,ii,jj))
      END DO
   END DO
END IF

! Off-diagonal contribution from Y
IF (j /= ys) THEN
   DO jj = ys, j, incy
      DO ii = xs, i, incx
         jeq = ((ii + (jj-1)*nx) - 1)*ordsq
         jmat((ieq+1):(ieq+ordsq),(jeq+1):(jeq+ordsq)) = jmat((ieq+1):(ieq+ordsq),(jeq+1):(jeq+ordsq)) &
                                                          + w(n)*MATMUL(gax,yold(:,:,ii,jj))
      END DO
   END DO
END IF

! Now check to see if the computational cell is at end of X or Y dimension,
! Place an xmat or ymat set of values if it is
! First get the y-direction coefficients
IF (j == ye) THEN
   ! Reset equation index
   ieq = (n-1)*nx*order + (i-1)*order
   DO jj = ys, j, incy
      DO ii = xs, i, incx
         jeq = ((ii-1) + (jj-1)*nx)*ordsq
         jpsi((ieq+1):(ieq+order),(jeq+1):(jeq+ordsq),quad) = ymat(i,:,:,ii,jj)
      END DO
   END DO
END IF
! Next get the x-direction coefficients
IF (i == xe) THEN
   ! Reset equation index
   ieq = apo*nx*order + (n-1)*ny*order + (j-1)*order
   DO jj = ys, j, incy
      DO ii = xs, i, incx
         jeq = ((ii-1) + (jj-1)*nx)*ordsq
         jpsi((ieq+1):(ieq+order),(jeq+1):(jeq+ordsq),quad) = xmat(:,:,ii,jj)
      END DO
   END DO
END IF

RETURN
END SUBROUTINE jima
