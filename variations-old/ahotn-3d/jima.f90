SUBROUTINE jima(i,j,k,xs,ys,zs,xe,ye,ze,incx,incy,incz,n,c,oct)

!-------------------------------------------------------------
!
!  Jacobian Iteration MAtrix
!  Directs the equations to set up the matrices that work
!   with the phi and source vectors. takes in parameters from
!   matsweep and solves for the updates to jmat and jpsi 
!
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER, INTENT(IN) :: i, j, k, xs, ys, zs, xe, ye, ze, incx, incy, incz, n, oct
INTEGER :: nxy, ii, jj, kk, ieq, jeq
REAL*8, INTENT(IN) :: c
REAL*8 :: wt

! Save multiplications
nxy = nx*ny

! Set the quad weight
wt = w(n)

! Save the X, Y, and Z matrices for updates
IF (i /= xs) xold = xmat
IF (j /= ys) yold(:,:,:,:,:) = ymat(:,:,:,:,:,i)
IF (k /= zs) zold(:,:,:,:,:) = zmat(:,:,:,:,:,i,j)

! Older version without angular flux out did not compute anything at xe, ye, ze
! Now compute that info and use it to construct jpsi

! Update X, Y, Z matrices with old X matrix
IF (i /= xs) THEN
   DO kk = zs, k, incz
      DO jj = ys, j, incy
         DO ii = xs, i, incx
            xmat(:,:,ii,jj,kk) = MATMUL(gyzyz,xold(:,:,ii,jj,kk))
            ymat(:,:,ii,jj,kk,i) = MATMUL(gxzyz,xold(:,:,ii,jj,kk))
            zmat(:,:,ii,jj,kk,i,j) = MATMUL(gxyyz,xold(:,:,ii,jj,kk))
         END DO
      END DO
   END DO
END IF

! Update X, Y, and Z matrices with old Y matrix
IF (j /= ys) THEN
   DO kk = zs, k, incz
      DO jj = ys, j, incy
         DO ii = xs, i, incx
            IF (i /= xs) THEN
               xmat(:,:,ii,jj,kk) = xmat(:,:,ii,jj,kk) + MATMUL(gyzxz,yold(:,:,ii,jj,kk))
               ymat(:,:,ii,jj,kk,i) = ymat(:,:,ii,jj,kk,i) + MATMUL(gxzxz,yold(:,:,ii,jj,kk))
               zmat(:,:,ii,jj,kk,i,j) = zmat(:,:,ii,jj,kk,i,j) + MATMUL(gxyxz,yold(:,:,ii,jj,kk))
            ELSE IF (i == xs) THEN
               xmat(:,:,ii,jj,kk) = MATMUL(gyzxz,yold(:,:,ii,jj,kk))
               ymat(:,:,ii,jj,kk,i) = MATMUL(gxzxz,yold(:,:,ii,jj,kk))
               zmat(:,:,ii,jj,kk,i,j) = MATMUL(gxyxz,yold(:,:,ii,jj,kk))
            END IF
         END DO
      END DO
   END DO
END IF

! Update X, Y, and Z matrices with old Z matrix
IF (k /= zs) THEN
   DO kk = zs, k, incz
      DO jj = ys, j, incy
         DO ii = xs, i, incx
            IF (i /= xs .OR. j /= ys) THEN
               xmat(:,:,ii,jj,kk) = xmat(:,:,ii,jj,kk) + MATMUL(gyzxy,zold(:,:,ii,jj,kk))
               ymat(:,:,ii,jj,kk,i) = ymat(:,:,ii,jj,kk,i) + MATMUL(gxzxy,zold(:,:,ii,jj,kk))
               zmat(:,:,ii,jj,kk,i,j) = zmat(:,:,ii,jj,kk,i,j) + MATMUL(gxyxy,zold(:,:,ii,jj,kk))
            ELSE IF (i == xs .AND. j == ys) THEN
               xmat(:,:,ii,jj,kk) = MATMUL(gyzxy,zold(:,:,ii,jj,kk))
               ymat(:,:,ii,jj,kk,i) = MATMUL(gxzxy,zold(:,:,ii,jj,kk))
               zmat(:,:,ii,jj,kk,i,j) = MATMUL(gxyxy,zold(:,:,ii,jj,kk))
            END IF
         END DO
      END DO
   END DO
END IF

! Append X, Y, and Z matrices
zmat(:,:,i,j,k,i,j) = c*gxya
ymat(:,:,i,j,k,i) = c*gxza
xmat(:,:,i,j,k) = c*gyza

!-----------------------------------------------------------------------------
! Start to formulate Jacobian matrix, diagonal blocks
! Depends on the 'tpose' flag how the indices are incremented
! tpose = 1 --> make transposed Jphi and Jpsi
IF (tpose == 1) THEN
   jeq = ((i-1) + (j-1)*nx + (k-1)*nxy)*ordcb
   jmat((jeq+1):(jeq+ordcb),(jeq+1):(jeq+ordcb)) = jmat((jeq+1):(jeq+ordcb),(jeq+1):(jeq+ordcb)) &
                                                   + wt*c*TRANSPOSE(gaa)

   ! Off-diagonal contribution from X
   IF (i /= xs) THEN
      DO kk = zs, k, incz
         DO jj = ys, j, incy
            DO ii = xs, i, incx
               ieq = ((ii-1) + (jj-1)*nx + (kk-1)*nxy)*ordcb
               jmat((ieq+1):(ieq+ordcb),(jeq+1):(jeq+ordcb)) = jmat((ieq+1):(ieq+ordcb),(jeq+1):(jeq+ordcb)) &
                                                             + w(n)*TRANSPOSE(MATMUL(gayz,xold(:,:,ii,jj,kk)))
            END DO
         END DO
      END DO
   END IF

   ! Off-diagonal contribution from Y
   IF (j /= ys) THEN
      DO kk = zs, k, incz
         DO jj = ys, j, incy
            DO ii = xs, i, incx
               ieq = ((ii-1) + (jj-1)*nx + (kk-1)*nxy)*ordcb
               jmat((ieq+1):(ieq+ordcb),(jeq+1):(jeq+ordcb)) = jmat((ieq+1):(ieq+ordcb),(jeq+1):(jeq+ordcb)) &
                                                                + w(n)*TRANSPOSE(MATMUL(gaxz,yold(:,:,ii,jj,kk)))
            END DO
         END DO
      END DO
   END IF

   ! Off-diagonal contribution from Z
   IF (k /= zs) THEN
      DO kk = zs, k, incz
         DO jj = ys, j, incy
            DO ii = xs, i, incx
               ieq = ((ii-1) + (jj-1)*nx + (kk-1)*nxy)*ordcb
               jmat((ieq+1):(ieq+ordcb),(jeq+1):(jeq+ordcb)) = jmat((ieq+1):(ieq+ordcb),(jeq+1):(jeq+ordcb)) &
                                                                + w(n)*TRANSPOSE(MATMUL(gaxy,zold(:,:,ii,jj,kk)))
            END DO
         END DO
      END DO
   END IF

   ! Now check to see if the computational cell is at end of X or Y dimension,
   ! Place an xmat or ymat set of values if it is
   ! First get the z-direction coefficients
   IF (k == ze) THEN
      ! Reset equation index
      jeq = (n-1)*xys + ((i-1) + (j-1)*nx)*ordsq
      DO kk = zs, k, incz
         DO jj = ys, j, incy
            DO ii = xs, i, incx
               ieq = ((ii-1) + (jj-1)*nx + (kk-1)*nxy)*ordcb
               jpsi((ieq+1):(ieq+ordcb),(jeq+1):(jeq+ordsq),oct) = TRANSPOSE(zmat(:,:,ii,jj,kk,i,j))
            END DO
         END DO
      END DO
   END IF
   IF (j == ye) THEN
      ! Reset equation index
      jeq = apo*xys + (n-1)*xzs + ((i-1) + (k-1)*nx)*ordsq
      DO kk = zs, k, incz
         DO jj = ys, j, incy
            DO ii = xs, i, incx
               ieq = ((ii-1) + (jj-1)*nx + (kk-1)*nxy)*ordcb
               jpsi((ieq+1):(ieq+ordcb),(jeq+1):(jeq+ordsq),oct) = TRANSPOSE(ymat(:,:,ii,jj,kk,i))
            END DO
         END DO
      END DO
   END IF
   IF (i == xe) THEN
      ! Reset equation index
      jeq = apo*(xys+xzs) + (n-1)*yzs + ((j-1) + (k-1)*ny)*ordsq
      DO kk = zs, k, incz
         DO jj = ys, j, incy
            DO ii = xs, i, incx
               ieq = ((ii-1) + (jj-1)*nx + (kk-1)*nxy)*ordcb
               jpsi((ieq+1):(ieq+ordcb),(jeq+1):(jeq+ordsq),oct) = TRANSPOSE(xmat(:,:,ii,jj,kk))
            END DO
         END DO
      END DO
   END IF
!-----------------------------------------------------------------------
! If tpose = 0 --> make normal Jphi and Jpsi
ELSE
   ieq = ((i-1) + (j-1)*nx + (k-1)*nxy)*ordcb
   jmat((ieq+1):(ieq+ordcb),(ieq+1):(ieq+ordcb)) = jmat((ieq+1):(ieq+ordcb),(ieq+1):(ieq+ordcb)) + wt*c*gaa

   ! Off-diagonal contribution from X
   IF (i /= xs) THEN
      DO kk = zs, k, incz
         DO jj = ys, j, incy
            DO ii = xs, i, incx
               jeq = ((ii-1) + (jj-1)*nx + (kk-1)*nxy)*ordcb
               jmat((ieq+1):(ieq+ordcb),(jeq+1):(jeq+ordcb)) = jmat((ieq+1):(ieq+ordcb),(jeq+1):(jeq+ordcb)) &
                                                             + w(n)*MATMUL(gayz,xold(:,:,ii,jj,kk))
            END DO
         END DO
      END DO
   END IF

   ! Off-diagonal contribution from Y
   IF (j /= ys) THEN
      DO kk = zs, k, incz
         DO jj = ys, j, incy
            DO ii = xs, i, incx
               jeq = ((ii-1) + (jj-1)*nx + (kk-1)*nxy)*ordcb
               jmat((ieq+1):(ieq+ordcb),(jeq+1):(jeq+ordcb)) = jmat((ieq+1):(ieq+ordcb),(jeq+1):(jeq+ordcb)) &
                                                                + w(n)*MATMUL(gaxz,yold(:,:,ii,jj,kk))
            END DO
         END DO
      END DO
   END IF

   ! Off-diagonal contribution from Z
   IF (k /= zs) THEN
      DO kk = zs, k, incz
         DO jj = ys, j, incy
            DO ii = xs, i, incx
               jeq = ((ii-1) + (jj-1)*nx + (kk-1)*nxy)*ordcb
               jmat((ieq+1):(ieq+ordcb),(jeq+1):(jeq+ordcb)) = jmat((ieq+1):(ieq+ordcb),(jeq+1):(jeq+ordcb)) &
                                                                + w(n)*MATMUL(gaxy,zold(:,:,ii,jj,kk))
            END DO
         END DO
      END DO
   END IF

   ! Now check to see if the computational cell is at end of X or Y dimension,
   ! Place an xmat or ymat set of values if it is
   ! First get the z-direction coefficients
   IF (k == ze) THEN
      ! Reset equation index
      ieq = (n-1)*xys + ((i-1) + (j-1)*nx)*ordsq
      DO kk = zs, k, incz
         DO jj = ys, j, incy
            DO ii = xs, i, incx
               jeq = ((ii-1) + (jj-1)*nx + (kk-1)*nxy)*ordcb
               jpsi((ieq+1):(ieq+ordsq),(jeq+1):(jeq+ordcb),oct) = zmat(:,:,ii,jj,kk,i,j)
            END DO
         END DO
      END DO
   END IF
   IF (j == ye) THEN
      ! Reset equation index
      ieq = apo*xys + (n-1)*xzs + ((i-1) + (k-1)*nx)*ordsq
      DO kk = zs, k, incz
         DO jj = ys, j, incy
            DO ii = xs, i, incx
               jeq = ((ii-1) + (jj-1)*nx + (kk-1)*nxy)*ordcb
               jpsi((ieq+1):(ieq+ordsq),(jeq+1):(jeq+ordcb),oct) = ymat(:,:,ii,jj,kk,i)
            END DO
         END DO
      END DO
   END IF
   IF (i == xe) THEN
      ! Reset equation index
      ieq = apo*(xys+xzs) + (n-1)*yzs + ((j-1) + (k-1)*ny)*ordsq
      DO kk = zs, k, incz
         DO jj = ys, j, incy
            DO ii = xs, i, incx
               jeq = ((ii-1) + (jj-1)*nx + (kk-1)*nxy)*ordcb
               jpsi((ieq+1):(ieq+ordsq),(jeq+1):(jeq+ordcb),oct) = xmat(:,:,ii,jj,kk)
            END DO
         END DO
      END DO
   END IF
END IF

RETURN
END SUBROUTINE jima
