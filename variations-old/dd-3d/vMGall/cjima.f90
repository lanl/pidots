SUBROUTINE cjima(i,j,k,xs,ys,zs,xe,ye,ze,incx,incy,incz,n,c,oct)

!-------------------------------------------------------------
!
!  Jacobian Iteration MAtrix
!  Directs the equations to set up the matrices that work
!   with the phi and source vectors. takes in parameters from
!   matsweep and solves for the updates to cjmat and cjpsi 
!
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER, INTENT(IN) :: i, j, k, xs, ys, zs, xe, ye, ze, incx, incy, incz, n, oct
REAL*8, INTENT(IN) :: c
REAL*8 :: wt
INTEGER :: ii, jj, kk, ieq, jeq

! Set the quad weight
wt = w(n)

! Save the X, Y, and Z matrices for updates
IF (i /= xs) THEN
   xold = xmat
END IF
IF (j /= ys) THEN
   yold(:,:,:) = ymat(:,:,:,i)
END IF
IF (k /= zs) THEN
   zold(:,:,:) = zmat(:,:,:,i,j)
END IF

! Older version without angular flux out did not compute anything at xe, ye, ze
! Now compute that info and use it to construct cjpsi

! Update X, Y, Z matrices with old X matrix
IF (i /= xs) THEN
   DO kk = zs, k, incz
      DO jj = ys, j, incy
         DO ii = xs, i, incx
            xmat(ii,jj,kk) = gyzyz*xold(ii,jj,kk)
            ymat(ii,jj,kk,i) = gxzyz*xold(ii,jj,kk)
            zmat(ii,jj,kk,i,j) = gxyyz*xold(ii,jj,kk)
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
               xmat(ii,jj,kk) = xmat(ii,jj,kk) + gyzxz*yold(ii,jj,kk)
               ymat(ii,jj,kk,i) = ymat(ii,jj,kk,i) + gxzxz*yold(ii,jj,kk)
               zmat(ii,jj,kk,i,j) = zmat(ii,jj,kk,i,j) + gxyxz*yold(ii,jj,kk)
            ELSE IF (i == xs) THEN
               xmat(ii,jj,kk) = gyzxz*yold(ii,jj,kk)
               ymat(ii,jj,kk,i) = gxzxz*yold(ii,jj,kk)
               zmat(ii,jj,kk,i,j) = gxyxz*yold(ii,jj,kk)
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
               xmat(ii,jj,kk) = xmat(ii,jj,kk) + gyzxy*zold(ii,jj,kk)
               ymat(ii,jj,kk,i) = ymat(ii,jj,kk,i) + gxzxy*zold(ii,jj,kk)
               zmat(ii,jj,kk,i,j) = zmat(ii,jj,kk,i,j) + gxyxy*zold(ii,jj,kk)
            ELSE IF (i == xs .AND. j == ys) THEN
               xmat(ii,jj,kk) = gyzxy*zold(ii,jj,kk)
               ymat(ii,jj,kk,i) = gxzxy*zold(ii,jj,kk)
               zmat(ii,jj,kk,i,j) = gxyxy*zold(ii,jj,kk)
            END IF
         END DO
      END DO
   END DO
END IF

! Append X, Y, and Z matrices
zmat(i,j,k,i,j) = c*gxya
ymat(i,j,k,i) = c*gxza
xmat(i,j,k) = c*gyza

! Start to formulate Jacobian matrix, diagonal blocks
! Depends on the 'tpose' flag how the indices are incremented
! tpose = 1 --> make transposed CJphi and CJpsi
IF (tpose == 1) THEN
   jeq = i + (j-1)*npx + (k-1)*cxys
   cjmat(jeq,jeq) = cjmat(jeq,jeq) + wt*c*gaa

   ! Off-diagonal contribution from X
   IF (i /= xs) THEN
      DO kk = zs, k, incz
         DO jj = ys, j, incy
            DO ii = xs, i, incx
               ieq = ii + (jj-1)*npx + (kk-1)*cxys
               cjmat(ieq,jeq) = cjmat(ieq,jeq) + wt*gayz*xold(ii,jj,kk)
            END DO
         END DO
      END DO
   END IF

   ! Off-diagonal contribution from Y
   IF (j /= ys) THEN
      DO kk = zs, k, incz
         DO jj = ys, j, incy
            DO ii = xs, i, incx
               ieq = ii + (jj-1)*npx + (kk-1)*cxys
               cjmat(ieq,jeq) = cjmat(ieq,jeq) + wt*gaxz*yold(ii,jj,kk)
            END DO
         END DO
      END DO
   END IF

   ! Off-diagonal contribution from Z
   IF (k /= zs) THEN
      DO kk = zs, k, incz
         DO jj = ys, j, incy
            DO ii = xs, i, incx
               ieq = ii + (jj-1)*npx + (kk-1)*cxys
               cjmat(ieq,jeq) = cjmat(ieq,jeq) + wt*gaxy*zold(ii,jj,kk)
            END DO
         END DO
      END DO
   END IF

   ! Now check to see if the computational cell is at end of X or Y dimension,
   ! Place an xmat or ymat set of values if it is
   ! First get the z-direction coefficients
   IF (k == ze) THEN
      ! Reset equation index
      jeq = (n-1)*cxys + i + (j-1)*npx
      DO kk = zs, k, incz
         DO jj = ys, j, incy
            DO ii = xs, i, incx
               ieq = ii + (jj-1)*npx + (kk-1)*cxys
               cjpsi(ieq,jeq,oct) = zmat(ii,jj,kk,i,j)
            END DO
         END DO
      END DO
   END IF
   IF (j == ye) THEN
      ! Reset equation index
      jeq = apo*cxys + (n-1)*cxzs + i + (k-1)*npx
      DO kk = zs, k, incz
         DO jj = ys, j, incy
            DO ii = xs, i, incx
               ieq = ii + (jj-1)*npx + (kk-1)*cxys
               cjpsi(ieq,jeq,oct) = ymat(ii,jj,kk,i)
            END DO
         END DO
      END DO
   END IF
   IF (i == xe) THEN
      ! Reset equation index
      jeq = apo*(cxys+cxzs) + (n-1)*cyzs + j + (k-1)*npy
      DO kk = zs, k, incz
         DO jj = ys, j, incy
            DO ii = xs, i, incx
               ieq = ii + (jj-1)*npx + (kk-1)*cxys
               cjpsi(ieq,jeq,oct) = xmat(ii,jj,kk)
            END DO
         END DO
      END DO
   END IF
!-----------------------------------------------------------------------
! If tpose = 0 --> make normal CJphi and CJpsi
ELSE
   ieq = i + (j-1)*npx + (k-1)*cxys
   cjmat(ieq,ieq) = cjmat(ieq,ieq) + wt*c*gaa

   ! Off-diagonal contribution from X
   IF (i /= xs) THEN
      DO kk = zs, k, incz
         DO jj = ys, j, incy
            DO ii = xs, i, incx
               jeq = ii + (jj-1)*npx + (kk-1)*cxys
               cjmat(ieq,jeq) = cjmat(ieq,jeq) + wt*gayz*xold(ii,jj,kk)
            END DO
         END DO
      END DO
   END IF

   ! Off-diagonal contribution from Y
   IF (j /= ys) THEN
      DO kk = zs, k, incz
         DO jj = ys, j, incy
            DO ii = xs, i, incx
               jeq = ii + (jj-1)*npx + (kk-1)*cxys
               cjmat(ieq,jeq) = cjmat(ieq,jeq) + wt*gaxz*yold(ii,jj,kk)
            END DO
         END DO
      END DO
   END IF

   ! Off-diagonal contribution from Z
   IF (k /= zs) THEN
      DO kk = zs, k, incz
         DO jj = ys, j, incy
            DO ii = xs, i, incx
               jeq = ii + (jj-1)*npx + (kk-1)*cxys
               cjmat(ieq,jeq) = cjmat(ieq,jeq) + wt*gaxy*zold(ii,jj,kk)
            END DO
         END DO
      END DO
   END IF

   ! Now check to see if the computational cell is at end of X or Y dimension,
   ! Place an xmat or ymat set of values if it is
   ! First get the z-direction coefficients
   IF (k == ze) THEN
      ! Reset equation index
      ieq = (n-1)*cxys + i + (j-1)*npx
      DO kk = zs, k, incz
         DO jj = ys, j, incy
            DO ii = xs, i, incx
               jeq = ii + (jj-1)*npx + (kk-1)*cxys
               cjpsi(ieq,jeq,oct) = zmat(ii,jj,kk,i,j)
            END DO
         END DO
      END DO
   END IF
   IF (j == ye) THEN
      ! Reset equation index
      ieq = apo*cxys + (n-1)*cxzs + i + (k-1)*npx
      DO kk = zs, k, incz
         DO jj = ys, j, incy
            DO ii = xs, i, incx
               jeq = ii + (jj-1)*npx + (kk-1)*cxys
               cjpsi(ieq,jeq,oct) = ymat(ii,jj,kk,i)
            END DO
         END DO
      END DO
   END IF
   IF (i == xe) THEN
      ! Reset equation index
      ieq = apo*(cxys+cxzs) + (n-1)*cyzs + j + (k-1)*npy
      DO kk = zs, k, incz
         DO jj = ys, j, incy
            DO ii = xs, i, incx
               jeq = ii + (jj-1)*npx + (kk-1)*cxys
               cjpsi(ieq,jeq,oct) = xmat(ii,jj,kk)
            END DO
         END DO
      END DO
   END IF
END IF

RETURN
END SUBROUTINE cjima
