SUBROUTINE conkn(n,i,j,k,incx,incy,incz,xs,ys,zs,ktmpz,ktmpy,ktmpx)

!-----------------------------------------------------------
!
! Constructs the ktmp row for a given octant and cell and 
!  returns the values to afcm.
!
! Takes in from afcm: All the computed values not in solvar
!     n,i,j,k,incx,incy,incz,xs,ys,zs,bcs
! Returns to afcm a temporary matrix ktmp that fills in the
!  particular ktmp#
!
!------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER, INTENT(IN) :: n, i, j, k 
INTEGER, INTENT(IN) :: incx, incy, incz, xs, ys, zs
INTEGER :: ii, jj, kk, jeq
REAL*8 :: wt
REAL*8, DIMENSION(ordcb,xys), INTENT(OUT) :: ktmpz
REAL*8, DIMENSION(ordcb,xzs), INTENT(OUT) :: ktmpy
REAL*8, DIMENSION(ordcb,yzs), INTENT(OUT) :: ktmpx

! Initialize ktmps
ktmpz = 0.0
ktmpy = 0.0
ktmpx = 0.0

! Set the weight
wt = w(n)

! All the bc matrices and gamma matrices needed come from solvar

! Perform the operations to update from the different BCs
! Start with zbc's
IF (k /= zs) THEN
   DO jj = ys, j, incy
      DO ii = xs, i-incx, incx
         jeq = ((ii-1)+(jj-1)*nx)*ordsq
         ktmpz(:,(jeq+1):(jeq+ordsq)) = ktmpz(:,(jeq+1):(jeq+ordsq)) + wt*MATMUL(gayz,zbcyzo(:,:,ii,jj))
      END DO
   END DO
   DO jj = ys, j-incy, incy
      DO ii = xs, i, incx
         jeq = ((ii-1)+(jj-1)*nx)*ordsq
         ktmpz(:,(jeq+1):(jeq+ordsq)) = ktmpz(:,(jeq+1):(jeq+ordsq)) + wt*MATMUL(gaxz,zbcxzo(:,:,ii,jj))
      END DO
   END DO
   DO jj = ys, j, incy
      DO ii = xs, i, incx
         jeq = ((ii-1)+(jj-1)*nx)*ordsq
         ktmpz(:,(jeq+1):(jeq+ordsq)) = ktmpz(:,(jeq+1):(jeq+ordsq)) + wt*MATMUL(gaxy,zbcxyo(:,:,ii,jj))
      END DO
   END DO
ELSE IF (k == zs) THEN
   DO jj = ys, j, incy
      DO ii = xs, i, incx
         jeq = ((ii-1)+(jj-1)*nx)*ordsq
         IF (ii /= i) THEN
            ktmpz(:,(jeq+1):(jeq+ordsq)) = ktmpz(:,(jeq+1):(jeq+ordsq)) + wt*MATMUL(gayz,zbcyzo(:,:,ii,jj))
         END IF
         IF (jj /= j) THEN
            ktmpz(:,(jeq+1):(jeq+ordsq)) = ktmpz(:,(jeq+1):(jeq+ordsq)) + wt*MATMUL(gaxz,zbcxzo(:,:,ii,jj))
         END IF
         IF (ii == i .AND. jj == j) THEN
            ktmpz(:,(jeq+1):(jeq+ordsq)) = ktmpz(:,(jeq+1):(jeq+ordsq)) + wt*gaxy
         END IF
      END DO
   END DO
END IF

! Update values from ybc's
IF (j /= ys) THEN
   DO kk = zs, k, incz
      DO ii = xs, i-incx, incx
         jeq = ((ii-1)+(kk-1)*nx)*ordsq
         ktmpy(:,(jeq+1):(jeq+ordsq)) = ktmpy(:,(jeq+1):(jeq+ordsq)) + wt*MATMUL(gayz,ybcyzo(:,:,ii,kk))
      END DO
   END DO
   DO kk = zs, k, incz
      DO ii = xs, i, incx
         jeq = ((ii-1)+(kk-1)*nx)*ordsq
         ktmpy(:,(jeq+1):(jeq+ordsq)) = ktmpy(:,(jeq+1):(jeq+ordsq)) + wt*MATMUL(gaxz,ybcxzo(:,:,ii,kk))
      END DO
   END DO
   DO kk = zs, k-incz, incz
      DO ii = xs, i, incx
         jeq = ((ii-1)+(kk-1)*nx)*ordsq
         ktmpy(:,(jeq+1):(jeq+ordsq)) = ktmpy(:,(jeq+1):(jeq+ordsq)) + wt*MATMUL(gaxy,ybcxyo(:,:,ii,kk))
      END DO
   END DO
ELSE IF (j == ys) THEN
   DO kk = zs, k, incz
      DO ii = xs, i, incx
         jeq = ((ii-1)+(kk-1)*nx)*ordsq
         IF (ii /= i) THEN
            ktmpy(:,(jeq+1):(jeq+ordsq)) = ktmpy(:,(jeq+1):(jeq+ordsq)) + wt*MATMUL(gayz,ybcyzo(:,:,ii,kk))
         END IF
         IF (kk /= k) THEN
            ktmpy(:,(jeq+1):(jeq+ordsq)) = ktmpy(:,(jeq+1):(jeq+ordsq)) + wt*MATMUL(gaxy,ybcxyo(:,:,ii,kk))
         END IF
         IF (ii == i .AND. kk == k) THEN
            ktmpy(:,(jeq+1):(jeq+ordsq)) = ktmpy(:,(jeq+1):(jeq+ordsq)) + wt*gaxz
         END IF
      END DO
   END DO
END IF

! Update values from xbc's
IF (i /= xs) THEN
   DO kk = zs, k, incz
      DO jj = ys, j, incy
         jeq = ((jj-1)+(kk-1)*ny)*ordsq
         ktmpx(:,(jeq+1):(jeq+ordsq)) = ktmpx(:,(jeq+1):(jeq+ordsq)) + wt*MATMUL(gayz,xbcyzo(:,:,jj,kk))
      END DO
   END DO
   DO kk = zs, k, incz
      DO jj = ys, j-incy, incy
         jeq = ((jj-1)+(kk-1)*ny)*ordsq
         ktmpx(:,(jeq+1):(jeq+ordsq)) = ktmpx(:,(jeq+1):(jeq+ordsq)) + wt*MATMUL(gaxz,xbcxzo(:,:,jj,kk))
      END DO
   END DO
   DO kk = zs, k-incz, incz
      DO jj = ys, j, incy
         jeq = ((jj-1)+(kk-1)*ny)*ordsq
         ktmpx(:,(jeq+1):(jeq+ordsq)) = ktmpx(:,(jeq+1):(jeq+ordsq)) + wt*MATMUL(gaxy,xbcxyo(:,:,jj,kk))
      END DO
   END DO
ELSE IF (i == xs) THEN
   DO kk = zs, k, incz
      DO jj = ys, j, incy
         jeq = ((jj-1)+(kk-1)*ny)*ordsq
         IF (jj /= j) THEN
            ktmpx(:,(jeq+1):(jeq+ordsq)) = ktmpx(:,(jeq+1):(jeq+ordsq)) + wt*MATMUL(gaxz,xbcxzo(:,:,jj,kk))
         END IF
         IF (kk /= k) THEN
            ktmpx(:,(jeq+1):(jeq+ordsq)) = ktmpx(:,(jeq+1):(jeq+ordsq)) + wt*MATMUL(gaxy,xbcxyo(:,:,jj,kk))
         END IF
         IF (jj == j .AND. kk == k) THEN
            ktmpx(:,(jeq+1):(jeq+ordsq)) = ktmpx(:,(jeq+1):(jeq+ordsq)) + wt*gayz
         END IF
      END DO
   END DO
END IF

! Finished updating ktmp for given cell
RETURN

END SUBROUTINE conkn
