SUBROUTINE conk(n,i,j,incx,incy,xs,ys,bcs,ktmp)

!------------------------------------------------------
!
! Constructs the ktmp row for a given quadrant and cell
!  and returns the values to jima.
!
! Takes in from jima: All the computed values not in 
!  solvar
!   n,i,j,incx,incy,xs,ys,bcs
! Returns to jima a temporary matrix ktmp that fills in
!  the particular ktmp#
!
!------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER, INTENT(IN) :: n, i, j, incx, incy, xs, ys, bcs
INTEGER :: ii, jj, jndx, jeq
REAL*8 :: wt
REAL*8, DIMENSION(ordsq,bcs), INTENT(OUT) :: ktmp

! Initialize ktmp
ktmp = 0.0

! Set the weight
wt = w(n)

! Set the initial jndx
jndx = (n-1)*nx*order

! All the bc matrices and gamma matrices needed come from solvar

! Perform the operations to update from the different BCs
! Start with ybc's
IF (j /= ys) THEN
   DO ii = xs, i, incx
      jeq = jndx + (ii-1)*order
      ktmp(:,(jeq+1):(jeq+order)) = ktmp(:,(jeq+1):(jeq+order)) + wt*MATMUL(gax,ybcxo(:,:,ii))
   END DO
   DO ii = xs, i-incx, incx
      jeq = jndx + (ii-1)*order
      ktmp(:,(jeq+1):(jeq+order)) = ktmp(:,(jeq+1):(jeq+order)) + wt*MATMUL(gay,ybcyo(:,:,ii))
   END DO
ELSE
   DO ii = xs, i, incx
      jeq = jndx + (ii-1)*order
      IF (ii /= i) THEN
         ktmp(:,(jeq+1):(jeq+order)) = ktmp(:,(jeq+1):(jeq+order)) + wt*MATMUL(gay,ybcyo(:,:,ii))
      ELSE
         ktmp(:,(jeq+1):(jeq+order)) = ktmp(:,(jeq+1):(jeq+order)) + wt*gax
      END IF
   END DO
END IF

jndx = apo*nx*order + (n-1)*ny*order
! Update values from xbc's
IF (i /= xs) THEN
   DO jj = ys, j, incy
      jeq = jndx + (jj-1)*order
      ktmp(:,(jeq+1):(jeq+order)) = ktmp(:,(jeq+1):(jeq+order)) + wt*MATMUL(gay,xbcyo(:,:,jj))
   END DO
   DO jj = ys, j-incy, incy
      jeq = jndx + (jj-1)*order
      ktmp(:,(jeq+1):(jeq+order)) = ktmp(:,(jeq+1):(jeq+order)) + wt*MATMUL(gax,xbcxo(:,:,jj))
   END DO
ELSE
   DO jj = ys, j, incy
      jeq = jndx + (jj-1)*order
      IF (jj /= j) THEN
         ktmp(:,(jeq+1):(jeq+order)) = ktmp(:,(jeq+1):(jeq+order)) + wt*MATMUL(gax,xbcxo(:,:,jj))
      ELSE
         ktmp(:,(jeq+1):(jeq+order)) = ktmp(:,(jeq+1):(jeq+order)) + wt*gay
      END IF
   END DO
END IF

! Finished updating ktmp for given cell
RETURN

END SUBROUTINE conk
