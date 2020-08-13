SUBROUTINE jima(i,j,k,xs,ys,zs,incx,incy,incz,n,c)

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
INTEGER, INTENT(IN) :: i, j, k, xs, ys, zs, incx, incy, incz, n
REAL*8, INTENT(IN) :: c
REAL*8 :: wt

! Set the quad weight
wt = w(n)

! Save the X, Y, and Z matrices for updates
IF (i /= xs) THEN
   xold = xmat
END IF
IF (j /= ys) THEN
   yold = ymat(i)
END IF
IF (k /= zs) THEN
   zold = zmat(i,j)
END IF

zmat(i,j) = c*gxya
ymat(i) = c*gxza
xmat = c*gyza

jmat(4,i,j,k) = jmat(4,i,j,k) + wt*c*gaa

IF (incx > 0) THEN
   jmat(1,i,j,k) = jmat(1,i,j,k) + wt*gayz*xold
ELSE
   jmat(5,i,j,k) = jmat(5,i,j,k) + wt*gayz*xold
END IF

IF (incy > 0) THEN
   jmat(2,i,j,k) = jmat(2,i,j,k) + wt*gaxz*yold
ELSE
   jmat(6,i,j,k) = jmat(6,i,j,k) + wt*gaxz*yold
END IF

IF (incz > 0) THEN
   jmat(3,i,j,k) = jmat(3,i,j,k) + wt*gaxy*zold
ELSE
   jmat(7,i,j,k) = jmat(7,i,j,k) + wt*gaxy*zold
END IF

RETURN
END SUBROUTINE jima
