SUBROUTINE matsweep(g,bcs)

!----------------------------------------------------------------
!
!  Directs the sweep across all computational cells for all
!  angles in all quadrants for constructing the matrices.
!
!  For each cell, calls the weight, conab, cong, gammas routines
!
!  Then calls jima and afcm to set construct the matrices
! 
!----------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER, INTENT(IN) :: g, bcs
INTEGER :: quad, xs, ys, xe, ye, incx, incy, n, i, j, sgm, sge, m, ord, neq
REAL*8 :: c, x, y, sig, mu, eta

! Size amat, bmat, gmat
neq = ordsq + 2*order
ALLOCATE(amat(neq,neq), bmat(neq,neq), gmat(neq,neq))

! Allocate all the gamma submatrices sizes
ALLOCATE(gaa(ordsq,ordsq), gax(ordsq,order), gay(ordsq,order))
ALLOCATE(gxa(order,ordsq), gxx(order,order), gxy(order,order))
ALLOCATE(gya(order,ordsq), gyx(order,order), gyy(order,order))

! Allocate the matrices for jima
ALLOCATE(xmat(order,ordsq,nx,ny), ymat(nx,order,ordsq,nx,ny))
ALLOCATE(xold(order,ordsq,nx,ny), yold(order,ordsq,nx,ny))

! Allocate the matrices for afcm
ALLOCATE(xbcy(order,order,ny), xbcyo(order,order,ny), xbcxo(order,order,ny))
ALLOCATE(xbcx(order,order,ny,nx))
ALLOCATE(ybcy(order,order,nx), ybcyo(order,order,nx), ybcxo(order,order,nx))
ALLOCATE(ybcx(order,order,nx,nx))

! Start the sweep
! Start the loops in the quadrants
DO quad = 1, 4
   IF (quad == 1) THEN
      xs = 1
      xe = nx
      incx = 1
      ys = 1
      ye = ny
      incy = 1
   ELSE IF (quad == 2) THEN
      xs = nx
      xe = 1
      incx = -1
      ys = 1
      ye = ny
      incy = 1
   ELSE IF (quad == 3) THEN
      xs = nx
      xe = 1
      incx = -1
      ys = ny
      ye = 1
      incy = -1
   ELSE IF (quad == 4) THEN
      xs = 1
      xe = nx
      incx = 1
      ys = ny
      ye = 1
      incy = -1
   END IF

   ! Start the loop in all direction
   DO n = 1, apo
      mu  = ang(n,1)
      eta = ang(n,2)   

      ! Reset the Y matrix at the beginning of a new sweep
      ymat = 0.0
      
      ! Reset the y-directed BC variables at start of new sweep
      xbcx = 0.0
      ybcx = 0.0

   ! Start the loops    
   DO j = ys, ye, incy
      sge = incy
      y = dy(j)
      
      ! Reset the X matrix at beginning of new row
      xmat = 0.0

      ! Reset the x-directed BC variables at start of a new row
      xbcy = 0.0
      ybcy = 0.0      
         
      ! Start the row sweeps
      DO i = xs, xe, incx
         sgm = incx
         ! Get the spatial weights
         x = dx(i)
         m = mat(i,j)
         sig = sigt(m,g)
         ord = lambda
         c = sigs(m,g,g)/sig     ! Scattering ratio
         ! Call for the calculation of the spatial weights
         CALL weight(ord,x,y,sig,mu,eta)
         
         ! Construct amat and bmat, A and B
         CALL conab(sgm,sge)
         ! Construct the gamma matrix
         CALL cong
         ! Unload the gamma matrix into useful sub-matrices
         CALL gammas

         ! Call jima to set up matrices for phi and S vectors
         CALL jima(i,j,xs,ys,xe,ye,incx,incy,n,c,quad)

         ! Call afcm to set up matrices for psi_in vector
         CALL afcm(i,j,xs,ys,xe,ye,incx,incy,n,quad,bcs)

         ! All matrices should be constructed, can exit sweep and solve

      ! End of the mesh sweep
      END DO   ! Rows
   END DO   ! Columns
   END DO   ! Angles
! End the loop over all quadrants
END DO

DEALLOCATE(amat,bmat,gmat)
DEALLOCATE(gaa,gax,gay,gxa,gxx,gxy,gya,gyx,gyy)
DEALLOCATE(xmat,xold,ymat,yold)
DEALLOCATE(xbcy,xbcyo,xbcx,xbcxo,ybcy,ybcyo,ybcx,ybcxo)

RETURN
END SUBROUTINE matsweep
