SUBROUTINE matsweep(g)

!-----------------------------------------------------------------
!
!  Directs the sweep across all computational cells for all
!   angles in all quadrants for constructing the matrices
!
!  For each cell, calls the weight, conab, cong, gammas routines
!
!  The calls jima and afcm to set construct the matrices
!
!-----------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER, INTENT(IN) :: g
INTEGER :: oct, xs, ys, xe, ye, incx, incy, i, j, n, m
REAL*8 :: x, y, mu, eta, omega, sig, pi, tmp
REAL*8, DIMENSION(0:anord) :: tsxs

tmp = 0.0
pi = 2.0*ACOS(tmp)

! Allocate X-xmat and Y-ymat and Z-zmat and their 'old' counterparts
ALLOCATE(xmat(nmom,nx,ny), ymat(nmom,nx,ny,nx))
ALLOCATE(xold(nmom,nx,ny), yold(nmom,nx,ny))

! Allocate Boundary Condition Variables
ALLOCATE(ybcx(nx,nx),xbcx(ny,nx))
ALLOCATE(ybcy(nx),ybcxo(nx),ybcyo(nx))
ALLOCATE(xbcy(ny),xbcxo(ny),xbcyo(ny))

! Allocate gamma sub-matrices
ALLOCATE(gaa(nmom), gxa(nmom), gya(nmom))

! Start the sweep
! Start the loops in the quadrants
DO oct = 1, 4
   IF (oct == 1) THEN
      xs = 1
      xe = nx
      incx = 1
      ys = 1
      ye = ny
      incy = 1
   ELSE IF (oct == 2) THEN
      xs = nx
      xe = 1
      incx = -1
      ys = 1
      ye = ny
      incy = 1
   ELSE IF (oct == 3) THEN
      xs = nx
      xe = 1
      incx = -1
      ys = ny
      ye = 1
      incy = -1
   ELSE IF (oct == 4) THEN
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
      omega = omg(n)
      IF (incy == -1) omega = pi - omega

      ! Reset the Y matrix at beginning of new angle
      ymat = 0.0

      ! Reset the y-directed BC variables at start of a new angle
      ybcx = 0.0
      xbcx = 0.0

      ! Start the plane sweeps
      DO j = ys, ye, incy
         y = dy(j)

         ! Reset the X matrix at beginning of new row
         xmat = 0.0

         ! Reset the x-directed BC variables at start of new row
         ybcy = 0.0
         xbcy = 0.0

         ! Start the sweep in each row
         DO i = xs, xe, incx
            x = dx(i)

            ! Set up the spatial parameters
            m = mat(i,j)
            sig = sigt(m,g)
            tsxs = sigs(:,m,g,g)
            ex = mu/x
            ey = eta/y

            ! Get the elements of gamma
            CALL gammas(incx,mu,sig,omega)

            ! Call jima to set up matrices for phi and S vectors
            CALL jima(i,j,xs,ys,xe,ye,incx,incy,n,oct,mu,omega,tsxs)

            ! Call afcm to set up matrices for psi_in vector
            CALL afcm(i,j,xs,ys,xe,ye,incx,incy,n,oct,mu,omega)

            ! All matrices should be constructed, can exit sweep and solve
         ! End of the mesh sweep
         END DO                ! Rows x
      END DO                   ! Columns y
   END DO                      ! Angles
END DO                         ! End the loop over all octants

DEALLOCATE(xmat,xold,ymat,yold)
DEALLOCATE(ybcx,ybcy,ybcxo,ybcyo)
DEALLOCATE(xbcx,xbcy,xbcxo,xbcyo)
DEALLOCATE(gaa,gxa,gya)

RETURN
END SUBROUTINE matsweep
