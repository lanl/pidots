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
INTEGER :: oct, xs, ys, zs, xe, ye, ze, incx, incy, incz, i, j, k, n, m
REAL*8 :: c, x, y, z, mu, eta, xi, sig

! Allocate X-xmat and Y-ymat and Z-zmat and their 'old' counterparts
ALLOCATE(ymat(nx), zmat(nx,ny))

! Start the sweep
! Start the loops in the quadrants
DO oct = 1, 8
   IF (oct == 1) THEN
      xs = 1
      xe = nx
      incx = 1
      ys = 1
      ye = ny
      incy = 1
      zs = 1
      ze = nz
      incz = 1
   ELSE IF (oct == 2) THEN
      xs = nx
      xe = 1
      incx = -1
      ys = 1
      ye = ny
      incy = 1
      zs = 1
      ze = nz
      incz = 1
   ELSE IF (oct == 3) THEN
      xs = nx
      xe = 1
      incx = -1
      ys = ny
      ye = 1
      incy = -1
      zs = 1
      ze = nz
      incz = 1
   ELSE IF (oct == 4) THEN
      xs = 1
      xe = nx
      incx = 1
      ys = ny
      ye = 1
      incy = -1
      zs = 1
      ze = nz
      incz = 1
   ELSE IF (oct == 5) THEN
      xs = 1
      xe = nx
      incx = 1
      ys = 1
      ye = ny
      incy = 1
      zs = nz
      ze = 1
      incz = -1
   ELSE IF (oct == 6) THEN
      xs = nx
      xe = 1
      incx = -1
      ys = 1
      ye = ny
      incy = 1
      zs = nz
      ze = 1
      incz = -1
   ELSE IF (oct == 7) THEN
      xs = nx
      xe = 1
      incx = -1
      ys = ny
      ye = 1
      incy = -1
      zs = nz
      ze = 1
      incz = -1
   ELSE IF (oct == 8) THEN
      xs = 1
      xe = nx
      incx = 1
      ys = ny
      ye = 1
      incy = -1
      zs = nz
      ze = 1
      incz = -1
   END IF

   ! Start the loop in all direction
   DO n = 1, apo
      mu  = ang(n,1)
      eta = ang(n,2)
      xi  = ang(n,3)   

      ! Reset the Z matrix at the beginning of a new sweep
      zmat = 0.0

      ! Start the loops    
      DO k = zs, ze, incz
         z = dz(k)
      
         ! Reset the Y matrix at beginning of new plane
         ymat = 0.0

         ! Start the plane sweeps
         DO j = ys, ye, incy
            y = dy(j)

            ! Reset the X matrix at beginning of new row
            xmat = 0.0

            ! Start the sweep in each row
            DO i = xs, xe, incx
               x = dx(i)

               ! Set up the spatial parameters
               m = mat(i,j,k)
               sig = sigt(m,g)
               c = sigs(m,g,g)/sig     ! Scattering ratio
               ex = mu/(sig*x)
               ey = eta/(sig*y)
               ez = xi/(sig*z)

               ! Get the elements of gamma
               CALL gammas

               ! Call jima to set up matrices for phi and S vectors
               CALL jima(i,j,k,xs,ys,zs,incx,incy,incz,n,c)

            ! End of the mesh sweep
            END DO                ! Rows x
         END DO                   ! Columns y
      END DO                      ! Planes z
   END DO                         ! Angles
END DO                            ! End the loop over all octants

DEALLOCATE(ymat,zmat)

RETURN
END SUBROUTINE matsweep
