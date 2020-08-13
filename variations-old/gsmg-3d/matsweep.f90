SUBROUTINE matsweep(rsx,rsy,rsz,spx,spy,spz,jp,kp,js,ks)

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
INTEGER, INTENT(IN) :: rsx, rsy, rsz, spx, spy, spz
INTEGER :: oct, xs, ys, zs, xe, ye, ze, incx, incy, incz, i, j, k, n, m
REAL*8 :: c, x, y, z, mu, eta, xi, sig, wt
REAL*8, DIMENSION(neq,neq), INTENT(OUT) :: jp
REAL*8, DIMENSION(bcs,neq,8), INTENT(OUT) :: kp
REAL*8, DIMENSION(neq,bcs,8), INTENT(OUT) :: js
REAL*8, DIMENSION(bcs2,bcs2,apo,8), INTENT(OUT) :: ks

! Allocate X-xmat and Y-ymat and Z-zmat and their 'old' counterparts
ALLOCATE(xmat(sdnx,sdny,sdnz), ymat(sdnx,sdny,sdnz,sdnx), zmat(sdnx,sdny,sdnz,sdnx,sdny))
ALLOCATE(xold(sdnx,sdny,sdnz), yold(sdnx,sdny,sdnz), zold(sdnx,sdny,sdnz))

! Allocate Boundary Condition Variables
ALLOCATE(zbcxy(sdnx,sdny,sdnx,sdny),ybcxy(sdnx,sdnz,sdnx,sdny),xbcxy(sdny,sdnz,sdnx,sdny))
ALLOCATE(zbcxz(sdnx,sdny,sdnx),ybcxz(sdnx,sdnz,sdnx),xbcxz(sdny,sdnz,sdnx))
ALLOCATE(zbcyz(sdnx,sdny),zbcxyo(sdnx,sdny),zbcxzo(sdnx,sdny),zbcyzo(sdnx,sdny))
ALLOCATE(ybcyz(sdnx,sdnz),ybcxyo(sdnx,sdnz),ybcxzo(sdnx,sdnz),ybcyzo(sdnx,sdnz))
ALLOCATE(xbcyz(sdny,sdnz),xbcxyo(sdny,sdnz),xbcxzo(sdny,sdnz),xbcyzo(sdny,sdnz))

! Start the sweep
! Start the loops in the quadrants
DO oct = 1, 8
   IF (oct == 1) THEN
      xs = 1
      xe = sdnx
      incx = 1
      ys = 1
      ye = sdny
      incy = 1
      zs = 1
      ze = sdnz
      incz = 1
   ELSE IF (oct == 2) THEN
      xs = sdnx
      xe = 1
      incx = -1
      ys = 1
      ye = sdny
      incy = 1
      zs = 1
      ze = sdnz
      incz = 1
   ELSE IF (oct == 3) THEN
      xs = sdnx
      xe = 1
      incx = -1
      ys = sdny
      ye = 1
      incy = -1
      zs = 1
      ze = sdnz
      incz = 1
   ELSE IF (oct == 4) THEN
      xs = 1
      xe = sdnx
      incx = 1
      ys = sdny
      ye = 1
      incy = -1
      zs = 1
      ze = sdnz
      incz = 1
   ELSE IF (oct == 5) THEN
      xs = 1
      xe = sdnx
      incx = 1
      ys = 1
      ye = sdny
      incy = 1
      zs = sdnz
      ze = 1
      incz = -1
   ELSE IF (oct == 6) THEN
      xs = sdnx
      xe = 1
      incx = -1
      ys = 1
      ye = sdny
      incy = 1
      zs = sdnz
      ze = 1
      incz = -1
   ELSE IF (oct == 7) THEN
      xs = sdnx
      xe = 1
      incx = -1
      ys = sdny
      ye = 1
      incy = -1
      zs = sdnz
      ze = 1
      incz = -1
   ELSE IF (oct == 8) THEN
      xs = 1
      xe = sdnx
      incx = 1
      ys = sdny
      ye = 1
      incy = -1
      zs = sdnz
      ze = 1
      incz = -1
   END IF

   ! Start the loop in all direction
   DO n = 1, apo
      mu  = ang(n,1)
      eta = ang(n,2)
      xi  = ang(n,3)
      wt = w(n) 

      ! Reset the Z matrix at the beginning of a new sweep
      zmat = 0.0

      ! Reset the z-directed BC variables at the beginning of a new sweep
      zbcxy = 0.0
      ybcxy = 0.0
      xbcxy = 0.0

      ! Start the loops    
      DO k = zs, ze, incz
         z = spz*dz(k+rsz)
      
         ! Reset the Y matrix at beginning of new plane
         ymat = 0.0

         ! Reset the y-directed BC varaibles at start of a new plane
         zbcxz = 0.0
         ybcxz = 0.0
         xbcxz = 0.0

         ! Start the plane sweeps
         DO j = ys, ye, incy
            y = spy*dy(j+rsy)

            ! Reset the X matrix at beginning of new row
            xmat = 0.0

            ! Reset the x-directed BC variables at start of new row
            zbcyz = 0.0
            ybcyz = 0.0
            xbcyz = 0.0

            ! Start the sweep in each row
            DO i = xs, xe, incx
               x = spx*dx(i+rsx)

               ! Set up the spatial parameters
               m = mat(i+rsx,j+rsy,k+rsz)
               sig = sigt(m)
               c = sigs(m)/sig     ! Scattering ratio
               ex = mu/(sig*x)
               ey = eta/(sig*y)
               ez = xi/(sig*z)

               ! Get the elements of gamma
               CALL gammas

               ! Call jima to set up matrices for phi and S vectors
               CALL jima(i,j,k,xs,ys,zs,xe,ye,ze,incx,incy,incz,n,c,wt,oct,jp,js)

               ! Call afcm to set up matrices for psi_in vector
               CALL afcm(i,j,k,xs,ys,zs,xe,ye,ze,incx,incy,incz,n,wt,oct,kp,ks)
               
               ! All matrices should be constructed, can exit sweep and solve
            ! End of the mesh sweep
            END DO                ! Rows x
         END DO                   ! Columns y
      END DO                      ! Planes z
   END DO                         ! Angles
END DO                            ! End the loop over all octants

DEALLOCATE(xmat,xold,ymat,yold,zmat,zold)
DEALLOCATE(zbcxy,zbcxyo,zbcxz,zbcxzo,zbcyz,zbcyzo)
DEALLOCATE(ybcxy,ybcxyo,ybcxz,ybcxzo,ybcyz,ybcyzo)
DEALLOCATE(xbcxy,xbcxyo,xbcxz,xbcxzo,xbcyz,xbcyzo)

RETURN
END SUBROUTINE matsweep
