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
INTEGER :: oct, xs, ys, zs, xe, ye, ze, incx, incy, incz, i, j, k, n, sgm, sge, sgx, m, ord, eqs
REAL*8 :: c, x, y, z, mu, eta, xi, sig

! Size amat, bmat, gmat
eqs = ordcb + 3*ordsq
ALLOCATE(amat(eqs,eqs), bmat(eqs,eqs), gmat(eqs,eqs))

! Allocate all the gamma submatrices sizes
ALLOCATE(gaa(ordcb,ordcb),  gaxy(ordcb,ordsq),  gaxz(ordcb,ordsq),  gayz(ordcb,ordsq))
ALLOCATE(gxya(ordsq,ordcb), gxyxy(ordsq,ordsq), gxyxz(ordsq,ordsq), gxyyz(ordsq,ordsq))
ALLOCATE(gxza(ordsq,ordcb), gxzxy(ordsq,ordsq), gxzxz(ordsq,ordsq), gxzyz(ordsq,ordsq))
ALLOCATE(gyza(ordsq,ordcb), gyzxy(ordsq,ordsq), gyzxz(ordsq,ordsq), gyzyz(ordsq,ordsq))

! Allocate X-xmat and Y-ymat and Z-zmat and their 'old' counterparts
ALLOCATE(xmat(ordsq,ordcb,nx,ny,nz), ymat(ordsq,ordcb,nx,ny,nz,nx), zmat(ordsq,ordcb,nx,ny,nz,nx,ny))
ALLOCATE(xold(ordsq,ordcb,nx,ny,nz), yold(ordsq,ordcb,nx,ny,nz), zold(ordsq,ordcb,nx,ny,nz))

! Allocate Boundary Condition Variables
ALLOCATE(zbcxy(ordsq,ordsq,nx,ny,nx,ny),ybcxy(ordsq,ordsq,nx,nz,nx,ny),xbcxy(ordsq,ordsq,ny,nz,nx,ny))
ALLOCATE(zbcxz(ordsq,ordsq,nx,ny,nx),ybcxz(ordsq,ordsq,nx,nz,nx),xbcxz(ordsq,ordsq,ny,nz,nx))
ALLOCATE(zbcyz(ordsq,ordsq,nx,ny),zbcxyo(ordsq,ordsq,nx,ny),zbcxzo(ordsq,ordsq,nx,ny),zbcyzo(ordsq,ordsq,nx,ny))
ALLOCATE(ybcyz(ordsq,ordsq,nx,nz),ybcxyo(ordsq,ordsq,nx,nz),ybcxzo(ordsq,ordsq,nx,nz),ybcyzo(ordsq,ordsq,nx,nz))
ALLOCATE(xbcyz(ordsq,ordsq,ny,nz),xbcxyo(ordsq,ordsq,ny,nz),xbcxzo(ordsq,ordsq,ny,nz),xbcyzo(ordsq,ordsq,ny,nz))

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

      ! Reset the z-directed BC variables at the beginning of a new sweep
      zbcxy = 0.0
      ybcxy = 0.0
      xbcxy = 0.0

      ! Start the loops    
      DO k = zs, ze, incz
         sgx = incz
         z = dz(k)
      
         ! Reset the Y matrix at beginning of new plane
         ymat = 0.0

         ! Reset the y-directed BC varaibles at start of a new plane
         zbcxz = 0.0
         ybcxz = 0.0
         xbcxz = 0.0

         ! Start the plane sweeps
         DO j = ys, ye, incy
            sge = incy
            y = dy(j)

            ! Reset the X matrix at beginning of new row
            xmat = 0.0

            ! Reset the x-directed BC variables at start of new row
            zbcyz = 0.0
            ybcyz = 0.0
            xbcyz = 0.0

            ! Start the sweep in each row
            DO i = xs, xe, incx
               sgm = incx

               ! Get the spatial weights
               x = dx(i)
               m = mat(i,j,k)
               sig = sigt(m,g)
               ord = lambda
               c = sigs(m,g,g)/sig     ! Scattering ratio

               ! Call for the calculation of the spatial weights
               CALL weight(ord,x,y,z,sig,mu,eta,xi)

               ! Construct amat and bmat, A and B
               CALL conab(sgm,sge,sgx)
               ! Construct the gamma matrix
               CALL cong
               ! Unload the gamma matrix into useful sub-matrices
               CALL gammas

               ! Call jima to set up matrices for phi and S vectors
               CALL jima(i,j,k,xs,ys,zs,xe,ye,ze,incx,incy,incz,n,c,oct)

               ! Call afcm to set up matrices for psi_in vector
               CALL afcm(i,j,k,xs,ys,zs,xe,ye,ze,incx,incy,incz,n,oct)
               
               ! All matrices should be constructed, can exit sweep and solve
            ! End of the mesh sweep
            END DO                ! Rows x
         END DO                   ! Columns y
      END DO                      ! Planes z
   END DO                         ! Angles
END DO                            ! End the loop over all octants

DEALLOCATE(amat,bmat,gmat)
DEALLOCATE(gaa,gaxy,gaxz,gayz,gxya,gxza,gyza)
DEALLOCATE(gxyxy,gxyxz,gxyyz,gxzxy,gxzxz,gxzyz,gyzxy,gyzxz,gyzyz)
DEALLOCATE(xmat,xold,ymat,yold,zmat,zold)
DEALLOCATE(zbcxy,zbcxyo,zbcxz,zbcxzo,zbcyz,zbcyzo)
DEALLOCATE(ybcxy,ybcxyo,ybcxz,ybcxzo,ybcyz,ybcyzo)
DEALLOCATE(xbcxy,xbcxyo,xbcxz,xbcxzo,xbcyz,xbcyzo)

RETURN
END SUBROUTINE matsweep
