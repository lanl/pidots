SUBROUTINE sweep(g,e)

!-------------------------------------------------------------
!
!  Sweeps across the 3-D matrix
!   Starts at top, far, right corner (mu, eta, xi < 0), then sweeps
!   down all planes and rows, accounting for reflection if necessary. Then
!   sweeps up all planes and rows for xi>0
! 
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER, INTENT(IN) :: g
INTEGER :: xs, xe, ys, ye, zs, ze, incx, incy, incz, nfy, nfz
INTEGER :: i, j, k, m, n, ydir, xdir, zdir, indx, jndx, oct
REAL*8 :: sig, mu, eta, xi, x, y, z, c, a, b, fx
REAL*8, DIMENSION(nx, 2) :: fy
REAL*8, DIMENSION(nx, ny, 2, 2) :: fz
REAL*8, DIMENSION(nx,ny,nz), INTENT(IN) :: e

! Initialize the flux solution to zero
f(:,:,:,g) = 0.0

! Start with loop over all angles
DO n = 1, apo
   ! Set up the angles
   mu  = ang(n,1)
   eta = ang(n,2)
   xi  = ang(n,3)
  
   ! Initialize the incoming z-flux for z-hi bc (O5-8)
   indx = (n-1)*(nx*ny)
   DO j = 1, ny
      DO i = 1, nx
         jndx = indx + i + (j-1)*nx
         fz(i,j,1,1) = psii(jndx,7)
         fz(i,j,2,1) = psii(jndx,8)
         fz(i,j,1,2) = psii(jndx,6)
         fz(i,j,2,2) = psii(jndx,5)
      END DO
   END DO
               
   ! Loop over xi<0 then xi>0
   DO zdir = 1, 2
      IF (zdir == 1) THEN
         zs = nz
         ze = 1
         incz = -1
      ELSE IF (zdir == 2) THEN
         zs = 1
         ze = nz
         incz = 1
         ! Reset the incoming z-flux for the z-lo bc (octants 1-4)
         indx = (n-1)*(nx*ny)
         DO j = 1, ny
            DO i = 1, nx
               jndx = indx + i + (j-1)*nx
               fz(i,j,1,1) = psii(jndx,3)
               fz(i,j,2,1) = psii(jndx,4)
               fz(i,j,1,2) = psii(jndx,2)
               fz(i,j,2,2) = psii(jndx,1)
            END DO
         END DO
      END IF

   ! Start the loop in the negative z-direction, then do positive z-direction
   DO k = zs, ze, incz
      z = dz(k)
 
      ! Reset the index for ybc's
      indx = apo*(nx*ny) + (n-1)*(nx*nz)
      ! Set the initial incoming y-flux depending on the z-dir
      IF (zdir == 1) THEN
         DO i = 1, nx
            jndx = indx + i + (k-1)*nx
            fy(i,1) = psii(jndx,7)
            fy(i,2) = psii(jndx,8)
         END DO
      ELSE IF (zdir == 2) THEN
         DO i = 1, nx
            jndx = indx + i + (k-1)*nx
            fy(i,1) = psii(jndx,3)
            fy(i,2) = psii(jndx,4)
         END DO
      END IF

      ! Loop over eta<0 then eta>0
      DO ydir = 1, 2
         IF (ydir == 1) THEN
            ys = ny
            ye = 1
            incy = -1
            nfz = 1
         ELSE IF (ydir == 2) THEN
            ys = 1
            ye = ny
            incy = 1
            nfz = 2
            ! Reset incoming y-flux for the y-lo bc
            indx = apo*(nx*ny) + (n-1)*(nx*nz)
            IF (zdir == 1) THEN
               DO i = 1, nx
                  jndx = indx + i + (k-1)*nx
                  fy(i,1) = psii(jndx,6)
                  fy(i,2) = psii(jndx,5)
               END DO
            ELSE IF (zdir == 2) THEN
               DO i = 1, nx
                  jndx = indx + i + (k-1)*nx
                  fy(i,1) = psii(jndx,2)
                  fy(i,2) = psii(jndx,1)
               END DO
            END IF
         END IF

      ! Start the loop in negative y-direction, then do positve y-direction
      DO j = ys, ye, incy
         y = dy(j)
         
         ! Reset the index for the xbc's
         indx = apo*(nx*ny+nx*nz) + (n-1)*(ny*nz)
         ! Set the incoming x-flux with the x-hi bc (octants 2,3,6,7)
         ! Depends on both zdir and ydir
         IF (zdir == 1) THEN
            IF (ydir == 1) THEN
               jndx = indx + j + (k-1)*ny
               fx = psii(jndx,7)
            ELSE IF (ydir == 2) THEN
               jndx = indx + j + (k-1)*ny
               fx = psii(jndx,6)
            END IF
         ELSE IF (zdir == 2) THEN
            IF (ydir == 1) THEN
               jndx = indx + j + (k-1)*ny
               fx = psii(jndx,3)
            ELSE IF (ydir == 2) THEN
               jndx = indx + j + (k-1)*ny
               fx = psii(jndx,2)
            END IF
         END IF

         ! Perform two loops, one in negative x-direction, then positive
         DO xdir = 1, 2
            IF (xdir == 1) THEN
               xs = nx
               xe = 1
               incx = -1
               nfy = 1
               IF (ydir == 1 .AND. zdir == 1) oct = 7
               IF (ydir == 1 .AND. zdir == 2) oct = 3
               IF (ydir == 2 .AND. zdir == 1) oct = 6
               IF (ydir == 2 .AND. zdir == 2) oct = 2
            ELSE IF (xdir == 2) THEN
               xs = 1
               xe = nx
               incx = 1
               nfy = 2
               IF (ydir == 1 .AND. zdir == 1) oct = 8
               IF (ydir == 1 .AND. zdir == 2) oct = 4
               IF (ydir == 2 .AND. zdir == 1) oct = 5
               IF (ydir == 2 .AND. zdir == 2) oct = 1
               ! Reset the incoming x-flux for the x-lo bc
               indx = apo*(nx*ny+nx*nz) + (n-1)*(ny*nz)
               IF (zdir == 1) THEN
                  IF (ydir == 1) THEN
                     jndx = indx + j + (k-1)*ny
                     fx = psii(jndx,8)
                  ELSE IF (ydir == 2) THEN
                     jndx = indx + j + (k-1)*ny
                     fx = psii(jndx,5)
                  END IF
               ELSE IF (zdir == 2) THEN
                  IF (ydir == 1) THEN
                     jndx = indx + j + (k-1)*ny
                     fx = psii(jndx,4)
                  ELSE IF (ydir == 2) THEN
                     jndx = indx + j + (k-1)*ny
                     fx = psii(jndx,1)
                  END IF
               END IF
            END IF

         ! Start the loop in the negative x-direction
         DO i = xs, xe, incx
            x = dx(i)
            m = mat(i,j,k)
            sig = sigt(m,g)
            c = sigs(m,g,g)/sig     ! Scattering ratio

            ! Get spatial parameters
            ex = mu/(sig*x)
            ey = eta/(sig*y)
            ez = xi/(sig*z)
         
            ! Initialize 'a' coefficient
            a = 0.0

            ! Contributions from DD outgoing substitutions
            a = a + 2.0*(ex + ey + ez)
            ! Contribution from total interaction
            a = a + 1.0

            ! Initialize rhs coefficient
            b = 0.0

            ! Contribution from incoming fluxes
            b = b + 2.0*(ex*fx + ey*fy(i,nfy) + ez*fz(i,j,nfy,nfz))
            ! Contribution from scattering plus fixed sources
            b = b + c*e(i,j,k) + s(i,j,k,g)/sig            

            ! Solve for the flux
            b = b/a

            ! Update the scalar flux
            f(i,j,k,g) = f(i,j,k,g) + w(n)*b
            
            ! Compute the outgoing fluxes with the DD equations
            fx = 2.0*b - fx
            fy(i,nfy) = 2.0*b - fy(i,nfy)
            fz(i,j,nfy,nfz) = 2.0*b - fz(i,j,nfy,nfz)

            ! Store the new angular fluxes at ze ye xe
            IF (k == ze) THEN
               indx = (n-1)*(nx*ny) + i + (j-1)*nx
               psio(indx,oct,g) = fz(i,j,nfy,nfz)
            END IF
            IF (j == ye) THEN
               indx = apo*(nx*ny) + (n-1)*(nx*nz) + i + (k-1)*nx
               psio(indx,oct,g) = fy(i,nfy)
            END IF
            IF (i == xe) THEN
               indx = apo*(nx*ny+nx*nz) + (n-1)*(ny*nz) + j + (k-1)*ny
               psio(indx,oct,g) = fx
            END IF

         ! End loop over x cells
         END DO
         ! End loop over negative and positive x-directions
         END DO
   
      ! End loop over y cells
      END DO
      ! End loop over negative and positve y-directions
      END DO
   
   ! End loop over z cells
   END DO
   ! End loop over negative and positive z-directions
   END DO

! End loop over angles
END DO         

RETURN
END SUBROUTINE sweep
