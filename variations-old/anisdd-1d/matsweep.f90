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
INTEGER :: oct, xs, xe, incx, i, n, m
REAL*8 :: x, mu, sig
REAL*8, DIMENSION(0:anord) :: tsxs
REAL*8, DIMENSION(2,2) :: amat
REAL*8, DIMENSION(2,(sord+1)) :: bmat, gmat

! Allocate the two multi-element pieces of gamma
ALLOCATE(gaa(0:anord), gxa(0:anord))

! Allocate X-xmat and Y-ymat and Z-zmat and their 'old' counterparts
ALLOCATE(xmat(0:anord,nx), xold(0:anord,nx))

! Start the sweep
DO oct = 1, 2
   IF (oct == 1) THEN
      xs = 1
      xe = nx
      incx = 1
   ELSE IF (oct == 2) THEN
      xs = nx
      xe = 1
      incx = -1
   END IF

   ! Start the loop in all direction
   DO n = 1, apo
      mu  = ang(n)

      ! Reset the X matrix at beginning of new sweep
      xmat = 0.0

      ! Reset the x-directed BC variable at start of new sweep
      xbc = 0.0

      ! Start the sweep
      DO i = xs, xe, incx
         x = dx(i)

         ! Set up the spatial parameters
         m = mat(i)
         sig = sigt(m,g)
         tsxs = sigs(m,:,g,g)
         ex = mu/x

         ! Get the elements of gamma
         CALL gammas(incx,mu,sig,amat,bmat,gmat)

         ! Call jima to set up matrices for phi and S vectors
         CALL jima(i,xs,xe,incx,n,oct,mu,tsxs)
         
         ! Call afcm to set up matrices for psi_in vector
         CALL afcm(i,xs,xe,incx,n,oct,mu)
               
         ! All matrices should be constructed, can exit sweep and solve
      ! End of the mesh sweep
      END DO                ! Cells
   END DO                   ! Angles
END DO                      ! End the loop over all octants

DEALLOCATE(gaa,gxa)
DEALLOCATE(xmat,xold)

RETURN
END SUBROUTINE matsweep
