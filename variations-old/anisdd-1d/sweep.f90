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
INTEGER :: xs, xe, incx, i, m, n, xdir, l, ieq
REAL*8 :: sig, mu, x, a, b, fx, plm1, pl, ea
REAL*8, DIMENSION(0:anord) :: tsxs
REAL*8, DIMENSION(nx*sord), INTENT(IN) :: e
REAL*8, DIMENSION(nx*sord) :: et

! Initialize the flux solution to zero
f(:,g) = 0.0

! Start with loop over all angles
DO n = 1, apo
   ! Set up the angles
   mu  = ang(n)
  
   ! Set the incoming x-flux 
   fx = psii(n,1)

   ! Perform two loops, one in positive x-direction, then negative
   DO xdir = 1, 2
      IF (xdir == 1) THEN
         xs = 1
         xe = nx
         incx = 1
      ELSE IF (xdir == 2) THEN
         xs = nx
         xe = 1
         incx = -1
         ! Reset the incoming x-flux for the x-hi bc
         fx = psii(n,2)
      END IF

      ! Start the loop in the positive x-direction
      DO i = xs, xe, incx
         x = dx(i)
         m = mat(i)
         sig = sigt(m,g)
         tsxs = sigs(m,:,g,g)

         ! Get spatial parameters
         ex = mu/x
       
         ! Initialize 'a' coefficient
         a = 0.0

         ! Contributions from DD outgoing substitutions
         a = a + 2.0*ex
         ! Contribution from total interaction
         a = a + sig

         ! Initialize rhs coefficient
         b = 0.0

         ! Contribution from incoming fluxes
         b = b + 2.0*ex*fx
         ! Contribution from scattering plus fixed sources
         ! Initialize Legendre polynomials
         plm1 = 0.0
         pl   = 1.0
         ea   = 0.0
         DO l = 0, anord
            ieq = (i-1)*sord + 1 + l
            et(ieq) = tsxs(l)*e(ieq) + sm(ieq,g)
            ea = ea + (2.0*l+1.0)*pl*et(ieq)
            CALL legpoly(l,incx,mu,plm1,pl)
         END DO
         b = b + ea

         ! Solve for the flux
         b = b/a

         ! Update the scalar flux
         plm1 = 0.0
         pl   = 1.0
         DO l = 0, anord
            ieq = (i - 1)*sord + 1 + l
            f(ieq,g) = f(ieq,g) + w(n)*pl*b
            CALL legpoly(l,incx,mu,plm1,pl)
         END DO

         ! Compute the outgoing fluxes with the DD equations
         fx = 2.0*b - fx

         ! Store the new angular fluxes at xe
         IF (i == xe) THEN
            psio(n,xdir,g) = fx
         END IF

      ! End loop over x cells
      END DO
   ! End loop over negative and positive x-directions
   END DO

! End loop over angles
END DO         

RETURN
END SUBROUTINE sweep
