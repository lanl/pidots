SUBROUTINE sweep(g)

!-------------------------------------------------------------
!
!  Sweeps across the 2-D matrix
!   Starts at top right corner (mu, eta < 0), then sweeps
!   down all rows, accounting for reflection if necessary. Then
!   sweeps up all rows for eta>0
! 
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER, INTENT(IN) :: g
INTEGER :: xs, xe, ys, ye, incx, incy, ord, nfy, quad
INTEGER :: i, j, k, l, m, n, sgm, sge, ydir, xdir, mltx, mlty
INTEGER :: ieq, ll, col, indx, jndx, info
REAL*8, DIMENSION(ordsq, ordsq) :: a
REAL*8, DIMENSION(ordsq) :: b
REAL*8, DIMENSION(0:lambda) :: fx
REAL*8, DIMENSION(nx, 0:lambda, 2) :: fy
REAL*8 :: sig, mu, eta, x, y, c, sgn, factor
INTEGER, DIMENSION(ordsq) :: piv
REAL*8, DIMENSION(ordsq) :: wrk

! Initialize the flux solution to zero
f(:,:,:,:,g) = 0.0

! Start with loop over all angles
DO n = 1, apo
   ! Set up the angles
   mu  = ang(n,1)
   eta = ang(n,2)
  
   ! Initialize the incoming y-flux for top bc (Q3 and Q4)
   DO k = 0,lambda
      DO i = 1, nx
         indx = (n-1)*nx*order + (i-1)*order + 1 + k
         fy(i,k,1) = psi3(indx)
         fy(i,k,2) = psi4(indx)
      END DO
   END DO
   
   ! Loop over eta<0 then eta>0
   DO ydir = 1, 2
      IF (ydir == 1) THEN
         ys = ny
         ye = 1
         incy = -1
      ELSE IF (ydir == 2) THEN
         ys = 1
         ye = ny
         incy = 1
         ! Reset incoming y-flux for the bottom bc (Q1 and Q2)
         DO k = 0, lambda
            DO i = 1, nx
               indx = (n-1)*nx*order + (i-1)*order + 1 + k
               fy(i,k,1) = psi2(indx)
               fy(i,k,2) = psi1(indx)
            END DO
         END DO
      END IF

   ! Start the loop in negative y-direction, then do positve y-direction
   DO j = ys, ye, incy
      y = dy(j)

      ! Set the incoming x-flux with the right bc depending on ydir
      IF (ydir == 1) THEN
         ! Quadrant 3
         DO l = 0, lambda
            indx = apo*nx*order + (n-1)*ny*order + (j-1)*order + 1 + l
            fx(l) = psi3(indx)
         END DO
      ELSE IF (ydir == 2) THEN
         ! Quadrant 2
         DO l = 0, lambda
            indx = apo*nx*order + (n-1)*ny*order + (j-1)*order + 1 + l
            fx(l) = psi2(indx)
         END DO
      END IF

      ! Perform two loops, one in negative x-direction, then positive
      DO xdir = 1, 2
         IF (xdir == 1) THEN
            xs = nx
            xe = 1
            incx = -1
            nfy = 1
            IF (ydir == 1) quad = 3
            IF (ydir == 2) quad = 2
         ELSE IF (xdir == 2) THEN
            xs = 1
            xe = nx
            incx = 1
            nfy = 2
            IF (ydir == 1) quad = 4
            IF (ydir == 2) quad = 1
            ! Reset the incoming x-flux for the left bc
            IF (ydir == 1) THEN
               ! Quadrant 1
               DO l = 0, lambda
                  indx = apo*nx*order + (n-1)*ny*order + (j-1)*order + 1 + l
                  fx(l) = psi4(indx)
               END DO
            ELSE IF (ydir == 2) THEN
               ! Quadrant 2
               DO l = 0, lambda
                  indx = apo*nx*order + (n-1)*ny*order + (j-1)*order + 1 + l
                  fx(l) = psi1(indx)
               END DO
            END IF
         END IF
         
      ! Start the loop in the negative x-direction
      DO i = xs, xe, incx
         x = dx(i)
         m = mat(i,j)
         sig = sigt(m,g)
         ord = lambda
         c = sigs(m,g,g)/sig     ! Scattering ratio
         
         ! Call for the calculation of the spatial weights
         CALL weight(ord,x,y,sig,mu,eta)
         
         ! Begin constructing Matrix Equation
         sgm = incx
         sge = incy
         ieq = 0
         
         ! Initialize 'a' matrix
         a = 0.0

         DO k = 0, lambda
            mltx = sgm**k
            
            DO l = 0, lambda
               mlty = sge**l
               ieq = ieq + 1

               ! Contributions from outgoing fluxes
               ! Even summations
               DO ll = 0, lambda, 2
                  ! x-constant surface
                  col = order*ll + l + 1
                  a(ieq,col) = a(ieq,col) + mltx*(2.0*ll+1.0)/(ex*(1.0+alpha))
                  ! y-constant surface
                  col = order*k + ll + 1
                  a(ieq,col) = a(ieq,col) + mlty*(2.0*ll+1.0)/(ey*(1.0+beta))
               END DO
               
               ! Odd summations
               DO ll = 1, lambda, 2
                  ! x-constant surface
                  col = order*ll + l + 1
                  a(ieq,col) = a(ieq,col) + mltx*(2.0*ll+1.0)*sgm*alpha/(ex*(1.0+alpha))
                  ! y-constant surface
                  col = order*k + ll + 1
                  a(ieq,col) = a(ieq,col) + mlty*(2.0*ll+1.0)*sge*beta/(ey*(1.0+beta))
               END DO
               
               ! Contributions from two summations
               ! x-summations
               DO ll = MOD((k+1),2), (k-1), 2
                  col = order*ll + l + 1
                  a(ieq,col) = a(ieq,col) - sgm*(2.0*ll+1.0)/ex
               END DO
               ! y-summations
               DO ll = MOD((l+1),2), (l-1), 2
                  col = order*k + ll + 1
                  a(ieq,col) = a(ieq,col) - sge*(2.0*ll+1.0)/ey
               END DO
               
               ! Contribution along the diagonal -- total interaction
               a(ieq,ieq) = a(ieq,ieq) + 1.0
               ! Finished calculating the 'a' matrix LHS
         
               ! Begin composing the RHS vector
               ! Initially set b to the scattering + fixed source
               b(ieq) = c*e(i,j,k,l) + s(i,j,k,l,g)/sig
               ! Add contributions from incoming fluxes due to elimination
               b(ieq) = b(ieq) + mltx*(1.0-alpha)*fx(l)/(2.0*ex*(1.0+alpha))
               b(ieq) = b(ieq) + mlty*(1.0-beta)*fy(i,k,nfy)/(2.0*ey*(1.0+beta))
               ! Add contributions from incoming fluxes
               b(ieq) = b(ieq) + mltx*((-1)**k)*fx(l)/(2.0*ex)
               b(ieq) = b(ieq) + mlty*((-1)**l)*fy(i,k,nfy)/(2.0*ey)
               ! Finished calculating the b vector, RHS
            
            END DO
         END DO

         ! Make the matrix symmetric
         IF (order /= 1) THEN
            DO indx = 1, order
               sgn = 1.0
               jndx = indx/2
               jndx = indx - 2*jndx
               IF (jndx == 0) sgn = -1.0
               DO jndx = 1, order
                  factor = sgn*(2.0*indx-1.0)*(2.0*jndx-1.0)
                  ieq = jndx + (indx - 1)*order
                  sgn = -sgn
                  DO col = 1, ordsq
                     a(ieq,col) = factor*a(ieq,col)
                  END DO
                  b(ieq) = factor*b(ieq)
               END DO
            END DO
         END IF
        
         ! Asymmetric solver
           ! CALL dgesv(ordsq,1,a,ordsq,ipiv,b,ordsq,info)
         ! Symmetric solver
         ! Need to use different lapack solver for lambda>0 because of positive
         ! definite problem
         IF (lambda == 0) THEN
            CALL dposv('U',ordsq,1,a,ordsq,b,ordsq,info)
            IF (info /= 0) THEN
               WRITE (8,'(//,1X,A)') "WARNING: matrix either has illegal value or is singular."
               warn = warn + 1
               IF (info > 0) THEN
                  WRITE (8,'(2X,A,/)') "a-Matrix may not be positive definite, use LAPACK routine for symmetric indefinite."
                  CALL dsysv('U',ordsq,1,a,ordsq,piv,b,ordsq,wrk,ordsq,info)
                  IF (info /= 0) THEN
                     WRITE (8,'(/,1X,A)') "ERROR: Unable to solve system of equations in sweep for current cell."
                     STOP
                  END IF
               ELSE
                  WRITE (8,'(//,1X,A)') "ERROR: matrix has unresolved error and problem not solved."
                  STOP
               END IF
            END IF
         ELSE IF (lambda > 0) THEN
            CALL dsysv('U',ordsq,1,a,ordsq,piv,b,ordsq,wrk,ordsq,info)
            IF (info /= 0) THEN
               WRITE (8,'(/,1X,A)') "ERROR: matrix either has illegal value or is singular."
               STOP
            END IF
         END IF
 
         ! Update the scalar flux solution
         DO k = 0, lambda
            DO l = 0, lambda
               indx = order*k + l + 1
               f(i,j,k,l,g) = f(i,j,k,l,g) + w(n)*b(indx)
            END DO
         END DO
         
         ! Compute the outgoing fluxes with the WDD equations
         ! Outgoing flux moments in x-dir
         DO l = 0, lambda
            ! Contribution from incoming flux
            fx(l) = -((1.0 - alpha)/(1.0 + alpha))*fx(l)
            ! Contribution from even summation
            DO ll = 0, lambda, 2
               indx = order*ll + l + 1
               fx(l) = fx(l) + 2.0*(2.0*ll + 1.0)*b(indx)/(1.0+alpha)
               ! Contribution from odd summation
               IF ((ll+1) <= lambda) THEN
                  fx(l) = fx(l) + 2.0*(2.0*ll+3.0)*sgm*alpha*b(indx+order)/(1.0+alpha)
               END IF
            END DO
         END DO
         
         ! Outgoing flux moments in y-dir
         DO k = 0, lambda
            ! Contribution from incoming flux
            fy(i,k,nfy) = -((1.0-beta)/(1.0+beta))*fy(i,k,nfy)
            ! Contribution from even summation
            DO ll = 0, lambda, 2
               indx = order*k + ll + 1
               fy(i,k,nfy) = fy(i,k,nfy) + 2.0*(2.0*ll+1.0)*b(indx)/(1.0+beta)
               ! Contribution from odd summation
               IF ((ll+1) <= lambda) THEN
                  fy(i,k,nfy) = fy(i,k,nfy) + 2.0*(2.0*ll+3.0)*sge*beta*b(indx+1)/(1.0+beta)
               END IF
            END DO
         END DO

         ! Store the new angular fluxes at ye and xe
         IF (j == ye) THEN
            indx = (n-1)*nx*order + (i-1)*order
            psio((indx+1):(indx+order),quad,g) = fy(i,:,nfy)
         END IF
         IF (i == xe) THEN
            indx = apo*nx*order + (n-1)*ny*order + (j-1)*order
            psio((indx+1):(indx+order),quad,g) = fx(:)
         END IF

      ! End loop over x cells
      END DO
      ! End loop over negative and positive x-directions
      END DO
   
   ! End loop over y cells
   END DO
   
   ! End loop over negative and positve y-directions
   END DO

! End loop over angles
END DO      

RETURN
END SUBROUTINE sweep
