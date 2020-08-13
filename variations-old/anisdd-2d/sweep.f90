SUBROUTINE sweep(g)

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
INTEGER :: xs, xe, ys, ye, incx, incy, nfy
INTEGER :: i, j, m, n, ydir, xdir, indx, jndx, oct
INTEGER :: l, ll, lmm, lpm, fact, ieq, t
REAL*8 :: sig, mu, eta, omega, x, y, a, b, fx, pi, tmp
REAL*8 :: clm, sh, plm2, plm1, pl, plm1m, plm, plmp1, ea, et, p
REAL*8, DIMENSION(0:anord) :: tsxs, al, oal
REAL*8, DIMENSION(nx, 2) :: fy

tmp = 0.0
pi = 2.0*ACOS(tmp)

! Initialize the flux solution to zero
f(:,:,:,g) = 0.0

! Start with loop over all angles
DO n = 1, apo
   ! Set up the angles
   mu  = ang(n,1)
   eta = ang(n,2)

   ! Set the index for starting a new angle
   indx = (n-1)*nx

   ! Set the initial incoming y-flux
   DO i = 1, nx
      jndx = indx + i
      fy(i,1) = psii(jndx,3)
      fy(i,2) = psii(jndx,4)
   END DO

   ! Loop over eta<0 then eta>0
   DO ydir = 1, 2
      IF (ydir == 1) THEN
         ys = ny
         ye = 1
         incy = -1
         omega = pi - omg(n)
      ELSE IF (ydir == 2) THEN
         ys = 1
         ye = ny
         incy = 1
         omega = omg(n)
         ! Reset incoming y-flux for the y-lo bc
         indx = (n-1)*nx
         DO i = 1, nx
            jndx = indx + i
            fy(i,1) = psii(jndx,2)
            fy(i,2) = psii(jndx,1)
         END DO
      END IF

      ! Start the loop in negative y-direction, then do positve y-direction
      DO j = ys, ye, incy
         y = dy(j)
         
         ! Reset the index for the xbc's
         indx = apo*nx + (n-1)*ny

         ! Set the incoming x-flux, depends on ydir
         IF (ydir == 1) THEN
            jndx = indx + j
            fx = psii(jndx,3)
         ELSE IF (ydir == 2) THEN
            jndx = indx + j
            fx = psii(jndx,2)
         END IF

         ! Perform two loops, one in negative x-direction, then positive
         DO xdir = 1, 2
            IF (xdir == 1) THEN
               xs = nx
               xe = 1
               incx = -1
               nfy = 1
               IF (ydir == 1) oct = 3
               IF (ydir == 2) oct = 2
            ELSE IF (xdir == 2) THEN
               xs = 1
               xe = nx
               incx = 1
               nfy = 2
               IF (ydir == 1) oct = 4
               IF (ydir == 2) oct = 1
               ! Reset the incoming x-flux for the x-lo bc
               indx = apo*nx + (n-1)*ny
               IF (ydir == 1) THEN
                  jndx = indx + j
                  fx = psii(jndx,4)
               ELSE IF (ydir == 2) THEN
                  jndx = indx + j
                  fx = psii(jndx,1)
               END IF
            END IF

            ! Start the loop in the negative x-direction
            DO i = xs, xe, incx
               x = dx(i)
               m = mat(i,j)
               sig = sigt(m,g)
               tsxs = sigs(:,m,g,g)

               ! Formula for source index
               t = ((i-1) + (j-1)*nx)*nmom

               ! Get spatial parameters
               ex = mu/x
               ey = eta/y

               ! Initialize 'a' coefficient
               a = 0.0

               ! Contributions from DD outgoing substitutions
               a = a + 2.0*(ex + ey)
               ! Contribution from total interaction
               a = a + sig

               ! Initialize rhs coefficient
               b = 0.0

               ! Contribution from incoming fluxes
               b = b + 2.0*(ex*fx + ey*fy(i,nfy))
               ! Contribution from scattering plus fixed sources
               plm1 = 0.0
               pl   = 1.0
               ea   = 0.0
               DO l = 0, anord
                  al = 0.0      ! Initialize
                  IF (l /= 0) CALL legpoly(l,incx,mu,plm2,plm1,pl)
                  al(0) = pl
                  DO ll = 0, l
                     IF (ll /= 0) THEN
                        plm1m = oal(ll-1)
                        plm   = al(ll-1)
                        CALL aslegf(l,ll,incx,mu,plm1m,plm,plmp1)
                        al(ll) = plmp1
                     END IF
                     lmm = l - ll
                     lpm = l + ll
                     clm = (2.0*l + 1.0)*REAL(fact(lmm))/REAL(fact(lpm))
                     sh = SQRT(clm)*al(ll)*COS(ll*omega)
                     indx = l*(l+1)/2 + ll + 1
                     ieq = indx + t
                     et = tsxs(l)*e(indx,i,j) + sm(ieq,g)
                     p = 1.0
                     IF (ll > 0) p = 2.0
                     ea = ea + p*sh*et
                  END DO
                  ! Reset values
                  plm2 = plm1
                  plm1 = pl
                  oal = al
               END DO
               b = b + ea
               
               ! Solve for the flux
               b = b/a

               ! Update the scalar flux
               plm1 = 0.0
               pl   = 1.0
               ea   = 0.0
               DO l = 0, anord
                  al = 0.0      ! Initialize
                  IF (l /= 0) CALL legpoly(l,incx,mu,plm2,plm1,pl)
                  al(0) = pl
                  DO ll = 0, l
                     IF (ll /= 0) THEN
                        plm1m = oal(ll-1)
                        plm   = al(ll-1)
                        CALL aslegf(l,ll,incx,mu,plm1m,plm,plmp1)
                        al(ll) = plmp1
                     END IF
                     lmm = l - ll
                     lpm = l + ll
                     clm = (2.0*l + 1.0)*REAL(fact(lmm))/REAL(fact(lpm))
                     sh = SQRT(clm)*al(ll)*COS(ll*omega)
                     indx = l*(l+1)/2 + ll + 1
                     f(indx,i,j,g) = f(indx,i,j,g) + w(n)*sh*b
                  END DO
                  ! Reset values
                  plm2 = plm1
                  plm1 = pl
                  oal = al
               END DO

               ! Compute the outgoing fluxes with the DD equations
               fx = 2.0*b - fx
               fy(i,nfy) = 2.0*b - fy(i,nfy)

               ! Store the new angular fluxes at ye xe
               IF (j == ye) THEN
                  indx = (n-1)*nx + i
                  psio(indx,oct,g) = fy(i,nfy)
               END IF
               IF (i == xe) THEN
                  indx = apo*nx + (n-1)*ny + j
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
   
! End loop over angles
END DO         

RETURN
END SUBROUTINE sweep
