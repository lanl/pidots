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
INTEGER :: xs, xe, ys, ye, zs, ze, incx, incy, incz, nfy, nfz
INTEGER :: i, j, k, m, n, ydir, xdir, zdir, indx, jndx, oct
INTEGER :: l, ll, lmm, lpm, fact, ieq, t
REAL*8 :: sig, mu, eta, xi, omega, x, y, z, a, b, fx, pi, tmp
REAL*8 :: clm, she, sho, plm2, plm1, pl, plm1m, plm, plmp1, ea, et
REAL*8, DIMENSION(0:anord) :: tsxs, al, oal
REAL*8, DIMENSION(nx,2) :: fy
REAL*8, DIMENSION(nx,ny,2,2) :: fz

tmp = 0.0
pi = 2.0*ACOS(tmp)

! Initialize the flux solution to zero
f(:,:,:,:,g) = 0.0

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
            IF (incz == 1) THEN
               omega = pi - omg(n)
            ELSE
               omega = pi + omg(n)
            END IF 
         ELSE IF (ydir == 2) THEN
            ys = 1
            ye = ny
            incy = 1
            nfz = 2
            IF (incz == 1) THEN
               omega = omg(n)
            ELSE
               omega = -omg(n)
            END IF
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
            tsxs = sigs(:,m,g,g)

            ! Formula for source index
            t = ((i-1) + (j-1)*nx + (k-1)*nx*ny)*nmom

            ! Get spatial parameters
            ex = mu/x
            ey = eta/y
            ez = xi/z
         
            ! Initialize 'a' coefficient
            a = 0.0

            ! Contributions from DD outgoing substitutions
            a = a + 2.0*(ex + ey + ez)
            ! Contribution from total interaction
            a = a + sig

            ! Initialize rhs coefficient
            b = 0.0

            ! Contribution from incoming fluxes
            b = b + 2.0*(ex*fx + ey*fy(i,nfy) + ez*fz(i,j,nfy,nfz))
            ! Contribution from scattering plus fixed sources
            plm1 = 0.0
            pl = 1.0
            ea = 0.0
            DO l = 0, anord
               al = 0.0      ! Initialize
               IF (l /= 0) CALL legpoly(l,incx,mu,plm2,plm1,pl)
               al(0) = pl
               indx = l**2 + 1
               ieq = indx + t
               et = tsxs(l)*e(indx,i,j,k) + sm(ieq,g)
               ea = ea + SQRT(2.0*l+1.0)*pl*et
               DO ll = 1, l
                  plm1m = oal(ll-1)
                  plm   = al(ll-1)
                  CALL aslegf(l,ll,incx,mu,plm1m,plm,plmp1)
                  al(ll) = plmp1
                  lmm = l - ll
                  lpm = l + ll
                  clm = (2.0*l + 1.0)*REAL(fact(lmm))/REAL(fact(lpm))
                  she = 2.0*SQRT(clm)*al(ll)*COS(ll*omega)
                  sho = 2.0*SQRT(clm)*al(ll)*SIN(ll*omega)
                  indx = l**2 + 2*ll
                  ieq = indx + t
                  et = tsxs(l)*e(indx,i,j,k) + sm(ieq,g)
                  ea = ea + she*et
                  et = tsxs(l)*e(indx+1,i,j,k) + sm(ieq+1,g)
                  ea = ea + sho*et
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
            pl = 1.0
            DO l = 0, anord
               al = 0.0      ! Initialize
               IF (l /= 0) CALL legpoly(l,incx,mu,plm2,plm1,pl)
               al(0) = pl
               indx = l**2 + 1
               f(indx,i,j,k,g) = f(indx,i,j,k,g) + w(n)*SQRT(2.0*l+1.0)*pl*b
               DO ll = 1, l
                  plm1m = oal(ll-1)
                  plm   = al(ll-1)
                  CALL aslegf(l,ll,incx,mu,plm1m,plm,plmp1)
                  al(ll) = plmp1
                  lmm = l - ll
                  lpm = l + ll
                  clm = (2.0*l + 1.0)*REAL(fact(lmm))/REAL(fact(lpm))
                  she = SQRT(clm)*al(ll)*COS(ll*omega)
                  sho = SQRT(clm)*al(ll)*SIN(ll*omega)
                  indx = l**2 + 2*ll
                  f(indx,i,j,k,g) = f(indx,i,j,k,g) + w(n)*she*b
                  f(indx+1,i,j,k,g) = f(indx+1,i,j,k,g) + w(n)*sho*b
               END DO
               ! Reset values
               plm2 = plm1
               plm1 = pl
               oal = al
            END DO
            
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
