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
INTEGER :: xs, xe, ys, ye, zs, ze, incx, incy, incz, ord, nfy, nfz
INTEGER :: i, j, k, t, u, v, m, n, sgm, sge, sgx, ydir, xdir, zdir, mltx, mlty, mltz
INTEGER :: ieq, tt, col, indx, jndx, kndx, info, tmp1, tmp2, oct
INTEGER, DIMENSION(ordcb) :: pv2
REAL*8 :: sig, mu, eta, xi, x, y, z, c, sgn, factor
REAL*8, DIMENSION(ordcb, ordcb) :: a
REAL*8, DIMENSION(ordcb) :: b, wrk
REAL*8, DIMENSION(0:lambda,0:lambda) :: fx
REAL*8, DIMENSION(nx,0:lambda,0:lambda,2) :: fy
REAL*8, DIMENSION(nx,ny,0:lambda,0:lambda,2,2) :: fz
REAL*8, DIMENSION(nx,ny,nz,0:lambda,0:lambda,0:lambda), INTENT(IN) :: e
! Initialize the flux solution to zero
f(:,:,:,:,:,:,g) = 0.0

! Start with loop over all angles
DO n = 1, apo
   ! Set up the angles
   mu  = ang(n,1)
   eta = ang(n,2)
   xi  = ang(n,3)
  
   ! Initialize the incoming z-flux for z-hi bc (O5-8)
   indx = (n-1)*(nx*ny)*ordsq
   DO t = 0, lambda
      DO u = 0, lambda
         DO j = 1, ny
            DO i = 1, nx
               jndx = indx + ((i-1)+(j-1)*nx)*ordsq + t*order + u + 1
               fz(i,j,t,u,1,1) = psii(jndx,7)
               fz(i,j,t,u,2,1) = psii(jndx,8)
               fz(i,j,t,u,1,2) = psii(jndx,6)
               fz(i,j,t,u,2,2) = psii(jndx,5)
            END DO
         END DO
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
         indx = (n-1)*(nx*ny)*ordsq
         DO t = 0, lambda
            DO u = 0, lambda
               DO j = 1, ny
                  DO i = 1, nx
                     jndx = indx + ((i-1)+(j-1)*nx)*ordsq + t*order + u + 1
                     fz(i,j,t,u,1,1) = psii(jndx,3)
                     fz(i,j,t,u,2,1) = psii(jndx,4)
                     fz(i,j,t,u,1,2) = psii(jndx,2)
                     fz(i,j,t,u,2,2) = psii(jndx,1)
                  END DO
               END DO
            END DO
         END DO
      END IF

   ! Start the loop in the negative z-direction, then do positive z-direction
   DO k = zs, ze, incz
      z = dz(k)
 
      ! Reset the index for ybc's
      indx = apo*(nx*ny)*ordsq + (n-1)*(nx*nz)*ordsq
      ! Set the initial incoming y-flux depending on the z-dir
      IF (zdir == 1) THEN
         DO t = 0, lambda
            DO v = 0, lambda
               DO i = 1, nx
                  jndx = indx + ((i-1)+(k-1)*nx)*ordsq + t*order + v + 1
                  fy(i,t,v,1) = psii(jndx,7)
                  fy(i,t,v,2) = psii(jndx,8)
               END DO
            END DO
         END DO
      ELSE IF (zdir == 2) THEN
         DO t = 0, lambda
            DO v = 0, lambda
               DO i = 1, nx
                  jndx = indx + ((i-1)+(k-1)*nx)*ordsq + t*order + v + 1
                  fy(i,t,v,1) = psii(jndx,3)
                  fy(i,t,v,2) = psii(jndx,4)
               END DO
            END DO
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
            indx = apo*(nx*ny)*ordsq + (n-1)*(nx*nz)*ordsq
            IF (zdir == 1) THEN
               DO t = 0, lambda
                  DO v = 0, lambda
                     DO i = 1, nx
                        jndx = indx + ((i-1)+(k-1)*nx)*ordsq + t*order + v + 1
                        fy(i,t,v,1) = psii(jndx,6)
                        fy(i,t,v,2) = psii(jndx,5)
                     END DO
                  END DO
               END DO
            ELSE IF (zdir == 2) THEN
               DO t = 0, lambda
                  DO v = 0, lambda
                     DO i = 1, nx
                        jndx = indx + ((i-1)+(k-1)*nx)*ordsq + t*order + v + 1
                        fy(i,t,v,1) = psii(jndx,2)
                        fy(i,t,v,2) = psii(jndx,1)
                     END DO
                  END DO
               END DO
            END IF
         END IF

      ! Start the loop in negative y-direction, then do positve y-direction
      DO j = ys, ye, incy
         y = dy(j)
         
         ! Reset the index for the xbc's
         indx = apo*(nx*ny+nx*nz)*ordsq + (n-1)*(ny*nz)*ordsq
         ! Set the incoming x-flux with the x-hi bc (octants 2,3,6,7)
         ! Depends on both zdir and ydir
         IF (zdir == 1) THEN
            IF (ydir == 1) THEN
               DO u = 0, lambda
                  DO v = 0, lambda
                     jndx = indx + ((j-1)+(k-1)*ny)*ordsq + u*order + v + 1
                     fx(u,v) = psii(jndx,7)
                  END DO
               END DO
            ELSE IF (ydir == 2) THEN
               DO u = 0, lambda
                  DO v = 0, lambda
                     jndx = indx + ((j-1)+(k-1)*ny)*ordsq + u*order + v + 1
                     fx(u,v) = psii(jndx,6)
                  END DO
               END DO
            END IF
         ELSE IF (zdir == 2) THEN
            IF (ydir == 1) THEN
               DO u = 0, lambda
                  DO v = 0, lambda
                     jndx = indx + ((j-1)+(k-1)*ny)*ordsq + u*order + v + 1
                     fx(u,v) = psii(jndx,3)
                  END DO
               END DO
            ELSE IF (ydir == 2) THEN
               DO u = 0, lambda
                  DO v = 0, lambda
                     jndx = indx + ((j-1)+(k-1)*ny)*ordsq + u*order + v + 1
                     fx(u,v) = psii(jndx,2)
                  END DO
               END DO
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
               indx = apo*(nx*ny+nx*nz)*ordsq + (n-1)*(ny*nz)*ordsq
               IF (zdir == 1) THEN
                  IF (ydir == 1) THEN
                     DO u = 0, lambda
                        DO v = 0, lambda
                           jndx = indx + ((j-1)+(k-1)*ny)*ordsq + u*order + v + 1
                           fx(u,v) = psii(jndx,8)
                        END DO
                     END DO
                  ELSE IF (ydir == 2) THEN
                     DO u = 0, lambda
                        DO v = 0, lambda
                           jndx = indx + ((j-1)+(k-1)*ny)*ordsq + u*order + v + 1
                           fx(u,v) = psii(jndx,5)
                        END DO
                     END DO
                  END IF
               ELSE IF (zdir == 2) THEN
                  IF (ydir == 1) THEN
                     DO u = 0, lambda
                        DO v = 0, lambda
                           jndx = indx + ((j-1)+(k-1)*ny)*ordsq + u*order + v + 1
                           fx(u,v) = psii(jndx,4)
                        END DO
                     END DO
                  ELSE IF (ydir == 2) THEN
                     DO u = 0, lambda
                        DO v = 0, lambda
                           jndx = indx + ((j-1)+(k-1)*ny)*ordsq + u*order + v + 1
                           fx(u,v) = psii(jndx,1)
                        END DO
                     END DO
                  END IF
               END IF
            END IF

         ! Start the loop in the negative x-direction
         DO i = xs, xe, incx
            x = dx(i)
            m = mat(i,j,k)
            sig = sigt(m,g)
            ord = lambda
            c = sigs(m,g,g)/sig     ! Scattering ratio
         
            ! Call for the calculation of the spatial weights
            CALL weight(ord,x,y,z,sig,mu,eta,xi)
            
            ! Begin constructing Matrix Equation
            sgm = incx
            sge = incy
            sgx = incz
            ieq = 0
         
            ! Initialize 'a' matrix
            a = 0.0

            DO t = 0, lambda
               mltx = sgm**t
            
               DO u = 0, lambda
                  mlty = sge**u
                  
                  DO v = 0, lambda
                     mltz = sgx**v
                     ieq = ieq + 1

                     ! Contributions from outgoing fluxes
                     ! Even summations
                     DO tt = 0, lambda, 2
                        ! x-constant surface
                        col = ordsq*tt + order*u + v + 1
                        a(ieq,col) = a(ieq,col) + mltx*(2.0*tt+1.0)/(ex*(1.0+alpha))
                        ! y-constant surface
                        col = ordsq*t + order*tt + v + 1
                        a(ieq,col) = a(ieq,col) + mlty*(2.0*tt+1.0)/(ey*(1.0+beta))
                        ! z-constant surface
                        col = ordsq*t + order*u + tt + 1
                        a(ieq,col) = a(ieq,col) + mltz*(2.0*tt+1.0)/(ez*(1.0+gamma))
                     END DO
               
                     ! Odd summations
                     DO tt = 1, lambda, 2
                        ! x-constant surface
                        col = ordsq*tt + order*u + v + 1
                        a(ieq,col) = a(ieq,col) + mltx*(2.0*tt+1.0)*sgm*alpha/(ex*(1.0+alpha))
                        ! y-constant surface
                        col = ordsq*t + order*tt + v + 1
                        a(ieq,col) = a(ieq,col) + mlty*(2.0*tt+1.0)*sge*beta/(ey*(1.0+beta))
                        ! z-constant surface
                        col = ordsq*t + order*u + tt + 1
                        a(ieq,col) = a(ieq,col) + mltz*(2.0*tt+1.0)*sgx*gamma/(ez*(1.0+gamma))
                     END DO
               
                     ! Contributions from three summations
                     ! x-summations
                     DO tt = MOD((t+1),2), (t-1), 2
                        col = ordsq*tt + order*u + v + 1
                        a(ieq,col) = a(ieq,col) - sgm*(2.0*tt+1.0)/ex
                     END DO
                     ! y-summations
                     DO tt = MOD((u+1),2), (u-1), 2
                        col = ordsq*t + order*tt + v + 1
                        a(ieq,col) = a(ieq,col) - sge*(2.0*tt+1.0)/ey
                     END DO
                     ! z-summations
                     DO tt = MOD((v+1),2), (v-1), 2
                        col = ordsq*t + order*u + tt + 1
                        a(ieq,col) = a(ieq,col) - sgx*(2.0*tt+1.0)/ez
                     END DO
               
                     ! Contribution along the diagonal -- total interaction
                     a(ieq,ieq) = a(ieq,ieq) + 1.0
                     ! Finished calculating the 'a' matrix LHS
         
                     ! Begin composing the RHS vector
                     ! Initially set b to the scattering + fixed source
                     b(ieq) = c*e(i,j,k,t,u,v) + s(i,j,k,t,u,v,g)/sig
                     ! Add contributions from incoming fluxes due to elimination
                     b(ieq) = b(ieq) + mltx*(1.0-alpha)*fx(u,v)/(2.0*ex*(1.0+alpha))
                     b(ieq) = b(ieq) + mlty*(1.0-beta)*fy(i,t,v,nfy)/(2.0*ey*(1.0+beta))
                     b(ieq) = b(ieq) + mltz*(1.0-gamma)*fz(i,j,t,u,nfy,nfz)/(2.0*ez*(1.0+gamma))
                     ! Add contributions from incoming fluxes
                     b(ieq) = b(ieq) + mltx*((-1)**t)*fx(u,v)/(2.0*ex)
                     b(ieq) = b(ieq) + mlty*((-1)**u)*fy(i,t,v,nfy)/(2.0*ey)
                     b(ieq) = b(ieq) + mltz*((-1)**v)*fz(i,j,t,u,nfy,nfz)/(2.0*ez)
                     ! Finished calculating the b vector, RHS
                  END DO
               END DO
            END DO

            ! Make the matrix symmetric
            IF (order /= 1) THEN
               DO indx = 1, order
                  tmp1 = MOD(indx,2)
                  DO jndx = 1, order
                     tmp2 = MOD(jndx,2)
                     IF (tmp1 == tmp2) sgn =  1.0
                     IF (tmp1 /= tmp2) sgn = -1.0
                     DO kndx = 1, order
                        factor = sgn*(2.0*indx-1.0)*(2.0*jndx-1.0)*(2.0*kndx-1.0)
                        ieq = (indx-1)*ordsq + (jndx-1)*order + kndx
                        sgn = -sgn
                        DO col = 1, ordcb
                           a(ieq,col) = factor*a(ieq,col)
                        END DO
                        b(ieq) = factor*b(ieq)
                     END DO
                  END DO
               END DO
            END IF

            ! Solve the matrix
            ! Lapack symmetric solver
            IF (lambda == 0) THEN
               CALL dposv('U',ordcb,1,a,ordcb,b,ordcb,info)
               IF (info /= 0) THEN
                  WRITE (8,'(//,1X,A)') "WARNING: matrix either has illegal value or is singular."
                  warn = warn + 1
                  IF (info > 0) THEN
                     WRITE (8,'(2X,A,/)') "a-Matrix may not be positive definite, use LAPACK routine symmetric indefinite"
                     CALL dsysv('U',ordcb,1,a,ordcb,pv2,b,ordcb,wrk,ordcb,info)
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
               CALL dsysv('U',ordcb,1,a,ordcb,pv2,b,ordcb,wrk,ordcb,info)
               IF (info /= 0) THEN
                  WRITE (8,'(/,1X,A)') "ERROR: matrix either has illegal value or is singular."
                  STOP
               END IF
            END IF

            ! Update the scalar flux solution
            DO t = 0, lambda
               DO u = 0, lambda
                  DO v = 0, lambda
                     indx = ordsq*t + order*u + v + 1
                     f(i,j,k,t,u,v,g) = f(i,j,k,t,u,v,g) + w(n)*b(indx)
                  END DO
               END DO
            END DO
            
            ! Compute the outgoing fluxes with the WDD equations
            ! Outgoing flux moments in x-dir
            DO u = 0, lambda
               DO v = 0, lambda
                  ! Contribution from incoming flux
                  fx(u,v) = -((1.0 - alpha)/(1.0 + alpha))*fx(u,v)
                  ! Contribution from even summation
                  DO tt = 0, lambda, 2
                     indx = ordsq*tt + order*u + v + 1
                     fx(u,v) = fx(u,v) + 2.0*(2.0*tt + 1.0)*b(indx)/(1.0+alpha)
                     ! Contribution from odd summation
                     IF ((tt+1) <= lambda) THEN
                        fx(u,v) = fx(u,v) + 2.0*(2.0*tt+3.0)*sgm*alpha*b(indx+ordsq)/(1.0+alpha)
                     END IF
                  END DO
               END DO
            END DO
            
            ! Outgoing flux moments in y-dir
            DO t = 0, lambda
               DO v = 0, lambda
                  ! Contribution from incoming flux
                  fy(i,t,v,nfy) = -((1.0-beta)/(1.0+beta))*fy(i,t,v,nfy)
                  ! Contribution from even summation
                  DO tt = 0, lambda, 2
                     indx = ordsq*t + order*tt + v + 1
                     fy(i,t,v,nfy) = fy(i,t,v,nfy) + 2.0*(2.0*tt+1.0)*b(indx)/(1.0+beta)
                     ! Contribution from odd summation
                     IF ((tt+1) <= lambda) THEN
                        fy(i,t,v,nfy) = fy(i,t,v,nfy) + 2.0*(2.0*tt+3.0)*sge*beta*b(indx+order)/(1.0+beta)
                     END IF
                  END DO
               END DO
            END DO
         
            ! Outgoing flux moments in z-dir
            DO t = 0, lambda
               DO u = 0, lambda
                  ! Contribution from incoming flux
                  fz(i,j,t,u,nfy,nfz) = -((1.0-gamma)/(1.0+gamma))*fz(i,j,t,u,nfy,nfz)
                  ! Contribution from even summation
                  DO tt = 0, lambda, 2
                     indx = ordsq*t + order*u + tt + 1
                     fz(i,j,t,u,nfy,nfz) = fz(i,j,t,u,nfy,nfz) + 2.0*(2.0*tt+1.0)*b(indx)/(1.0+gamma)
                     ! Contribution from odd summation
                     IF ((tt+1) <= lambda) THEN
                        fz(i,j,t,u,nfy,nfz) = fz(i,j,t,u,nfy,nfz) + 2.0*(2.0*tt+3.0)*sgx*gamma*b(indx+1)/(1.0+gamma)
                     END IF
                  END DO
               END DO
            END DO

            ! Store the new angular fluxes at ze ye xe
            IF (k == ze) THEN
               indx = (n-1)*(nx*ny)*ordsq + ((i-1) + (j-1)*nx)*ordsq
               DO t = 0, lambda
                  DO u = 0, lambda
                     jndx = indx + t*order + u + 1
                     psio(jndx,oct,g) = fz(i,j,t,u,nfy,nfz)
                  END DO
               END DO
            END IF
            IF (j == ye) THEN
               indx = apo*(nx*ny)*ordsq + (n-1)*(nx*nz)*ordsq + ((i-1) + (k-1)*nx)*ordsq
               DO t = 0, lambda
                  DO v = 0, lambda
                     jndx = indx + t*order + v + 1
                     psio(jndx,oct,g) = fy(i,t,v,nfy)
                  END DO
               END DO
            END IF
            IF (i == xe) THEN
               indx = apo*(nx*ny+nx*nz)*ordsq + (n-1)*(ny*nz)*ordsq + ((j-1) + (k-1)*ny)*ordsq
               DO u = 0, lambda
                  DO v = 0, lambda
                     jndx = indx + u*order + v + 1
                     psio(jndx,oct,g) = fx(u,v)
                  END DO
               END DO
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
