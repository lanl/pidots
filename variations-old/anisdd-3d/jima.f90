SUBROUTINE jima(i,j,k,xs,ys,zs,xe,ye,ze,incx,incy,incz,n,oct,mu,omega,tsxs)

!-------------------------------------------------------------
!
!  Jacobian Iteration MAtrix
!  Directs the equations to set up the matrices that work
!   with the phi and source vectors. takes in parameters from
!   matsweep and solves for the updates to jmat and jpsi 
!
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER, INTENT(IN) :: i, j, k, xs, ys, zs, xe, ye, ze, incx, incy, incz, n, oct
INTEGER :: ii, jj, kk, ieq, jeq, indx, fact, lmm, lpm, l, ll
REAL*8, INTENT(IN) :: mu, omega
REAL*8 :: wt, clm, she, sho, plm2, plm1, pl, plm1m, plm, plmp1
REAL*8, DIMENSION(0:anord), INTENT(IN) :: tsxs
REAL*8, DIMENSION(0:anord) :: al, oal

! Set the quad weight
wt = w(n)

! Save the X, Y, and Z matrices for updates
IF (i /= xs) THEN
   xold = xmat
END IF
IF (j /= ys) THEN
   yold(:,:,:,:) = ymat(:,:,:,:,i)
END IF
IF (k /= zs) THEN
   zold(:,:,:,:) = zmat(:,:,:,:,i,j)
END IF

! Older version without angular flux out did not compute anything at xe, ye, ze
! Now compute that info and use it to construct jpsi

! Update X, Y, Z matrices with old X matrix
IF (i /= xs) THEN
   DO kk = zs, k, incz
      DO jj = ys, j, incy
         DO ii = xs, i, incx
            xmat(:,ii,jj,kk) = gyzyz*xold(:,ii,jj,kk)
            ymat(:,ii,jj,kk,i) = gxzyz*xold(:,ii,jj,kk)
            zmat(:,ii,jj,kk,i,j) = gxyyz*xold(:,ii,jj,kk)
         END DO
      END DO
   END DO
END IF

! Update X, Y, and Z matrices with old Y matrix
IF (j /= ys) THEN
   DO kk = zs, k, incz
      DO jj = ys, j, incy
         DO ii = xs, i, incx
            IF (i /= xs) THEN
               xmat(:,ii,jj,kk) = xmat(:,ii,jj,kk) + gyzxz*yold(:,ii,jj,kk)
               ymat(:,ii,jj,kk,i) = ymat(:,ii,jj,kk,i) + gxzxz*yold(:,ii,jj,kk)
               zmat(:,ii,jj,kk,i,j) = zmat(:,ii,jj,kk,i,j) + gxyxz*yold(:,ii,jj,kk)
            ELSE IF (i == xs) THEN
               xmat(:,ii,jj,kk) = gyzxz*yold(:,ii,jj,kk)
               ymat(:,ii,jj,kk,i) = gxzxz*yold(:,ii,jj,kk)
               zmat(:,ii,jj,kk,i,j) = gxyxz*yold(:,ii,jj,kk)
            END IF
         END DO
      END DO
   END DO
END IF

! Update X, Y, and Z matrices with old Z matrix
IF (k /= zs) THEN
   DO kk = zs, k, incz
      DO jj = ys, j, incy
         DO ii = xs, i, incx
            IF (i /= xs .OR. j /= ys) THEN
               xmat(:,ii,jj,kk) = xmat(:,ii,jj,kk) + gyzxy*zold(:,ii,jj,kk)
               ymat(:,ii,jj,kk,i) = ymat(:,ii,jj,kk,i) + gxzxy*zold(:,ii,jj,kk)
               zmat(:,ii,jj,kk,i,j) = zmat(:,ii,jj,kk,i,j) + gxyxy*zold(:,ii,jj,kk)
            ELSE IF (i == xs .AND. j == ys) THEN
               xmat(:,ii,jj,kk) = gyzxy*zold(:,ii,jj,kk)
               ymat(:,ii,jj,kk,i) = gxzxy*zold(:,ii,jj,kk)
               zmat(:,ii,jj,kk,i,j) = gxyxy*zold(:,ii,jj,kk)
            END IF
         END DO
      END DO
   END DO
END IF

! Append X, Y, and Z matrices
DO l = 0, anord
   indx = l**2 + 1
   ! l0 moments
   zmat(indx,i,j,k,i,j) = gxya(indx)*tsxs(l)
   ymat(indx,i,j,k,i) = gxza(indx)*tsxs(l)
   xmat(indx,i,j,k) = gyza(indx)*tsxs(l)
   ! lm moments
   DO ll = 1, l
      indx = l**2 + 2*ll
      zmat(indx,i,j,k,i,j) = gxya(indx)*tsxs(l)
      zmat(indx+1,i,j,k,i,j) = gxya(indx+1)*tsxs(l)
      ymat(indx,i,j,k,i) = gxza(indx)*tsxs(l)
      ymat(indx+1,i,j,k,i) = gxza(indx+1)*tsxs(l)
      xmat(indx,i,j,k) = gyza(indx)*tsxs(l)
      xmat(indx+1,i,j,k) = gyza(indx+1)*tsxs(l)
   END DO
END DO

! Reset gaa as gaa*sigma_s before updating jmat
DO l = 0, anord
   indx = l**2 + 1
   gaa(indx) = gaa(indx)*tsxs(l)
   DO ll = 1, l
      indx = l**2 + 2*ll
      gaa(indx) = gaa(indx)*tsxs(l)
      gaa(indx+1) = gaa(indx+1)*tsxs(l)
   END DO
END DO

!------------------------------------------------------------------------
! Start to formulate Jacobian matrix, diagonal blocks
! Depends on the 'tpose' flag how the indices are incremented
! tpose = 1 --> make transposed Jphi and Jpsi
IF (tpose == 1) THEN
   jeq = ((i-1) + (j-1)*nx + (k-1)*nx*ny)*nmom
   ! Diagonal block
   plm1 = 0.0
   pl = 1.0
   DO l = 0, anord
      al = 0.0      ! Initialize
      IF (l /= 0) CALL legpoly(l,incx,mu,plm2,plm1,pl)
      al(0) = pl
      indx = l**2 + 1
      jmat(jeq+1:jeq+nmom,jeq+indx) = jmat(jeq+1:jeq+nmom,jeq+indx) + wt*SQRT(2.0*l+1.0)*pl*gaa
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
         jmat(jeq+1:jeq+nmom,jeq+indx) = jmat(jeq+1:jeq+nmom,jeq+indx) + wt*she*gaa
         jmat(jeq+1:jeq+nmom,jeq+indx+1) = jmat(jeq+1:jeq+nmom,jeq+indx+1) + wt*sho*gaa
      END DO
      ! Reset values
      plm2 = plm1
      plm1 = pl
      oal = al
   END DO

   DO kk = zs, k, incz
      DO jj = ys, j, incy
         DO ii = xs, i, incx
            ieq = ((ii-1) + (jj-1)*nx + (kk-1)*nx*ny)*nmom
            plm1 = 0.0
            pl = 1.0
            DO l = 0, anord
               al = 0.0      ! Initialize
               IF (l /= 0) CALL legpoly(l,incx,mu,plm2,plm1,pl)
               al(0) = pl
               indx = l**2 + 1
               IF (i /= xs) jmat(ieq+1:ieq+nmom,jeq+indx) = jmat(ieq+1:ieq+nmom,jeq+indx) &
                                                          + wt*SQRT(2.0*l+1.0)*pl*gayz*xold(:,ii,jj,kk)
               IF (j /= ys) jmat(ieq+1:ieq+nmom,jeq+indx) = jmat(ieq+1:ieq+nmom,jeq+indx) &
                                                          + wt*SQRT(2.0*l+1.0)*pl*gaxz*yold(:,ii,jj,kk)
               IF (k /= zs) jmat(ieq+1:ieq+nmom,jeq+indx) = jmat(ieq+1:ieq+nmom,jeq+indx) &
                                                          + wt*SQRT(2.0*l+1.0)*pl*gaxy*zold(:,ii,jj,kk)
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
                  IF (i /= xs) THEN
                     jmat(ieq+1:ieq+nmom,jeq+indx) = jmat(ieq+1:ieq+nmom,jeq+indx) + wt*she*gayz*xold(:,ii,jj,kk)
                     jmat(ieq+1:ieq+nmom,jeq+indx+1) = jmat(ieq+1:ieq+nmom,jeq+indx+1) + wt*sho*gayz*xold(:,ii,jj,kk)
                  END IF
                  IF (j /= ys) THEN
                     jmat(ieq+1:ieq+nmom,jeq+indx) = jmat(ieq+1:ieq+nmom,jeq+indx) + wt*she*gaxz*yold(:,ii,jj,kk)
                     jmat(ieq+1:ieq+nmom,jeq+indx+1) = jmat(ieq+1:ieq+nmom,jeq+indx+1) + wt*sho*gaxz*yold(:,ii,jj,kk)
                  END IF
                  IF (k /= zs) THEN
                     jmat(ieq+1:ieq+nmom,jeq+indx) = jmat(ieq+1:ieq+nmom,jeq+indx) + wt*she*gaxy*zold(:,ii,jj,kk)
                     jmat(ieq+1:ieq+nmom,jeq+indx+1) = jmat(ieq+1:ieq+nmom,jeq+indx+1) + wt*sho*gaxy*zold(:,ii,jj,kk)
                  END IF
               END DO
               plm2 = plm1
               plm1 = pl
               oal = al
            END DO
         END DO
      END DO
   END DO

   ! Now check to see if the computational cell is at end of X or Y dimension,
   ! Place an xmat or ymat set of values if it is
   ! First get the z-direction coefficients
   IF (k == ze) THEN
      ! Reset equation index
      jeq = (n-1)*(nx*ny) + i + (j-1)*nx
      DO kk = zs, k, incz
         DO jj = ys, j, incy
            DO ii = xs, i, incx
               ieq = ((ii-1) + (jj-1)*nx + (kk-1)*nx*ny)*nmom
               jpsi(ieq+1:ieq+nmom,jeq,oct) = zmat(:,ii,jj,kk,i,j)
            END DO
         END DO
      END DO
   END IF
   IF (j == ye) THEN
      ! Reset equation index
      jeq = apo*(nx*ny) + (n-1)*(nx*nz) + i + (k-1)*nx
      DO kk = zs, k, incz
         DO jj = ys, j, incy
            DO ii = xs, i, incx
               ieq = ((ii-1) + (jj-1)*nx + (kk-1)*nx*ny)*nmom
               jpsi(ieq+1:ieq+nmom,jeq,oct) = ymat(:,ii,jj,kk,i)
            END DO
         END DO
      END DO
   END IF
   IF (i == xe) THEN
      ! Reset equation index
      jeq = apo*(nx*ny+nx*nz) + (n-1)*(ny*nz) + j + (k-1)*ny
      DO kk = zs, k, incz
         DO jj = ys, j, incy
            DO ii = xs, i, incx
               ieq = ((ii-1) + (jj-1)*nx + (kk-1)*nx*ny)*nmom
               jpsi(ieq+1:ieq+nmom,jeq,oct) = xmat(:,ii,jj,kk)
            END DO
         END DO
      END DO
   END IF
!-----------------------------------------------------------------------
! If tpose = 0 --> make normal Jphi and Jpsi
ELSE
   ieq = ((i-1) + (j-1)*nx + (k-1)*nx*ny)*nmom
   ! Diagonal block
   plm1 = 0.0
   pl = 1.0
   DO l = 0, anord
      al = 0.0      ! Initialize
      IF (l /= 0) CALL legpoly(l,incx,mu,plm2,plm1,pl)
      al(0) = pl
      indx = l**2 + 1
      jmat(ieq+indx,ieq+1:ieq+nmom) = jmat(ieq+indx,ieq+1:ieq+nmom) + wt*SQRT(2.0*l+1.0)*pl*gaa
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
         jmat(ieq+indx,ieq+1:ieq+nmom) = jmat(ieq+indx,ieq+1:ieq+nmom) + wt*she*gaa
         jmat(ieq+indx+1,ieq+1:ieq+nmom) = jmat(ieq+indx+1,ieq+1:ieq+nmom) + wt*sho*gaa
      END DO
      ! Reset values
      plm2 = plm1
      plm1 = pl
      oal = al
   END DO

   DO kk = zs, k, incz
      DO jj = ys, j, incy
         DO ii = xs, i, incx
            jeq = ((ii-1) + (jj-1)*nx + (kk-1)*nx*ny)*nmom
            plm1 = 0.0
            pl = 1.0
            DO l = 0, anord
               al = 0.0      ! Initialize
               IF (l /= 0) CALL legpoly(l,incx,mu,plm2,plm1,pl)
               al(0) = pl
               indx = l**2 + 1
               IF (i /= xs) jmat(ieq+indx,jeq+1:jeq+nmom) = jmat(ieq+indx,jeq+1:jeq+nmom) &
                                                          + wt*SQRT(2.0*l+1.0)*pl*gayz*xold(:,ii,jj,kk)
               IF (j /= ys) jmat(ieq+indx,jeq+1:jeq+nmom) = jmat(ieq+indx,jeq+1:jeq+nmom) &
                                                          + wt*SQRT(2.0*l+1.0)*pl*gaxz*yold(:,ii,jj,kk)
               IF (k /= zs) jmat(ieq+indx,jeq+1:jeq+nmom) = jmat(ieq+indx,jeq+1:jeq+nmom) &
                                                          + wt*SQRT(2.0*l+1.0)*pl*gaxy*zold(:,ii,jj,kk)
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
                  IF (i /= xs) THEN
                     jmat(ieq+indx,jeq+1:jeq+nmom) = jmat(ieq+indx,jeq+1:jeq+nmom) + wt*she*gayz*xold(:,ii,jj,kk)
                     jmat(ieq+indx+1,jeq+1:jeq+nmom) = jmat(ieq+indx+1,jeq+1:jeq+nmom) + wt*sho*gayz*xold(:,ii,jj,kk)
                  END IF
                  IF (j /= ys) THEN
                     jmat(ieq+indx,jeq+1:jeq+nmom) = jmat(ieq+indx,jeq+1:jeq+nmom) + wt*she*gaxz*yold(:,ii,jj,kk)
                     jmat(ieq+indx+1,jeq+1:jeq+nmom) = jmat(ieq+indx+1,jeq+1:jeq+nmom) + wt*sho*gaxz*yold(:,ii,jj,kk)
                  END IF
                  IF (k /= zs) THEN
                     jmat(ieq+indx,jeq+1:jeq+nmom) = jmat(ieq+indx,jeq+1:jeq+nmom) + wt*she*gaxy*zold(:,ii,jj,kk)
                     jmat(ieq+indx+1,jeq+1:jeq+nmom) = jmat(ieq+indx+1,jeq+1:jeq+nmom) + wt*sho*gaxy*zold(:,ii,jj,kk)
                  END IF
               END DO
               plm2 = plm1
               plm1 = pl
               oal = al
            END DO
         END DO
      END DO
   END DO

   ! Now check to see if the computational cell is at end of X or Y dimension,
   ! Place an xmat or ymat set of values if it is
   ! First get the z-direction coefficients
   IF (k == ze) THEN
      ! Reset equation index
      ieq = (n-1)*(nx*ny) + i + (j-1)*nx
      DO kk = zs, k, incz
         DO jj = ys, j, incy
            DO ii = xs, i, incx
               jeq = ((ii-1) + (jj-1)*nx + (kk-1)*nx*ny)*nmom
               jpsi(ieq,jeq+1:jeq+nmom,oct) = zmat(:,ii,jj,kk,i,j)
            END DO
         END DO
      END DO
   END IF
   IF (j == ye) THEN
      ! Reset equation index
      ieq = apo*(nx*ny) + (n-1)*(nx*nz) + i + (k-1)*nx
      DO kk = zs, k, incz
         DO jj = ys, j, incy
            DO ii = xs, i, incx
               jeq = ((ii-1) + (jj-1)*nx + (kk-1)*nx*ny)*nmom
               jpsi(ieq,jeq+1:jeq+nmom,oct) = ymat(:,ii,jj,kk,i)
            END DO
         END DO
      END DO
   END IF
   IF (i == xe) THEN
      ! Reset equation index
      ieq = apo*(nx*ny+nx*nz) + (n-1)*(ny*nz) + j + (k-1)*ny
      DO kk = zs, k, incz
         DO jj = ys, j, incy
            DO ii = xs, i, incx
               jeq = ((ii-1) + (jj-1)*nx + (kk-1)*nx*ny)*nmom
               jpsi(ieq,jeq+1:jeq+nmom,oct) = xmat(:,ii,jj,kk)
            END DO
         END DO
      END DO
   END IF
END IF

RETURN
END SUBROUTINE jima
