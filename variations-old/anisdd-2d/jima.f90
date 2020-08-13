SUBROUTINE jima(i,j,xs,ys,xe,ye,incx,incy,n,oct,mu,omega,tsxs)

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
INTEGER, INTENT(IN) :: i, j, xs, ys, xe, ye, incx, incy, n, oct
INTEGER :: ii, jj, ieq, jeq, indx, fact, lmm, lpm, l, ll
REAL*8, INTENT(IN) :: mu, omega
REAL*8, DIMENSION(0:anord), INTENT(IN) :: tsxs
REAL*8 :: wt, clm, sh, plm2, plm1, pl, plm1m, plm, plmp1
REAL*8, DIMENSION(0:anord) :: al, oal

! Set the quad weight
wt = w(n)

! Save the X, Y, and Z matrices for updates
IF (i /= xs) THEN
   xold = xmat
END IF
IF (j /= ys) THEN
   yold(:,:,:) = ymat(:,:,:,i)
END IF

! Older version without angular flux out did not compute anything at xe, ye, ze
! Now compute that info and use it to construct jpsi

! Update X, Y matrices with old X matrix
IF (i /= xs) THEN
   DO jj = ys, j, incy
      DO ii = xs, i, incx
         xmat(:,ii,jj) = gyy*xold(:,ii,jj)
         ymat(:,ii,jj,i) = gxy*xold(:,ii,jj)
      END DO
   END DO
END IF

! Update X, Y matrices with old Y matrix
IF (j /= ys) THEN
   DO jj = ys, j, incy
      DO ii = xs, i, incx
         IF (i /= xs) THEN
            xmat(:,ii,jj) = xmat(:,ii,jj) + gyx*yold(:,ii,jj)
            ymat(:,ii,jj,i) = ymat(:,ii,jj,i) + gxx*yold(:,ii,jj)
         ELSE IF (i == xs) THEN
            xmat(:,ii,jj) = gyx*yold(:,ii,jj)
            ymat(:,ii,jj,i) = gxx*yold(:,ii,jj)
         END IF
      END DO
   END DO
END IF

! Append X, Y matrices
DO l = 0, anord
   DO ll = 0, l
      indx = l*(l+1)/2 + ll + 1
      ymat(indx,i,j,i) = gxa(indx)*tsxs(l)
      xmat(indx,i,j)   = gya(indx)*tsxs(l)
   END DO
END DO

! Reset gaa as gaa*sigma_s before updating jmat
DO l = 0, anord
   DO ll = 0, l
      indx = l*(l+1)/2 + ll + 1
      gaa(indx) = gaa(indx)*tsxs(l)
   END DO
END DO

!-------------------------------------------------------------------------------------
! Start to formulate Jacobian matrix, diagonal blocks
! Depends on the 'tpose' flag how the indices are incremented
! tpose = 1 --> make transposed Jphi and Jpsi
IF (tpose == 1) THEN
   jeq = ((i-1) + (j-1)*nx)*nmom
   ! Diagonal Block
   plm1 = 0.0
   pl   = 1.0
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
         jmat(jeq+1:jeq+nmom,jeq+indx) = jmat(jeq+1:jeq+nmom,jeq+indx) + wt*sh*gaa
      END DO
      ! Reset values
      plm2 = plm1
      plm1 = pl
      oal = al
   END DO

   ! Off-diagonal contribution from X
   IF (i /= xs) THEN
      DO jj = ys, j, incy
         DO ii = xs, i, incx
            ieq = ((ii-1) + (jj-1)*nx)*nmom
            plm1 = 0.0
            pl   = 1.0
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
                  jmat(ieq+1:ieq+nmom,jeq+indx) = jmat(ieq+1:ieq+nmom,jeq+indx) + wt*sh*gay*xold(:,ii,jj)
               END DO
               ! Reset values
               plm2 = plm1
               plm1 = pl
               oal = al
            END DO
         END DO
      END DO
   END IF

   ! Off-diagonal contribution from Y
   IF (j /= ys) THEN
      DO jj = ys, j, incy
         DO ii = xs, i, incx
            ieq = ((ii-1) + (jj-1)*nx)*nmom
            plm1 = 0.0
            pl   = 1.0
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
                  jmat(ieq+1:ieq+nmom,jeq+indx) = jmat(ieq+1:ieq+nmom,jeq+indx) + wt*sh*gax*yold(:,ii,jj)
               END DO
               ! Reset values
               plm2 = plm1
               plm1 = pl
               oal = al
            END DO
         END DO
      END DO
   END IF

   ! Now check to see if the computational cell is at end of X or Y dimension,
   ! Place an xmat or ymat set of values if it is
   IF (j == ye) THEN
      ! Reset equation index
      jeq = (n-1)*nx + i
      DO jj = ys, j, incy
         DO ii = xs, i, incx
            ieq = ((ii-1) + (jj-1)*nx)*nmom
            jpsi(ieq+1:ieq+nmom,jeq,oct) = ymat(:,ii,jj,i)
         END DO
      END DO
   END IF
   IF (i == xe) THEN
      ! Reset equation index
      jeq = apo*nx + (n-1)*ny + j
      DO jj = ys, j, incy
         DO ii = xs, i, incx
            ieq = ((ii-1) + (jj-1)*nx)*nmom
            jpsi(ieq+1:ieq+nmom,jeq,oct) = xmat(:,ii,jj)
         END DO
      END DO
   END IF
!-----------------------------------------------------------------------
! If tpose = 0 --> make normal Jphi and Jpsi
ELSE
   ieq = ((i-1) + (j-1)*nx)*nmom
   ! Diagonal Block
   plm1 = 0.0
   pl   = 1.0
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
         jmat(ieq+indx,ieq+1:ieq+nmom) = jmat(ieq+indx,ieq+1:ieq+nmom) + wt*sh*gaa
      END DO
      ! Reset values
      plm2 = plm1
      plm1 = pl
      oal = al
   END DO

   ! Off-diagonal contribution from X
   IF (i /= xs) THEN
      DO jj = ys, j, incy
         DO ii = xs, i, incx
            jeq = ((ii-1) + (jj-1)*nx)*nmom
            plm1 = 0.0
            pl   = 1.0
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
                  jmat(ieq+indx,jeq+1:jeq+nmom) = jmat(ieq+indx,jeq+1:jeq+nmom) + wt*sh*gay*xold(:,ii,jj)
               END DO
               ! Reset values
               plm2 = plm1
               plm1 = pl
               oal = al
            END DO
         END DO
      END DO
   END IF

   ! Off-diagonal contribution from Y
   IF (j /= ys) THEN
      DO jj = ys, j, incy
         DO ii = xs, i, incx
            jeq = ((ii-1) + (jj-1)*nx)*nmom
            plm1 = 0.0
            pl   = 1.0
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
                  jmat(ieq+indx,jeq+1:jeq+nmom) = jmat(ieq+indx,jeq+1:jeq+nmom) + wt*sh*gax*yold(:,ii,jj)
               END DO
               ! Reset values
               plm2 = plm1
               plm1 = pl
               oal = al
            END DO
         END DO
      END DO
   END IF

   ! Now check to see if the computational cell is at end of X or Y dimension,
   ! Place an xmat or ymat set of values if it is
   IF (j == ye) THEN
      ! Reset equation index
      ieq = (n-1)*nx + i
      DO jj = ys, j, incy
         DO ii = xs, i, incx
            jeq = ((ii-1) + (jj-1)*nx)*nmom
            jpsi(ieq,jeq+1:jeq+nmom,oct) = ymat(:,ii,jj,i)
         END DO
      END DO
   END IF
   IF (i == xe) THEN
      ! Reset equation index
      ieq = apo*nx + (n-1)*ny + j
      DO jj = ys, j, incy
         DO ii = xs, i, incx
            jeq = ((ii-1) + (jj-1)*nx)*nmom
            jpsi(ieq,jeq+1:jeq+nmom,oct) = xmat(:,ii,jj)
         END DO
      END DO
   END IF
END IF

RETURN
END SUBROUTINE jima
