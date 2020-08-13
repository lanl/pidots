SUBROUTINE conkt(n,i,j,incx,incy,xs,ys,mu,omega,ktmpy,ktmpx)

!-----------------------------------------------------------
!
! Constructs the ktmp row for a given octant and cell and 
!  returns the values to afcm.
!
! Takes in from afcm: All the computed values not in solvar
!     n,i,j,k,incx,incy,incz,xs,ys,zs,bcs
! Returns to afcm a temporary matrix ktmp that fills in the
!  particular ktmp#
!
!------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER, INTENT(IN) :: n, i, j
INTEGER, INTENT(IN) :: incx, incy, xs, ys
INTEGER :: ii, jj, indx, fact, lmm, lpm, l, ll
REAL*8, INTENT(IN) :: mu, omega
REAL*8 :: wt, clm, sh, plm2, plm1, pl, plm1m, plm, plmp1
REAL*8, DIMENSION(nx,nmom), INTENT(OUT) :: ktmpy
REAL*8, DIMENSION(ny,nmom), INTENT(OUT) :: ktmpx
REAL*8, DIMENSION(0:anord) :: al, oal

! Initialize ktmps
ktmpy = 0.0
ktmpx = 0.0

! Set the weight
wt = w(n)

! All the bc matrices and gamma matrices needed come from solvar

! Perform the operations to update from the different BCs
! Update values from ybc's
IF (j /= ys) THEN
   DO ii = xs, i-incx, incx
      ktmpy(ii,1) = ktmpy(ii,1) + wt*gay*ybcyo(ii)
   END DO
   DO ii = xs, i, incx
      ktmpy(ii,1) = ktmpy(ii,1) + wt*gax*ybcxo(ii)
   END DO
ELSE IF (j == ys) THEN
   DO ii = xs, i, incx
      IF (ii /= i) THEN
         ktmpy(ii,1) = ktmpy(ii,1) + wt*gay*ybcyo(ii)
      END IF
      IF (ii == i) THEN
         ktmpy(ii,1) = ktmpy(ii,1) + wt*gax
      END IF
   END DO
END IF

! Update values from xbc's
IF (i /= xs) THEN
   DO jj = ys, j, incy
      ktmpx(jj,1) = ktmpx(jj,1) + wt*gay*xbcyo(jj)
   END DO
   DO jj = ys, j-incy, incy
      ktmpx(jj,1) = ktmpx(jj,1) + wt*gax*xbcxo(jj)
   END DO
ELSE IF (i == xs) THEN
   DO jj = ys, j, incy
      IF (jj /= j) THEN
         ktmpx(jj,1) = ktmpx(jj,1) + wt*gax*xbcxo(jj)
      END IF
      IF (jj == j) THEN
         ktmpx(jj,1) = ktmpx(jj,1) + wt*gay
      END IF
   END DO
END IF

! Having found ktmpx, ktmpy base values, can formulate values for all l, m
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
      ktmpx(:,indx) = sh*ktmpx(:,1)
      ktmpy(:,indx) = sh*ktmpy(:,1)
   END DO
   ! Reset values
   plm2 = plm1
   plm1 = pl
   oal = al
END DO

! Finished updating ktmp for given cell
RETURN

END SUBROUTINE conkt
