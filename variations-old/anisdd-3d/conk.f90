SUBROUTINE conk(n,i,j,k,incx,incy,incz,xs,ys,zs,mu,omega,ktmpz,ktmpy,ktmpx)

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
INTEGER, INTENT(IN) :: n, i, j, k 
INTEGER, INTENT(IN) :: incx, incy, incz, xs, ys, zs
INTEGER :: ii, jj, kk, ieq, indx, fact, lmm, lpm, l, ll
REAL*8, INTENT(IN) :: mu, omega
REAL*8 :: wt, clm, she, sho, plm2, plm1, pl, plm1m, plm, plmp1
REAL*8, DIMENSION(nmom,xys), INTENT(OUT) :: ktmpz
REAL*8, DIMENSION(nmom,xzs), INTENT(OUT) :: ktmpy
REAL*8, DIMENSION(nmom,yzs), INTENT(OUT) :: ktmpx
REAL*8, DIMENSION(0:anord) :: al, oal

! Initialize ktmps
ktmpz = 0.0
ktmpy = 0.0
ktmpx = 0.0

! Set the weight
wt = w(n)

! All the bc matrices and gamma matrices needed come from solvar

! Perform the operations to update from the different BCs
! Start with zbc's
IF (k /= zs) THEN
   DO jj = ys, j, incy
      DO ii = xs, i-incx, incx
         ieq = ii + (jj-1)*nx
         ktmpz(1,ieq) = ktmpz(1,ieq) + wt*gayz*zbcyzo(ii,jj)
      END DO
   END DO
   DO jj = ys, j-incy, incy
      DO ii = xs, i, incx
         ieq = ii + (jj-1)*nx
         ktmpz(1,ieq) = ktmpz(1,ieq) + wt*gaxz*zbcxzo(ii,jj)
      END DO
   END DO
   DO jj = ys, j, incy
      DO ii = xs, i, incx
         ieq = ii + (jj-1)*nx
         ktmpz(1,ieq) = ktmpz(1,ieq) + wt*gaxy*zbcxyo(ii,jj)
      END DO
   END DO
ELSE IF (k == zs) THEN
   DO jj = ys, j, incy
      DO ii = xs, i, incx
         ieq = ii + (jj-1)*nx
         IF (ii /= i) THEN
            ktmpz(1,ieq) = ktmpz(1,ieq) + wt*gayz*zbcyzo(ii,jj)
         END IF
         IF (jj /= j) THEN
            ktmpz(1,ieq) = ktmpz(1,ieq) + wt*gaxz*zbcxzo(ii,jj)
         END IF
         IF (ii == i .AND. jj == j) THEN
            ktmpz(1,ieq) = ktmpz(1,ieq) + wt*gaxy
         END IF
      END DO
   END DO
END IF

! Update values from ybc's
IF (j /= ys) THEN
   DO kk = zs, k, incz
      DO ii = xs, i-incx, incx
         ieq = ii + (kk-1)*nx
         ktmpy(1,ieq) = ktmpy(1,ieq) + wt*gayz*ybcyzo(ii,kk)
      END DO
   END DO
   DO kk = zs, k, incz
      DO ii = xs, i, incx
         ieq = ii + (kk-1)*nx
         ktmpy(1,ieq) = ktmpy(1,ieq) + wt*gaxz*ybcxzo(ii,kk)
      END DO
   END DO
   DO kk = zs, k-incz, incz
      DO ii = xs, i, incx
         ieq = ii + (kk-1)*nx
         ktmpy(1,ieq) = ktmpy(1,ieq) + wt*gaxy*ybcxyo(ii,kk)
      END DO
   END DO
ELSE IF (j == ys) THEN
   DO kk = zs, k, incz
      DO ii = xs, i, incx
         ieq = ii + (kk-1)*nx
         IF (ii /= i) THEN
            ktmpy(1,ieq) = ktmpy(1,ieq) + wt*gayz*ybcyzo(ii,kk)
         END IF
         IF (kk /= k) THEN
            ktmpy(1,ieq) = ktmpy(1,ieq) + wt*gaxy*ybcxyo(ii,kk)
         END IF
         IF (ii == i .AND. kk == k) THEN
            ktmpy(1,ieq) = ktmpy(1,ieq) + wt*gaxz
         END IF
      END DO
   END DO
END IF

! Update values from xbc's
IF (i /= xs) THEN
   DO kk = zs, k, incz
      DO jj = ys, j, incy
         ieq = jj + (kk-1)*ny
         ktmpx(1,ieq) = ktmpx(1,ieq) + wt*gayz*xbcyzo(jj,kk)
      END DO
   END DO
   DO kk = zs, k, incz
      DO jj = ys, j-incy, incy
         ieq = jj + (kk-1)*ny
         ktmpx(1,ieq) = ktmpx(1,ieq) + wt*gaxz*xbcxzo(jj,kk)
      END DO
   END DO
   DO kk = zs, k-incz, incz
      DO jj = ys, j, incy
         ieq = jj + (kk-1)*ny
         ktmpx(1,ieq) = ktmpx(1,ieq) + wt*gaxy*xbcxyo(jj,kk)
      END DO
   END DO
ELSE IF (i == xs) THEN
   DO kk = zs, k, incz
      DO jj = ys, j, incy
         ieq = jj + (kk-1)*ny
         IF (jj /= j) THEN
            ktmpx(1,ieq) = ktmpx(1,ieq) + wt*gaxz*xbcxzo(jj,kk)
         END IF
         IF (kk /= k) THEN
            ktmpx(1,ieq) = ktmpx(1,ieq) + wt*gaxy*xbcxyo(jj,kk)
         END IF
         IF (jj == j .AND. kk == k) THEN
            ktmpx(1,ieq) = ktmpx(1,ieq) + wt*gayz
         END IF
      END DO
   END DO
END IF

! Having found ktmpx, kmpty, ktmpz base values, can formulate values for all l, m
plm1 = 0.0
pl = 1.0
DO l = 0, anord
   al = 0.0      ! Initialize
   IF (l /= 0) CALL legpoly(l,incx,mu,plm2,plm1,pl)
   al(0) = pl
   indx = l**2 + 1
   ktmpx(indx,:) = SQRT(2.0*l+1.0)*pl*ktmpx(1,:)
   ktmpy(indx,:) = SQRT(2.0*l+1.0)*pl*ktmpy(1,:)
   ktmpz(indx,:) = SQRT(2.0*l+1.0)*pl*ktmpz(1,:)
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
      ktmpx(indx,:) = she*ktmpx(1,:)
      ktmpy(indx,:) = she*ktmpy(1,:)
      ktmpz(indx,:) = she*ktmpz(1,:)
      ktmpx(indx+1,:) = sho*ktmpx(1,:)
      ktmpy(indx+1,:) = sho*ktmpy(1,:)
      ktmpz(indx+1,:) = sho*ktmpz(1,:)
   END DO
   ! Reset values
   plm2 = plm1
   plm1 = pl
   oal = al
END DO

! Finished updating ktmp for given cell
RETURN

END SUBROUTINE conk
