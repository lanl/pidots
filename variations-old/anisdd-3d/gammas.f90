SUBROUTINE gammas(incx,mu,sig,omega)

!-------------------------------------------------------------
!
! Compute gamma elements from known analytic form of gamma
! 
!-------------------------------------------------------------

USE solvar
IMPLICIT NONE
INTEGER, INTENT(IN) :: incx
INTEGER :: bigl, l, ll, lmm, lpm, indx, fact
REAL*8, INTENT(IN) :: mu, sig, omega
REAL*8 :: alf, plm2, plm1, pl, plm1m, plm, plmp1, she, sho, clm
REAL*8, DIMENSION(0:sord-1) :: al, oal
REAL*8, DIMENSION(4,4) :: amat
REAL*8, DIMENSION(4,(nmom+3)) :: bmat, gmat

bigl = sord - 1

! Initialize matrix values
amat = 0.0
bmat = 0.0
gmat = 0.0

! Set the elements of bmat
! Elements from the balance equation
plm1 = 0.0
pl = 1.0
DO l = 0, bigl
   al = 0.0      ! Initialize
   IF (l /= 0) CALL legpoly(l,incx,mu,plm2,plm1,pl)
   al(0) = pl
   indx = l**2 + 1
   bmat(1,indx) = SQRT(2.0*l + 1.0)*pl       ! l0 moment
   DO ll = 1, l                              ! lm moments
      plm1m = oal(ll-1)
      plm   = al(ll-1)
      CALL aslegf(l,ll,incx,mu,plm1m,plm,plmp1)
      al(ll) = plmp1
      lmm = l - ll
      lpm = l + ll
      clm = (2.0*l + 1.0)*REAL(fact(lmm))/REAL(fact(lpm))
      she = 2.0*SQRT(clm)*al(ll)*COS(ll*omega)
      sho = 2.0*SQRT(clm)*al(ll)*SIN(ll*omega)
      ! Now have the odd and even values, can set in bmat
      indx = l**2 + 2*ll
      bmat(1,indx) = she
      bmat(1,indx+1) = sho
   END DO
   ! Reset values
   plm2 = plm1
   plm1 = pl
   oal = al
END DO

! Elements corresponding to the incoming
bmat(1,nmom+1) = ez
bmat(1,nmom+2) = ey
bmat(1,nmom+3) = ex

! Elements corresponding to the DD relations
bmat(2,nmom+1) = 1.0
bmat(3,nmom+2) = 1.0
bmat(4,nmom+3) = 1.0

! Set the elements of the inverted amat directly since it's only 4x4
alf = 1.0/(sig+2.0*ez+2.0*ey+2.0*ex)
amat(1,1) = alf
amat(2,1) = alf*2.0
amat(3,1) = amat(2,1)
amat(4,1) = amat(2,1)

amat(1,2) = alf*ez
amat(2,2) = alf*(-sig-2.0*ey-2.0*ex)
amat(3,2) = 2.0*amat(1,2)
amat(4,2) = amat(3,2)

amat(1,3) = alf*ey
amat(2,3) = 2.0*amat(1,3)
amat(3,3) = alf*(-sig-2.0*ez-2.0*ex)
amat(4,3) = amat(2,3)

amat(1,4) = alf*ex
amat(2,4) = 2.0*amat(1,4)
amat(3,4) = amat(2,4)
amat(4,4) = alf*(-sig-2.0*ez-2.0*ey)

! Form gmat as the product of inverted amat and bmat
gmat = MATMUL(amat,bmat)

! Store the gmat elements
! Set the gaa from gmat, first equation, first nmom equations
gaa = gmat(1,1:nmom)
! Set gaxy, gaxz, & gayz from gmat, first equation, last 3 elements
gaxy = gmat(1,nmom+1)
gaxz = gmat(1,nmom+2)
gayz = gmat(1,nmom+3)
! Set gxya, gxza, & gyza from gmat, last 3 equations, first nmom equations of each
gxya = gmat(2,1:nmom)
gxza = gmat(3,1:nmom)
gyza = gmat(4,1:nmom)
! Set the gxy**, gxz**, gyz** from gmat, bottom right 3x3 corner
gxyxy = gmat(2,nmom+1)
gxzxy = gmat(3,nmom+1)
gyzxy = gmat(4,nmom+1)
gxyxz = gmat(2,nmom+2)
gxzxz = gmat(3,nmom+2)
gyzxz = gmat(4,nmom+2)
gxyyz = gmat(2,nmom+3)
gxzyz = gmat(3,nmom+3)
gyzyz = gmat(4,nmom+3)

RETURN
END SUBROUTINE gammas
