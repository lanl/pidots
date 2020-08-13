SUBROUTINE gammas(incx,mu,sig,omega)

!-------------------------------------------------------------
!
! Compute gamma elements from known analytic form of gamma
! 
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER, INTENT(IN) :: incx
INTEGER :: l, ll, lmm, lpm, indx, fact
REAL*8, INTENT(IN) :: mu, sig, omega
REAL*8 :: alf, plm2, plm1, pl, plm1m, plm, plmp1, sh, p, clm
REAL*8, DIMENSION(0:anord) :: al, oal
REAL*8, DIMENSION(3,3) :: amat
REAL*8, DIMENSION(3,(nmom+2)) :: bmat, gmat

! Initialize the amat and bmat values
amat = 0.0
bmat = 0.0
gmat = 0.0

! Set the elements of bmat
! Elements from the balance equation
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
      ! Now have the value, can set in bmat
      indx = l*(l+1)/2 + ll + 1
      p = 1.0
      IF (ll > 0) p = 2.0
      bmat(1,indx) = p*sh
   END DO
   ! Reset values
   plm2 = plm1
   plm1 = pl
   oal = al
END DO

! Elements corresponding to incoming flux
bmat(1,nmom+1) = ey
bmat(1,nmom+2) = ex

! Elements corresponding to DD relations
bmat(2,nmom+1) = 1.0
bmat(3,nmom+2) = 1.0

! Set the elements of the inverted amat directly since it's only 3x3
alf = 1.0/(sig+2.0*ey+2.0*ex)
amat(1,1) = alf
amat(2,1) = alf*2.0
amat(3,1) = amat(2,1)
amat(1,2) = alf*ey
amat(2,2) = alf*(-sig-2.0*ex)
amat(3,2) = alf*2.0*ey
amat(1,3) = alf*ex
amat(2,3) = alf*2.0*ex
amat(3,3) = alf*(-sig-2.0*ey)

! Form gmat as the product of inverted amat on bmat
gmat = MATMUL(amat,bmat)

! Store the gmat elements
! Set the gaa from gmat, first equation, first nmom equations
gaa = gmat(1,1:nmom)
! Set gax & gay from gmat, first equation, last 2 elements
gax = gmat(1,nmom+1)
gay = gmat(1,nmom+2)
! Set gxa and gya from gmat, second & third equations, first nmom equations
gxa = gmat(2,1:nmom)
gya = gmat(3,1:nmom)
! Set gxx, gxy, gyx, gyy from gmat, bottom right 2x2 corner
gxx = gmat(2,nmom+1)
gxy = gmat(2,nmom+2)
gyx = gmat(3,nmom+1)
gyy = gmat(3,nmom+2)

RETURN
END SUBROUTINE gammas
