SUBROUTINE gammas(incx,mu,sig,amat,bmat,gmat)

!-------------------------------------------------------------
!
!  Construct the A and B matrices that make the gamma matrix
!  A and B matrices specific for an ordinate (n) and cell (x,y)
!
!  Invert A, multiply that with B and set in Gamma (gmat)
!
!  Store gmat into gamma elements
! 
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER, INTENT(IN) :: incx
INTEGER :: l, indx
REAL*8 :: pl, plm1, alf
REAL*8, INTENT(IN) :: mu, sig
REAL*8, DIMENSION(2,2) :: amat
REAL*8, DIMENSION(2,(sord+1)) :: bmat, gmat

! Initialize the amat and bmat values
amat = 0.0
bmat = 0.0
gmat = 0.0

! Set elements of bmat
! Elements from balance equation
plm1 = 0.0
pl   = 1.0
DO l = 0, anord
   indx = l + 1
   bmat(1,indx) = (2.0*l+1.0)*pl
   CALL legpoly(l,incx,mu,plm1,pl)
END DO
bmat(1,(sord+1)) = ex
! Elements from balance equation (only on psi_in)
bmat(2,(sord+1)) = 1.0

! Set the elements of the inverted amat directly since it is only 2x2
alf = 1.0/(-sig - 2.0*ex)
amat(1,1) = alf*(-1.0)
amat(1,2) = alf*(-ex)
amat(2,1) = alf*(-2.0)
amat(2,2) = alf*(sig)

! Form gmat as the product of inverted amat on bmat
gmat = MATMUL(amat,bmat)

! Store the gmat elements
! Set gaa from gmat, first equation, first L+1 elements
gaa = gmat(1,1:sord)
! Set gax from gmat, first equation, last element
gax = gmat(1,(sord+1))
! Set gxa from gmat, second equation, first L+1 elements
gxa = gmat(2,1:sord)
! Set gxa from gmat, first equation, last element
gxx = gmat(2,(sord+1))

RETURN
END SUBROUTINE gammas
