SUBROUTINE legpoly(l,sgm,mu,plm2,plm1,pl)

! ---------------------------------------------
!
! Uses a recursive relationship to compute the
!  next Legendre polynomial
!
! ---------------------------------------------

IMPLICIT NONE

INTEGER, INTENT(IN) :: l, sgm
REAL*8, INTENT(IN) :: mu, plm2, plm1
REAL*8, INTENT(OUT) :: pl

pl = (1.0/REAL(l))*((2.0*l-1.0)*sgm*mu*plm1-(l-1.0)*plm2)

END SUBROUTINE legpoly
