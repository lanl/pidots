SUBROUTINE legpoly(l,sgm,mu,plm1,pl)

! ---------------------------------------------
!
! Uses a recursive relationship to compute the
!  next Legendre polynomial
!
! ---------------------------------------------

IMPLICIT NONE

INTEGER, INTENT(IN) :: l, sgm
REAL*8, INTENT(IN) :: mu
REAL*8, INTENT(INOUT) :: plm1, pl
REAL*8 :: plp1

plp1 = (1.0/(l+1.0))*((2.0*l+1.0)*sgm*mu*pl-l*plm1)
plm1 = pl
pl   = plp1

END SUBROUTINE legpoly
