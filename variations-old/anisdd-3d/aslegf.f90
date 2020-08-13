SUBROUTINE aslegf(l,m,sgm,mu,plm1m,plm,plmp1)

! ---------------------------------------------
!
! Uses a recursive relationship to compute the
!  next associated Legendre function
!
! ---------------------------------------------

IMPLICIT NONE

INTEGER, INTENT(IN) :: l, m, sgm
REAL*8, INTENT(IN) :: mu, plm1m, plm
REAL*8, INTENT(OUT) :: plmp1

plmp1 = (1.0/SQRT(1.0-mu**2))*((l-m+1)*sgm*mu*plm - (l+m-1)*plm1m)

END SUBROUTINE aslegf
