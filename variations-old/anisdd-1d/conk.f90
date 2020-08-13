SUBROUTINE conk(n,i,xs,incx,mu,ktmp)

!-----------------------------------------------------------
!
! Constructs the ktmp row for a given octant and cell and 
!  returns the values to afcm.
!
! Takes in from afcm: All the computed values not in solvar
!     n,i,xs,incx,mu
! Returns to afcm a temporary matrix ktmp that fills in the
!  particular ktmp#
!
!------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER, INTENT(IN) :: n, i, xs, incx
INTEGER :: l
REAL*8, INTENT(IN) :: mu
REAL*8 :: wt, plm1, pl
REAL*8, DIMENSION(sord), INTENT(OUT) :: ktmp

! Initialize ktmp
ktmp = 0.0

! Set the weight
wt = w(n)

! All the bc matrices and gamma matrices needed come from solvar

! Perform the operations to update from the different BCs
IF (i /= xs) THEN
   plm1 = 0.0
   pl   = 1.0
   DO l = 0, anord
      ktmp(l+1) = ktmp(l+1) + wt*pl*gax*xbco
      CALL legpoly(l,incx,mu,plm1,pl)
   END DO
ELSE IF (i == xs) THEN
   plm1 = 0.0
   pl   = 1.0
   DO l = 0, anord
      ktmp(l+1) = ktmp(l+1) + wt*pl*gax
      CALL legpoly(l,incx,mu,plm1,pl)
   END DO
END IF

! Finished updating ktmp for given cell
RETURN

END SUBROUTINE conk
