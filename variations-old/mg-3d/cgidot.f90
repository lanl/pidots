SUBROUTINE cgidot(v)

!-------------------------------------------------------------
!
!  Integral discrete ordinates transport
!
!  Add sv to kmat*psii to get the new phi
!  Compute psio with jpsi*(phi+src) + kpsi*psii
!
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER, INTENT(IN) :: v
INTEGER :: info

! Compute the solution directly
! Always use sv and product of kmat and psii to form RHS
phi(:,v) = sv(:,v)

! Use LAPACK solvers if matrix was factored only
CALL dgetrs('T',neq,1,jmat(:,:,v),neq,piv(:,v),phi(:,v),neq,info)
IF (info /= 0) THEN
   WRITE (8,'(//,1X,A)') "ERROR: dgetrs cannot solve."
   STOP
END IF

! Coarsest grid: don't need angular flux values.

RETURN
END SUBROUTINE cgidot
