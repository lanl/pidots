SUBROUTINE idotv(v,bit)

!-------------------------------------------------------------
!
!  Integral discrete ordinates transport
!
!  For v>1 stages, only performing 1 iteration
!  psii = 0.0 always
!
!  Add sv to kmat*psii to get the new phi
!  Compute psio with jpsi*(phi+src) + kpsi*psii
!
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER, INTENT(IN) :: v, bit
INTEGER :: i, j, info

! Compute the solution directly
phi(:,v) = sv(:,v)

! Use LAPACK solvers if matrix was factored only
CALL dgetrs('T',neq,1,jmat(:,:,v),neq,piv(:,v),phi(:,v),neq,info)
IF (info /= 0) THEN
   WRITE (8,'(//,1X,A)') "ERROR: dgetrs cannot solve."
   STOP
END IF

! Compute the new angular flux values at the boundaries from jpsi and kpsi
! Solve for psio
DO j = 1, 8
   DO i = 1, bcs
      psio(i,j,v) = DOT_PRODUCT(jpsi(:,i,j,v),phi(:,v))
   END DO
END DO

! Add in RHS vector
psio(:,:,v) = psio(:,:,v) + av(:,:,v)

RETURN
END SUBROUTINE idotv
