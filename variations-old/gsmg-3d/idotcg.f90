SUBROUTINE idotcg(fpiv,f,fsv,fjmt)

!-------------------------------------------------------------
!
!  Integral discrete ordinates transport
!
!  Add sv to kmat*psii to get the new phi
!  Compute psio with jpsi*(phi+src) + kpsi*psii
!
!-------------------------------------------------------------

USE invar
IMPLICIT NONE
INTEGER :: info
INTEGER, DIMENSION(neq), INTENT(IN) :: fpiv
REAL*8, DIMENSION(neq), INTENT(OUT) :: f
REAL*8, DIMENSION(neq), INTENT(IN) :: fsv
REAL*8, DIMENSION(neq,neq), INTENT(IN) :: fjmt

! Compute the solution directly
! Always use sv and product of kmat and psii to form RHS
f = fsv

! Use LAPACK solvers if matrix was factored only
CALL dgetrs('T',neq,1,fjmt,neq,fpiv,f,neq,info)
IF (info /= 0) THEN
   WRITE (8,'(//,1X,A)') "ERROR: dgetrs cannot solve."
   STOP
END IF

! Coarsest grid: don't need angular flux values.

RETURN
END SUBROUTINE idotcg
