SUBROUTINE idotf(f,fsv,fpiv,fpo,fav,fjmt,fjps)

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
IMPLICIT NONE
INTEGER :: i, j, info
INTEGER, DIMENSION(neq), INTENT(IN) :: fpiv
REAL*8, DIMENSION(neq), INTENT(OUT) :: f
REAL*8, DIMENSION(neq), INTENT(IN) :: fsv
REAL*8, DIMENSION(bcs,8), INTENT(OUT) :: fpo
REAL*8, DIMENSION(bcs,8), INTENT(IN) :: fav
REAL*8, DIMENSION(neq,neq), INTENT(IN) :: fjmt
REAL*8, DIMENSION(neq,bcs,8), INTENT(IN) :: fjps

! Compute the solution directly
f = fsv

! Use LAPACK solvers if matrix was factored only
CALL dgetrs('T',neq,1,fjmt,neq,fpiv,f,neq,info)
IF (info /= 0) THEN
   WRITE (8,'(//,1X,A)') "ERROR: dgetrs cannot solve."
   STOP
END IF

! Compute the new angular flux values at the boundaries from jpsi and kpsi
! Solve for psio
DO j = 1, 8
   DO i = 1, bcs
      fpo = DOT_PRODUCT(fjps(:,i,j),f)
   END DO
END DO

! Add in RHS vector
fpo = fpo + fav

RETURN
END SUBROUTINE idotf
