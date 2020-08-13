SUBROUTINE cidomats

!-------------------------------------------------------------
!
!   Coarse Integral discrete ordinates matrices
!
!   Control the construction of cjmat, ckmat, cjpsi, ckpsi
!
!-------------------------------------------------------------

USE solvar
IMPLICIT NONE
INTEGER ::  info, t

! Initialize the matrices
ckmat = 0.0
cjpsi = 0.0
ckpsi = 0.0

! Allocate the iteration Jacobian
ALLOCATE(cjmat(cneq,cneq))
cjmat = 0.0

!-------------------------------------------------------------------
! Can now start the single sweep for constructing all matrices
CALL cmatsweep

! Form the LHS matrix
cjmat = -cjmat
DO t = 1, cneq
   cjmat(t,t) = 1.0 + cjmat(t,t)
END DO

!------------------------------------------------------------------------
! Force the coarse problem to be asymmetric and solved directly
  
! Allocate the piv vector as it's now needed
ALLOCATE(cpiv(cneq))

! Always factor the cjmat matrix
CALL dgetrf(cneq,cneq,cjmat,cneq,cpiv,info)
IF (info /= 0) THEN
   WRITE (8,'(//,1X,A)') "ERROR: cjmat either has illegal value or is singular, cannot factor."
   STOP
END IF

! Have needed values cjmat, ckmat, cjpsi, ckpsi, cpiv

RETURN
END SUBROUTINE cidomats
