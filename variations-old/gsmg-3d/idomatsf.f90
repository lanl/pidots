SUBROUTINE idomatsf(spx,spy,spz,fjmt,fkmt,fjps,fkps,fpiv)

!-------------------------------------------------------------
!
!  Integral discrete ordinates matrices
!
!   Control the construction of jmat, kmat, jpsi, kpsi
!   Matrices nned to be made only once, then idot.f90 will
!    use them to solve the system of equations
!
!-------------------------------------------------------------

USE invar
IMPLICIT NONE
INTEGER, INTENT(IN) :: spx, spy, spz
INTEGER :: ieq, i, j, k, info, t, rsx, rsy, rsz
INTEGER, DIMENSION(neq), INTENT(OUT) :: fpiv
REAL*8, DIMENSION(neq,neq), INTENT(OUT) :: fjmt
REAL*8, DIMENSION(bcs,neq,8), INTENT(OUT) :: fkmt
REAL*8, DIMENSION(neq,bcs,8), INTENT(OUT) :: fjps
REAL*8, DIMENSION(bcs2,bcs2,apo,8), INTENT(OUT) :: fkps

! Initialize the matrices
fjmt = 0.0
fkmt = 0.0
fjps = 0.0
fkps = 0.0

!-------------------------------------------------------------------
! Can now start the single sweep for constructing all matrices
rsx = 0
rsy = 0
rsz = 0
CALL matsweep(rsx,rsy,rsz,spx,spy,spz,fjmt,fkmt,fjps,fkps)

! Form the LHS matrix
fjmt = -fjmt
DO t = 1, neq
   fjmt(t,t) = 1.0 + fjmt(t,t)
END DO

!------------------------------------------------------------------------
! Always factor the jmat matrix
CALL dgetrf(neq,neq,fjmt,neq,fpiv,info)
IF (info /= 0) THEN
   WRITE (8,'(//,1X,A)') "ERROR: jmat either has illegal value or is singular, cannot factor."
   STOP
END IF

! Have needed values fjmt, fkmt, fjps, fkps, fpiv

RETURN
END SUBROUTINE idomatsf
