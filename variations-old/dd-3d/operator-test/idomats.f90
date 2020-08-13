SUBROUTINE idomats

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
USE solvar
IMPLICIT NONE
INTEGER :: ieq, i, j, k, info, t

! Initialize the matrices
kmatz = 0.0
kmaty = 0.0
kmatx = 0.0

jpsiz = 0.0
jpsiy = 0.0
jpsix = 0.0

kpsizz = 0.0
kpsizy = 0.0
kpsizx = 0.0
kpsiyz = 0.0
kpsiyy = 0.0
kpsiyx = 0.0
kpsixz = 0.0
kpsixy = 0.0
kpsixx = 0.0


! Allocate the iteration Jacobian
ALLOCATE(jmat(neq,neq))
jmat = 0.0

!-------------------------------------------------------------------
! Can now start the single sweep for constructing all matrices
CALL matsweep

! Begin constructing the RHS
DO k = 1, nz
   DO j = 1, ny
      DO i = 1, nx
         ieq = i + (j-1)*nx + (k-1)*xys
         src(ieq) = s(i,j,k)/sigs(mat(i,j,k))
      END DO
   END DO
END DO

! Multiply the source vector by the Jacobi iteration matrix
DO j = 1, neq
   sv(j) = DOT_PRODUCT(jmat(:,j),src)
END DO

! Form the LHS matrix
jmat = -jmat
DO t = 1, neq
   jmat(t,t) = 1.0 + jmat(t,t)
END DO

!------------------------------------------------------------------------
! Need to factor jmat
   
! Allocate the piv vector as it's now needed
ALLOCATE(piv(neq))

! Always factor the jmat matrix
CALL dgetrf(neq,neq,jmat,neq,piv,info)
IF (info /= 0) THEN
   WRITE (8,'(//,1X,A)') "ERROR: jmat either has illegal value or is singular, cannot factor."
   STOP
END IF

! Have needed values src, sv, jmat, kmat, jpsi, kpsi, piv

RETURN
END SUBROUTINE idomats
