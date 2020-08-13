SUBROUTINE idomats(v,spx,spy,spz)

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
INTEGER, INTENT(IN) :: v, spx, spy, spz
INTEGER :: ieq, i, j, k, info, t
REAL*8, DIMENSION(neq) :: src

! Initialize the matrices
jmat(:,:,v)     = 0.0
kmat(:,:,:,v)   = 0.0
jpsi(:,:,:,v)   = 0.0
kpsi(:,:,:,:,v) = 0.0

!-------------------------------------------------------------------
! Can now start the single sweep for constructing all matrices
CALL matsweep(v,spx,spy,spz)

IF (v == 1) THEN
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
      sv(j,1) = DOT_PRODUCT(jmat(:,j,1),src)
   END DO

   ! Multiply the source vector by Jpsi to form the other saved vector
   DO j = 1, 8
      DO i = 1, bcs
         av(i,j,1) = DOT_PRODUCT(jpsi(:,i,j,1),src)
      END DO
   END DO
END IF

! Form the LHS matrix
jmat(:,:,v) = -jmat(:,:,v)
DO t = 1, neq
   jmat(t,t,v) = 1.0 + jmat(t,t,v)
END DO

!------------------------------------------------------------------------
! Always factor the jmat matrix
CALL dgetrf(neq,neq,jmat(:,:,v),neq,piv(:,v),info)
IF (info /= 0) THEN
   WRITE (8,'(//,1X,A)') "ERROR: jmat either has illegal value or is singular, cannot factor."
   STOP
END IF

! Have needed values sv, av, jmat, kmat, jpsi, kpsi, piv

RETURN
END SUBROUTINE idomats
