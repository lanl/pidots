SUBROUTINE check(iexit)

!-------------------------------------------------------------
!
! Check the read input that includes the problem size,
!  the angular quadrature, the boundary conditions, 
!  the iteration controls, the material map, and
!  the spatial moments
!
!-------------------------------------------------------------

USE invar
USE totvar
IMPLICIT NONE
INTEGER, INTENT(INOUT) :: iexit
INTEGER :: n
REAL*8, DIMENSION(apo) :: leng

! Before starting IF constructs, set up ordinate lengths to make sure they're 1 or less
DO n = 1, apo
   leng(n) = SQRT(ang(n,1)*ang(n,1) + ang(n,2)*ang(n,2))
END DO

! Spatial moment order
IF (lambda < 0) THEN
   WRITE(8,'(/,3x,A)') "ERROR: illegal entry for lambda. Must be zero or greater."
   iexit = iexit + 1

! Method of solution
ELSE IF (meth /= 0 .AND. meth /= 1) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Method value 'meth' must be 0 or 1."
   iexit = iexit + 1

! Cell number of cells, groups, materials
ELSE IF (nxt < 1) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal number of x cells. Must be positive."
   iexit = iexit + 1
ELSE IF (nyt < 1) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal number of y cells. Must be positive."
   iexit = iexit + 1
ELSE IF (npx*npy /= isize) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal number of P_i x P_j sub-domains."
   iexit = iexit + 1
ELSE IF (ng < 1) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal number of energy groups. Must be positive."
   iexit = iexit + 1
ELSE IF (nm < 1) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal number of materials. Must be positive."
   iexit = iexit + 1
   
! Cell sizes
ELSE IF (MINVAL(dxt) <= 0.) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal x cell dimension, dx. Must be positive."
   iexit = iexit + 1
ELSE IF (MINVAL(dyt) <= 0.) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal y cell dimension, dy. Must be positive."
   iexit = iexit + 1

! Material map
ELSE IF (MINVAL(matt) < 1 .OR. MAXVAL(matt) > nm) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal value in material map. Must be in [1, #materials]."
   iexit = iexit + 1
   
! BCs
ELSE IF (q1bc /= 0 .AND. q1bc /= 1) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal Q1 BC. Must be 0-Vac or 1-Input."
   iexit = iexit + 1
ELSE IF (q2bc /= 0 .AND. q2bc /= 1) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal Q2 BC. Must be 0-Vac or 1-Input."
   iexit = iexit + 1
ELSE IF (q3bc /= 0 .AND. q3bc /= 1) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal Q3 BC. Must be 0-Vac or 1-Input."
   iexit = iexit + 1
ELSE IF (q4bc /= 0 .AND. q4bc /= 1) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal Q4 BC. Must be 0-Vac or 1-Input."
   iexit = iexit + 1

! Iteration control parameters
ELSE IF (err <= 0.) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal convergence criterion. Must be positive."
   iexit = iexit + 1
ELSE IF (bitmx <= 0) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal max PBJ iterations, bitmx. Must be positive."
   iexit = iexit + 1
ELSE IF (itmx <= 0) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal max inner iterations, itmx. Must be positive."
   iexit = iexit + 1
ELSE IF (tolr <= 0.) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal tolerance setting, tolr. Must be positive."
   iexit = iexit + 1
ELSE IF (iall < 0 .OR. iall > lambda) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal value for moments to converge, iall. Must be in [0, lambda]."
   iexit = iexit + 1

! Checks on the solution
ELSE IF (ichk < 0) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal value for solution check flag, iall."
   WRITE(8,'(/,3x,A)') "iall is 0 (skip check) or positive integer"
   iexit = iexit + 1
ELSE IF (tchk < 0.) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal value for solution check tolerance, tchk. Must be positive."
   iexit = iexit + 1
   
! Checks on the angular quadrature
ELSE IF (MINVAL(ang) <= 0. .OR. MAXVAL(ang) >= 1.) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal value for direction cosine. Must be entered positive, less than 1."
   iexit = iexit + 1
ELSE IF (MINVAL(w) <= 0. .OR. MAXVAL(w) > 0.25) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal value for weight. Must be entered positive, less than 0.25."
   iexit = iexit + 1
ELSE IF (MAXVAL(leng) >= 1.) THEN
   WRITE(8,'(/,3x,A)') "ERROR: a discrete ordinate has length >= 1 based on mu, eta values."
   iexit = iexit + 1

! Checks on the editing input
ELSE IF (momp < 0 .OR. momp > lambda) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal value for max moment to print, momp. Must be between 0 and lambda."
   iexit = iexit + 1
ELSE IF (pmoaf /= 0 .AND. pmoaf /= 1 .AND. pmoaf /= 2) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal value for flag for printing outward angular flux at boundaries."
   WRITE(8,'(/,3X,A)') "       Must be 0/1/2 = none/zero moment only/all moments"
   iexit = iexit + 1

! Or there is nothing wrong, report that and continue
ELSE
   WRITE(8,'(/,3X,A,//)') "No input errors have been detected."

END IF
   

RETURN
END SUBROUTINE check
