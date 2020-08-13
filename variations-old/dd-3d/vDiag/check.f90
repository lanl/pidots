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
IMPLICIT NONE
INTEGER, INTENT(INOUT) :: iexit
INTEGER :: n
REAL*8, DIMENSION(apo) :: leng

! Before starting IF constructs, set up ordinate lengths to make sure they're 1 or less
DO n = 1, apo
   leng(n) = SQRT(ang(n,1)*ang(n,1) + ang(n,2)*ang(n,2) + ang(n,3)*ang(n,3))
END DO

! Cell number of cells, groups, materials
IF (nx < 1) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal number of x cells. Must be positive."
   iexit = iexit + 1
ELSE IF (ny < 1) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal number of y cells. Must be positive."
   iexit = iexit + 1
ELSE IF (nz < 1) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal number of z cells. Must be positive."
   iexit = iexit + 1
ELSE IF (ng < 1) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal number of energy groups. Must be positive."
   iexit = iexit + 1
ELSE IF (nm < 1) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal number of materials. Must be positive."
   iexit = iexit + 1
   
! Cell sizes
ELSE IF (MINVAL(dx) <= 0.) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal x cell dimension, dx. Must be positive."
   iexit = iexit + 1
ELSE IF (MINVAL(dy) <= 0.) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal y cell dimension, dy. Must be positive."
   iexit = iexit + 1
ELSE IF (MINVAL(dz) <= 0.) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal z cell dimension, dz. Must be positive."
   iexit = iexit + 1

! Material map
ELSE IF (MINVAL(mat) < 1 .OR. MAXVAL(mat) > nm) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal value in material map. Must be in [1, #materials]."
   iexit = iexit + 1
   
! BCs
ELSE IF (MINVAL(bc) < 0 .OR. MAXVAL(bc) > 2) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal BC. Must be 0-Vac, 1-Reflective, or 2-Input."
   iexit = iexit + 1

! Iteration control parameters
ELSE IF (err <= 0.) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal convergence criterion. Must be positive."
   iexit = iexit + 1
ELSE IF (itmx <= 0) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal max inner iterations, itmx. Must be positive."
   iexit = iexit + 1
ELSE IF (tolr <= 0.) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal tolerance setting, tolr. Must be positive."
   iexit = iexit + 1

! Checks on the angular quadrature
ELSE IF (MINVAL(ang) <= 0. .OR. MAXVAL(ang) >= 1.) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal value for direction cosine. Must be entered positive, less than 1."
   WRITE (8,*) MINVAL(ang), MINLOC(ang), MAXVAL(ang), MAXLOC(ang)
   iexit = iexit + 1
ELSE IF (MINVAL(w) <= 0. .OR. MAXVAL(w) > 0.125) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal value for weight. Must be entered positive, less than 0.125."
   iexit = iexit + 1
ELSE IF (MINVAL(leng) <= 0.99 .OR. MAXVAL(leng) >= 1.01) THEN
   WRITE(8,'(/,3X,A)') "ERROR: a discrete ordinate has length not in range 0.99<1.00<1.01 based on mu, eta, xi values."
   iexit = iexit + 1

! Or there is nothing wrong, report that and continue
ELSE
   WRITE(8,'(/,3X,A,//)') "No input errors have been detected."

END IF
   

RETURN
END SUBROUTINE check
