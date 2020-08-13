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
   leng(n) = SQRT(ang(n,1)*ang(n,1) + ang(n,2)*ang(n,2) + ang(n,3)*ang(n,3))
END DO

! Cell number of cells, groups, materials
IF (nxt < 1) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal number of x cells. Must be positive."
   iexit = iexit + 1
ELSE IF (nyt < 1) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal number of y cells. Must be positive."
   iexit = iexit + 1
ELSE IF (nzt < 1) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal number of z cells. Must be positive."
   iexit = iexit + 1
ELSE IF (npx*npy*npz /= isize) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal number of P_i x P_j x P_k sub-domains."
   iexit = iexit + 1
ELSE IF (npx > nxt .OR. npy > nyt .OR. npz > nzt) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Number processes in a dimension cannot be greater than number cells."
   iexit = iexit + 1
ELSE IF (nm < 1) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal number of materials. Must be positive."
   iexit = iexit + 1
ELSE IF (2**(ns-1) /= npx) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal number of V-Cycle stages."
   iexit = iexit + 1
  
! Cell sizes
ELSE IF (MINVAL(dxt) <= 0.) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal x cell dimension, dx. Must be positive."
   iexit = iexit + 1
ELSE IF (MINVAL(dyt) <= 0.) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal y cell dimension, dy. Must be positive."
   iexit = iexit + 1
ELSE IF (MINVAL(dzt) <= 0.) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal z cell dimension, dz. Must be positive."
   iexit = iexit + 1

! Material map
ELSE IF (MINVAL(matt) < 1 .OR. MAXVAL(matt) > nm) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal value in material map. Must be in [1, #materials]."
   iexit = iexit + 1
   
! BCs
ELSE IF (MINVAL(bc) < 0 .OR. MAXVAL(bc) > 1) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal BC. Must be 0-Vac, 1-Reflective."
   iexit = iexit + 1

! Iteration control parameters
ELSE IF (err <= 0.) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal convergence criterion. Must be positive."
   iexit = iexit + 1
ELSE IF (vitmx <= 0) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal max V-Cycles, vitmx. Must be positive."
   iexit = iexit + 1
ELSE IF (bitmx1 <= 0) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal max pre-correction PBJ iterations, bitmx1. Must be positive."
   iexit = iexit + 1
ELSE IF (bitmx2 <= 0) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal max post-correction PBJ iterations, bitmx2. Must be positive."
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

! Checks on the editing input
ELSE IF (sfp < 0 .OR. sfp > 2) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal value for scalar flux print flag, sfp. Must be 0/1/2=No/Yes/Root P only."
   iexit = iexit + 1
ELSE IF (itp /= 0 .AND. itp /= 1) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal value for iteration print flag, itp. Must be 0/1=No/Yes."
   iexit = iexit + 1

! Or there is nothing wrong, report that and continue
ELSE
   WRITE(8,'(/,3X,A,//)') "No input errors have been detected."

END IF
   

RETURN
END SUBROUTINE check
