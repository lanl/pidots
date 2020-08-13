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

! Method of solution
IF (meth /= 0 .AND. meth /= 1) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Method value 'meth' must be 0 or 1."
   iexit = iexit + 1
ELSE IF (tpose /= 0 .AND. tpose /= 1) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Transpose flag 'tpose' must be 0 or 1."
   iexit = iexit + 1
ELSE IF (sym /= 0 .AND. sym /= 1) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Symmetry flag 'sym' must be 0 or 1."
   iexit = iexit + 1
ELSE IF (idos /= 0 .AND. idos /= 1) THEN
   WRITE(8,'(/,3X,A)') "ERROR: IDO Solver value 'idos' must be 0 or 1."
   iexit = iexit + 1

! Anisotropic order
ELSE IF (anord < 0) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Anisotropic Order, 'anord' (L), must be non-negative."
   iexit = iexit + 1

! Cell number of cells, groups, materials
ELSE IF (nxt < 1) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal number of x cells. Must be positive."
   iexit = iexit + 1
ELSE IF (npx /= isize) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal number of P_i sub-domains."
   iexit = iexit + 1
ELSE IF (npx > nxt) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Number processes cannot be greater than number cells."
   iexit = iexit + 1
ELSE IF (ng < 1) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal number of energy groups. Must be positive."
   iexit = iexit + 1
ELSE IF (nm < 1) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal number of materials. Must be positive."
   iexit = iexit + 1
   
! Cell sizes
ELSE IF (MINVAL(dxt) <= 0.) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal x cell dimension, dx. Must be positive."
   iexit = iexit + 1

! Material map
ELSE IF (MINVAL(matt) < 1 .OR. MAXVAL(matt) > nm) THEN
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
ELSE IF (idos == 1 .AND. cgr <= 0.) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal CG residual conv. crit. for IDOS=1. Must be positive."
   iexit = iexit + 1
ELSE IF (bitmx <= 0) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal max PBJ iterations, bitmx. Must be positive."
   iexit = iexit + 1
ELSE IF (itmx <= 0) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal max inner iterations, itmx. Must be positive."
   iexit = iexit + 1
ELSE IF (tolr <= 0.) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal tolerance setting, tolr. Must be positive."
   iexit = iexit + 1

! Checks on the solution
ELSE IF (ichk < 0) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal value for solution check flag, iall."
   WRITE(8,'(/,3X,A)') "iall is 0 (skip check) or positive integer"
   iexit = iexit + 1
ELSE IF (tchk < 0.) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal value for solution check tolerance, tchk. Must be positive."
   iexit = iexit + 1
   
! Checks on the angular quadrature
ELSE IF (MINVAL(ang) <= 0. .OR. MAXVAL(ang) >= 1.) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal value for direction cosine. Must be entered positive, less than 1."
   WRITE (8,*) MINVAL(ang), MINLOC(ang), MAXVAL(ang), MAXLOC(ang)
   iexit = iexit + 1
ELSE IF (MINVAL(w) <= 0. .OR. MAXVAL(w) > 0.5) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal value for weight. Must be entered positive, less than 0.5."
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
