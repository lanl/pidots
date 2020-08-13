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
ELSE IF (MOD(nxt,npx) /= 0 .OR. MOD(nyt,npy) /= 0 .OR. MOD(nzt,npz) /= 0) THEN
   WRITE(8,'(/,3X,A)') "ERROR: The number of cells in any direction must be evenly divided by the &
                                                            number of processes in that direction."
   iexit = iexit + 1
ELSE IF (sdnx < 1 .OR. sdny < 1 .OR. sdnz < 1) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Number of cells in a sub-domain must be positive."
   iexit = iexit + 1
ELSE IF (MOD(nx,sdnx) /= 0 .OR. MOD(ny,sdny) /= 0 .OR. MOD(nz,sdnz) /= 0) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Number of cells in a sub-domain must divide evenly into the number &
                                                                       of cells on each processor."
   iexit = iexit + 1
ELSE IF (sdnx /= sdny .OR. sdnx /= sdnz .OR. sdny /= sdnz) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Restricting sub-domain sizes to be equal in each direction. Fix sdnx/y/z."
   iexit = iexit + 1
ELSE IF (nx /= ny .OR. nx /= nz .OR. ny /= nz) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Restricting processor # nodes to be equal in each direction. Fix nx/y/z."
   iexit = iexit + 1
ELSE IF (nxsd < 2 .OR. nysd < 2 .OR. nzsd < 2) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Must have at least two sub-domains per P each direction."
   iexit = iexit + 1
ELSE IF (MOD(nxsd,2) /= 0 .OR. MOD(nysd,2) /= 0 .OR. MOD(nzsd,2) /= 0) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Must have even number of sub-domains per P each direction."
   iexit = iexit + 1
ELSE IF (nm < 1) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal number of materials. Must be positive."
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
ELSE IF (MINVAL(bc) /= 0 .OR. MAXVAL(bc) /= 0) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal BC. Must be 0-Vac."
   iexit = iexit + 1

! Iteration control parameters
ELSE IF (err <= 0.) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal convergence criterion. Must be positive."
   iexit = iexit + 1
ELSE IF (vitmx <= 0) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal max v-cycle iterations, vitmx. Must be positive."
   iexit = iexit + 1
ELSE IF (bitmx1 <= 0) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal max PGS pre-correction iterations, bitmx1. Must be positive."
   iexit = iexit + 1
ELSE IF (bitmx2 <= 0) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal max PGS post-correction iterations, bitmx2. Must be positive."
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
