SUBROUTINE echo(infile, outfile, qdfile, xsfile, srcfile, mtfile, bcfile)

!-------------------------------------------------------------
!
!    Echo the input
!
!-------------------------------------------------------------

USE invar
IMPLICIT NONE
INTEGER :: i, j, k, g, m, n
! File Names
CHARACTER(8), INTENT(IN) :: infile, outfile, qdfile, xsfile, srcfile, mtfile, bcfile

! Start the echo
WRITE (8,*) "-------------------- BEGIN INPUT ECHO ------------------------------------"

! Write the title of the case and the basics of the problem setup
WRITE (8,'(//,1X, A)') title

WRITE (8,105) "Angular Order N = ", qdord
WRITE (8,105) "Number of x-cells = ", nx
WRITE (8,105) "Number of y-cells = ", ny
WRITE (8,105) "Number of z-cells = ", nz
WRITE (8,105) "Number of energy groups = ", ng
WRITE (8,105) "Number of materials = ", nm
105 FORMAT(1X,A,I5)
106 FORMAT(1X,A)

! Write the iteration controls
WRITE (8,'(/,1X,A)') "Iteration Control Parameters"
WRITE (8,105) "Maximum number of iterations = ", itmx
WRITE (8,107) "Pointwise convergence criterion = ", err
107 FORMAT(1X,A,ES10.3)

! Write the names of the files used
WRITE (8,'(/,1X,A,A8)') "Data read from input file: ", infile
WRITE (8,'(/,1X,A,A8)') "Material map read from file: ", mtfile
IF (qdtyp == 2) WRITE (8,'(A,A8)') "Quadrature data from file: ", qdfile
WRITE (8,'(1X,A,A8)') "Cross sections from file: ", xsfile
WRITE (8,'(1X,A,A8)') "Source data from file: ", srcfile
IF (MAXVAL(bc) == 1) THEN
   WRITE(8,'(1X,A,A8)') "BC data from file: ", bcfile
END IF
WRITE (8,'(1X,A,A8,/)') "Output written to file: ", outfile

! Write the angular quadrature information        
IF (qdtyp == 0) THEN
   WRITE (8,'(1X,A)') "TWOTRAN-type Discrete Ordinates/Octant"
ELSE IF (qdtyp == 1) THEN
   WRITE (8,'(1X,A)') "EQN-type Discrete Ordinates/Octant"
ELSE
   WRITE (8,'(1X,A)') "Read-In Discrete Ordinates/Octant"
END IF
WRITE (8, '(2X, A1, T8, A2, T16, A3, T26, A2, T35, A1)') "n", "mu", "eta", "xi", "w"
DO n = 1, apo
   WRITE (8, '(1X, I2, T6, F7.5, T15, F7.5, T24, F7.5, T33, F7.5)') n, ang(n, 1), ang(n, 2), ang(n,3), w(n)
END DO

! Write the boundary conditions
WRITE (8,*)
WRITE (8,*) "Boundary Conditions: 0-Vacuum, 1-Reflective, 2-Input (See file)"
WRITE (8, '(1X,A2,2X,A2,2X,A2,2X,A2,2X,A2,2X,A2)') "Z+", "Z-", "Y+", "Y-", "X+", "X-"
WRITE (8, '(T2,I1,3X,I1,3X,I1,3X,I1,3X,I1,3X,I1,/)') bc(1), bc(2), bc(3), bc(4), bc(5), bc(6)

! Write the computational cell data: dx, dy, dz
WRITE (8,*)
WRITE (8,105) "x-dimension of cells 1 to ", nx
WRITE (8,108) dx
WRITE (8,*)
WRITE (8,105) "y-dimension of cells 1 to ", ny
WRITE (8,108) dy
WRITE (8,*)
WRITE (8,105) "z-dimension of cells 1 to ", nz
WRITE (8,108) dz
108 FORMAT(2X, 8ES10.3)

! Write the cross section data
WRITE (8,'(/,1X,A)') "Cross Section Data"
WRITE (8,'(2X, A8, T15, A5, T22, A7, T35, A7)') "Material", "Group", "Sigma_t", "Sigma_s"
DO m = 1, nm
   DO g = 1, ng
      WRITE (8,'(T2, I5, T16, I3, T22, ES10.3, T35, ES10.3)') m, g, sigt(m,g), ssum(m,g)
   END DO
END DO

! Write the material map
! Loop over the z-planes
DO k = 1, nz
   WRITE (8,'(/,1X,A,I5)') "Material map for z-plane: ", k
   WRITE (8,'(1X,A)', ADVANCE = "NO") "Row "
   DO i = 1, nx
      WRITE (8,'(I3,1X)',ADVANCE = "NO") i
   END DO
   WRITE (8,*)
   DO j = ny, 1, -1
      WRITE (8,109) j, (mat(i,j,k), i = 1, nx)
   END DO
END DO
109 FORMAT(I3,1X,128I4)

! Write the source information
DO g = 1, ng
   DO k = 1, nz
      DO j = 1, ny
         WRITE (8,*)
         WRITE (8,110) "Source for Group: ", g, " z-plane(k): ", k, " Row(j): ", j
         WRITE (8,108) (s(i,j,k,g), i = 1, nx)
      END DO
   END DO
END DO
110 FORMAT(1X,A,I3,A,I3,A,I3)

! End the input echo
WRITE (8,*)
WRITE (8,*) "------------------------- END INPUT ECHO ---------------------------------"
WRITE (8,*)

RETURN
END SUBROUTINE echo
