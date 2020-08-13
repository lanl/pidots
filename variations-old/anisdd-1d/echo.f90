SUBROUTINE echo(infile, outfile, qdfile, xsfile, srcfile, mtfile, bcfile)

!-------------------------------------------------------------
!
!    Echo the input
!
!-------------------------------------------------------------

USE invar
USE totvar
IMPLICIT NONE
INTEGER :: i, g, m, n
! File Names
CHARACTER(8), INTENT(IN) :: infile, outfile, qdfile, xsfile, srcfile, mtfile, bcfile


! Start the echo
WRITE (8,*) "-------------------- BEGIN INPUT ECHO ------------------------------------"

! Write the title of the case and the basics of the problem setup
WRITE (8,'(//,1X, A)') title
IF (meth == 0) WRITE (8,106) "AHOT-N-SI-3D Method"
IF (meth == 1) WRITE (8,106) "AHOT-N-ITMM-3D Method"
IF (meth == 1 .AND. tpose == 1) WRITE(8,106) "ITMM operators constructed in transpose form."
IF (meth == 1 .AND. (sym==1 .OR. idos==1)) WRITE (8,106) "IDO matrix to be made symmetric"
WRITE (8,105) "Angular Order N = ", qdord
WRITE (8,105) "Max Scattering Expansion Order L = ", anord
WRITE (8,105) "Number of x-cells = ", nxt
WRITE (8,'(1X,A,I3,A1,I3,A1,I3)') "Process Grid ", npx
WRITE (8,105) "Number of energy groups = ", ng
WRITE (8,105) "Number of materials = ", nm
105 FORMAT(1X,A,I5)
106 FORMAT(1X,A)

! Write the iteration controls
WRITE (8,'(/,1X,A)') "Iteration Control Parameters"
WRITE (8,105) "Maximum number of parallel block iterations (angular flux) = ", bitmx
WRITE (8,105) "Maximum number of iterations = ", itmx
WRITE (8,107) "Pointwise convergence criterion = ", err
IF (idos == 1) WRITE (8,107) "CG residual convergence criterion = ", cgr
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
WRITE (8,'(1X,A,I,/)') "Scalar Flux Print (SFP) Flag 0/1/2=No/Yes/Root P Only: ", sfp

! Write the angular quadrature information        
IF (qdtyp == 0) THEN
   WRITE (8,'(1X,A)') "Legendre Quadrature for Slab-type Discrete Ordinates (positive x)"
ELSE
   WRITE (8,'(1X,A)') "Read-In Discrete Ordinates (positive x)"
END IF
WRITE (8, '(2X, A1, T8, A2, T16, A1)') "n", "mu", "w"
DO n = 1, apo
   WRITE (8, '(1X, I2, T6, F7.5, T15, F7.5)') n, ang(n), w(n)
END DO

! Write the boundary conditions
WRITE (8,*)
WRITE (8,*) "Boundary Conditions: 0-Vacuum, 1-Reflective, 2-Input (See file)"
WRITE (8, '(1X,A2,2X,A2)') "X+", "X-"
WRITE (8, '(T2,I1,3X,I1,/)') bc(1), bc(2)

! Write the computational cell data: dx
WRITE (8,*)
WRITE (8,105) "x-dimension of cells 1 to ", nxt
WRITE (8,108) dxt
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
WRITE (8,'(/,1X,A,I6)') "Material map cells 1 to ", nxt
WRITE (8,109) (matt(i), i = 1, nxt)
109 FORMAT(1X,40I4)

! Write the source information
DO g = 1, ng
   WRITE (8,*)
   WRITE (8,110) "(Isotropic) Source for Group: ", g, " Cells 1 to ", nxt
   WRITE (8,108) (st(i,g), i = 1, nxt)
END DO
110 FORMAT(1X,A,I3,A,I6)

! End the input echo
WRITE (8,*)
WRITE (8,*) "------------------------- END INPUT ECHO ---------------------------------"
WRITE (8,*)

RETURN
END SUBROUTINE echo
