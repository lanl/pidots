SUBROUTINE input(qdfile, xsfile, srcfile, mtfile, bcfile, iexit)

!-------------------------------------------------------------
!
!    Read the input from the input file
!
!    Comments below demonstrate order of the reading
!
!    Dependency: 
!           angle   = gets the angular quadrature data
!           readmt  = reads the material map from file
!           readxs  = reads the cross sections
!           readsrc = reads the source distribution
!           readbc  = read the boundary conditions from file
!           check   = input check on all the values
!
!    Allows for dynamic allocation. Uses module to hold all 
!      input variables: invar, totvar
!
!-------------------------------------------------------------

USE invar
IMPLICIT NONE
INTEGER :: i, j, k, n, ierr
INTEGER, INTENT(OUT) :: iexit
! File Names
CHARACTER(8), INTENT(OUT) :: qdfile, xsfile, srcfile, mtfile, bcfile
LOGICAL :: ex1, ex2, ex3, ex4, ex5
REAL*8 :: wtsum

! Read the title of the case
READ(7,103) title
103 FORMAT(A80)

! Read Problem Size Specification:
!   qdord => Angular quadrature order
!   qdtyp => Angular quadrature type = 0/1/2 = TWOTRAN/EQN/Read-in
!   nx   => Number of 'x' cells
!   ny   => Number of 'y' cells
!   nz   => Number of 'z' cells
!   ng    => Number of groups
!   nm    => Number of materials

READ(7,*) qdord, qdtyp

iexit = 0
warn = 0
! Check that the order given greater than zero and is even
IF (qdord <= 0) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal value for qdord. Must be greater than zero."
   iexit = iexit + 1
ELSE IF (MOD(qdord,2) /= 0) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal value for the quadrature order. Even #s only."
   iexit = iexit + 1
END IF

! Finish reading problem size
READ(7,*) nx, ny, nz
READ(7,*) ng, nm

! Set sizes for spatial mesh, read mesh sizes
!   dx  => cell x-dimension
!   dy  => cell y-dimension
!   dz  => cell z-dimension
ALLOCATE(dx(nx), dy(ny), dz(nz))

READ(7,*) (dx(i), i = 1, nx)
READ(7,*) (dy(j), j = 1, ny)
READ(7,*) (dz(k), k = 1, nz)

! Read the boundary conditions
!   *bc = 0/1/2 = Vacuum/Reflective/Input
READ(7,*) bc(1), bc(2)       ! Positive, Negative z-direction
READ(7,*) bc(3), bc(4)       ! Positive, Negative y-direction
READ(7,*) bc(5), bc(6)       ! Positive, Negative x-direction

! Read the names of files with cross sections and source distribution
READ(7,104) mtfile
READ(7,104) qdfile
READ(7,104) xsfile
READ(7,104) srcfile
READ(7,104) bcfile
104 FORMAT(A8)

! Perform quick checks on the files
INQUIRE(FILE = mtfile, EXIST = ex4)
INQUIRE(FILE = xsfile, EXIST = ex1)
INQUIRE(FILE = srcfile, EXIST = ex2)
IF (ex1 == .FALSE. .OR. ex2 == .FALSE. .OR. ex4 == .FALSE.) THEN
   WRITE(8,'(/,3X,A)') "ERROR: File does not exist for reading."
   iexit = iexit + 1
END IF
IF (MAXVAL(bc) == 2) THEN
   INQUIRE(FILE=bcfile, EXIST = ex5)
   IF (ex5 == .FALSE.) THEN
      WRITE(8,'(/,3X,A)') "ERROR: File does not exist for reading."
      iexit = iexit + 1
   END IF
END IF

! Read the iteration parameters
!   err    => Pointwise relative convergence criterion
!   itmx   => Maximum number of iterations
!   tolr   => Tolerence for Convergence Check: determine whether abs or rel difference
READ(7,*) err, itmx, tolr

! Set up the extra needed info from the read input
apo = (qdord*(qdord+2))/8

! Material map
CALL readmt(mtfile,iexit)
 
! Angular quadrature
ALLOCATE(ang(apo,3), w(apo))
IF (qdtyp == 2) THEN
   INQUIRE(FILE=qdfile, EXIST=ex3)
   IF (qdfile == '        ' .OR. ex3 == .FALSE.) THEN
      WRITE(8,'(/,3X,A)') "ERROR: illegal entry for the qdfile name."
      iexit = iexit + 1
   END IF
   OPEN(UNIT=10, FILE=qdfile)
   READ(10,*)
   READ(10,*) (ang(n,1),ang(n,2),ang(n,3),w(n),n=1,apo)
   ! Renormalize all the weights
   wtsum = SUM(w)
   DO n = 1, apo
      w(n) = w(n) * 0.125/wtsum
   END DO
ELSE
   CALL angle(iexit)
END IF

IF (qdtyp == 2) CLOSE(UNIT=10)

! Call for the input check
CALL check(iexit)
! Call to read the cross sections and source; do their own input check
CALL readxs(xsfile,iexit)
CALL readsrc(srcfile,iexit)

RETURN
END SUBROUTINE input
