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
USE totvar
IMPLICIT NONE
INTEGER :: i, n, ierr
INTEGER, INTENT(OUT) :: iexit
INTEGER, DIMENSION(16) :: ipak
! File Names
CHARACTER(8), INTENT(OUT) :: qdfile, xsfile, srcfile, mtfile, bcfile
LOGICAL :: ex1, ex2, ex3, ex4, ex5
REAL*8 :: wtsum
REAL*8, DIMENSION(4) :: rpak

INCLUDE 'mpif.h'

! Read the title of the case
READ(7,103) title
103 FORMAT(A80)

! Read Problem Size Specification:
!   meth  => = 0/1 = DD/DD-ITM
!   tpose => = 0/1 = Transpose Construction of Operators No/Yes
!   sym   => = 0/1 = No/Yes Symmetrize (I-J) matrix
!   idos  => = 0/1 = Direct/CG
!   qdord => Angular quadrature order
!   qdtyp => Angular quadrature type = 0/1 = LQn/Read-in
!   anord => Order of Anisotropy for Scattering Expansion, L
!   nxt   => Number of 'x' cells
!   npx   => Number of processes in x-dimension
!   ng    => Number of groups
!   nm    => Number of materials

READ(7,*) meth, tpose, sym, idos
READ(7,*) qdord, qdtyp, anord

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
READ(7,*) nxt
READ(7,*) npx
READ(7,*) ng, nm

! Set sizes for spatial mesh, read mesh sizes
!   dxt  => cell x-dimension
ALLOCATE(dxt(nxt))

READ(7,*) (dxt(i), i = 1, nxt)

! Read the boundary conditions
!   *bc = 0/1/2 = Vacuum/Reflective/Input
READ(7,*) bc(1), bc(2)       ! Positive, Negative x-direction

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
!   cgr    => CG method residual convergence criterion (Used if IDOS=1)
!   bitmx  => Parallel block Jacobi maximum number of iterations
!   itmx   => Maximum number of iterations
!   tolr   => Tolerence for Convergence Check: determine whether abs or rel difference
READ(7,*) err, cgr, bitmx, itmx, tolr

! Read the Solution Check Control:
!   ichk    => Frequency of check, 0 implies skip check
!   tchk    => Tolerence of solution check
READ(7,*) ichk, tchk

! Read the optional editing parameters
!   itp   => iteration data printing  0=no, 1=yes
!   sfp   => scalar flux print flag   0=no, 1=yes, 2=root process only
READ(7,*) itp
READ(7,*) sfp

! Set up the extra needed info from the read input
apo = qdord/2

! Material map
CALL readmt(mtfile,iexit)
 
! Angular quadrature
ALLOCATE(ang(apo), w(apo))
IF (qdtyp == 1) THEN
   INQUIRE(FILE=qdfile, EXIST=ex3)
   IF (qdfile == '        ' .OR. ex3 == .FALSE.) THEN
      WRITE(8,'(/,3X,A)') "ERROR: illegal entry for the qdfile name."
      iexit = iexit + 1
   END IF
   OPEN(UNIT=10, FILE=qdfile)
   READ(10,*)
   READ(10,*) (ang(n), w(n), n = 1, apo)
   ! Renormalize all the weights
   wtsum = SUM(w)
   DO n = 1, apo
      w(n) = w(n) * 0.125/wtsum
   END DO
ELSE
   CALL angle(iexit)
END IF

IF (qdtyp == 1) CLOSE(UNIT=10)

! Call for the input check
CALL check(iexit)
! Call to read the cross sections and source; do their own input check
CALL readxs(xsfile,iexit)
CALL readsrc(srcfile,iexit)

! Call to read the boundary condition input file; does own input check
CALL readbc(bcfile)

! Prepare the initial information to go to all the processes in the system
! Integer package
ipak(1) = meth
ipak(2) = tpose
ipak(3) = sym
ipak(4) = idos
ipak(5) = anord
ipak(6) = nxt
ipak(7) = npx
ipak(8) = ng
ipak(9) = nm
ipak(10) = bitmx
ipak(11) = itmx
ipak(12) = ichk
ipak(13) = apo
ipak(14) = iexit
ipak(15) = itp
ipak(16) = sfp
! Prepare the real package
rpak(1) = err
rpak(2) = cgr
rpak(3) = tolr
rpak(4) = tchk
! Broadcast the information
CALL MPI_BCAST(ipak,16,MPI_INTEGER,irank,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(bc,2,MPI_INTEGER,irank,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(rpak,4,MPI_DOUBLE_PRECISION,irank,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(ang,apo,MPI_DOUBLE_PRECISION,irank,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(w,apo,MPI_DOUBLE_PRECISION,irank,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(sigt,(ng*nm),MPI_DOUBLE_PRECISION,irank,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(sigs,(ng*ng*(anord+1)*nm),MPI_DOUBLE_PRECISION,irank,MPI_COMM_WORLD,ierr)

RETURN
END SUBROUTINE input
