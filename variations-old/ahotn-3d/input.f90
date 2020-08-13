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
INTEGER :: i, j, k, n, ierr
INTEGER, INTENT(OUT) :: iexit
INTEGER, DIMENSION(26) :: ipak
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
!   lambda => LAMBDA, the AHOT spatial order
!   meth  => = 0/1 = DD/DD-ITM
!   tpose => = 0/1 = Transpose Construction of Operators No/Yes
!   sym   => = 0/1 = No/Yes Symmetrize (I-J) matrix
!   idos  => = 0/1/2 = Direct/CG/SOR
!   invf  => = 0/1 = Factor/Invert (I-J) matrix
!   pcf   => = 0/1 = No/Yes Pre-conditioner for the CG method
!   pcty  => = 0/1/2 = Diagonal/Symmetric Gauss-Seidel/SSOR
!   sorw  => = (0,2) = SSOR weight - real double precision
!   qdord => Angular quadrature order
!   qdtyp => Angular quadrature type = 0/1/2 = TWOTRAN/EQN/Read-in
!   nxt   => Number of 'x' cells
!   nyt   => Number of 'y' cells
!   nzt   => Number of 'z' cells
!   npx   => Number of processes in x-dimension
!   npy   => Number of processes in y-dimension
!   npz   => Number of processes in z-dimension
!   ng    => Number of groups
!   nm    => Number of materials

READ(7,*) lambda, meth, tpose, sym, idos
READ(7,*) invf, pcf, pcty, sorw
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
READ(7,*) nxt, nyt, nzt
READ(7,*) npx, npy, npz
READ(7,*) ng, nm

! Set sizes for spatial mesh, read mesh sizes
!   dxt  => cell x-dimension
!   dyt  => cell y-dimension
!   dzt  => cell z-dimension
ALLOCATE(dxt(nxt), dyt(nyt), dzt(nzt))

READ(7,*) (dxt(i), i = 1, nxt)
READ(7,*) (dyt(j), j = 1, nyt)
READ(7,*) (dzt(k), k = 1, nzt)

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
!   cgr    => CG method residual convergence criterion (Used if IDOS=1)
!   bitmx  => Parallel block Jacobi maximum number of iterations
!   itmx   => Maximum number of iterations
!   iall   => Scalar Flux Spatial Moments Converged [0->LAMBDA]
!   tolr   => Tolerence for Convergence Check: determine whether abs or rel difference
READ(7,*) err, cgr, bitmx, itmx, iall, tolr

! Read the optional editing parameters
!   itp   => iteration data printing  0=no, 1=yes
!   sfp   => scalar flux print flag   0=no, 1=yes, 2=root process only
READ(7,*) itp
READ(7,*) sfp

! Set up the extra needed info from the read input
apo = (qdord*(qdord+2))/8
order = lambda+1
ordsq = order**2
ordcb = order**3

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

! Call to read the boundary condition input file; does own input check
CALL readbc(bcfile)

! Prepare the initial information to go to all the processes in the system
! Integer package
ipak(1) = lambda
ipak(2) = meth
ipak(3) = tpose
ipak(4) = sym
ipak(5) = idos
ipak(6) = invf
ipak(7) = pcf
ipak(8) = pcty
ipak(9) = nxt
ipak(10) = nyt
ipak(11) = nzt
ipak(12) = npx
ipak(13) = npy
ipak(14) = npz
ipak(15) = ng
ipak(16) = nm
ipak(17) = bitmx
ipak(18) = itmx
ipak(19) = iall
ipak(20) = apo
ipak(21) = order
ipak(22) = ordsq
ipak(23) = ordcb
ipak(24) = itp
ipak(25) = sfp
ipak(26) = iexit
! Prepare the real package
rpak(1) = sorw
rpak(2) = err
rpak(3) = cgr
rpak(4) = tolr
! Broadcast the information
CALL MPI_BCAST(ipak,26,MPI_INTEGER,irank,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(bc,6,MPI_INTEGER,irank,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(rpak,4,MPI_DOUBLE_PRECISION,irank,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(ang,(3*apo),MPI_DOUBLE_PRECISION,irank,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(w,apo,MPI_DOUBLE_PRECISION,irank,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(sigt,(ng*nm),MPI_DOUBLE_PRECISION,irank,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(sigs,(ng*ng*nm),MPI_DOUBLE_PRECISION,irank,MPI_COMM_WORLD,ierr)

RETURN
END SUBROUTINE input
