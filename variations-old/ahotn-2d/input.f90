SUBROUTINE input(qdfile, xsfile, srcfile, mtfile, bcfile, iexit)

!-------------------------------------------------------------
!
!    Read the input from the input file
!
!    Comments below demonstrate order of the reading
!
!    Dependency: 
!           angle   = gets the angular quadrature data
!           readxs  = reads the cross sections
!           readsrc = reads the source distribution
!           readmt  = read the material map
!           readbc  = read the boundary conditions
!           check   = input check on all the values
!
!    Allows for dynamic allocation. Uses module to hold all 
!      input variables: invar, totvar
!
!-------------------------------------------------------------

USE invar
USE totvar
IMPLICIT NONE
INTEGER :: i, j, n, ierr
INTEGER, INTENT(OUT) :: iexit
INTEGER, DIMENSION(16) :: ipak
! File Names
CHARACTER(8), INTENT(OUT) :: qdfile, xsfile, srcfile, mtfile, bcfile
LOGICAL :: ex1, ex2, ex3, ex4, ex5
REAL*8 :: wtsum
REAL*8, DIMENSION(3) :: rpak

INCLUDE 'mpif.h'

! Read the title of the case
READ(7,103) title
103 FORMAT(A80)

! Read Problem Size Specification:
!   lambda => LAMDBA, the AHOT spatial order
!   meth  => = 0/1 = AHOT-N/AHOT-N-ITM
!   qdord => Angular quadrature order
!   qdtyp => Angular quadrature type = 0/1/2 = TWOTRAN/EQN/Read-in
!   nxt    => Number of 'x' cells
!   nyt    => Number of 'y' cells
!   npx   => Number of processes in x-dimension
!   npy   => Number of processes in y-dimension
!   ng    => Number of groups
!   nm    => Number of materials

READ(7,*) lambda, meth
READ(7,*) qdord, qdtyp

iexit = 0

! Check that the order given greater than zero and is even
IF (qdord <= 0) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal value for qdord. Must be greater than zero."
   iexit = iexit + 1
ELSE IF (MOD(qdord,2) /= 0) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal value for the quadrature order. Even #s only."
   iexit = iexit + 1
END IF

! Finish reading problem size
READ(7,*) nxt, nyt
READ(7,*) npx, npy
READ(7,*) ng, nm

! Set sizes for spatial mesh, read mesh sizes
!   dxt  => cell x-dimension
!   dyt  => cell y-dimension
ALLOCATE(dxt(nxt), dyt(nyt))

READ(7,*) (dxt(i), i = 1, nxt)
READ(7,*) (dyt(j), j = 1, nyt)

! Read the boundary condition flags
!   *bc = 0/1 = Vacuum/Value from bcfile
READ(7,*) q1bc, q2bc, q3bc, q4bc

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
IF (q1bc==1 .OR. q2bc==1 .OR. q3bc==1 .OR. q4bc==1) THEN
   INQUIRE(FILE = bcfile, EXIST = ex5)
   IF (ex5 == .FALSE.) THEN
      WRITE(8,'(/,3X,A)') "ERROR: File does not exist for reading."
      iexit = iexit + 1
   END IF
END IF

! Read the iteration parameters
!   err    => Pointwise relative convergence criterion
!   bitmx  => Parallel block Jacobi maximum number of iterations
!   itmx   => Maximum number of iterations
!   iall   => Scalar Flux Spatial Moments Converged [0->LAMBDA]
!   tolr   => Tolerence for Convergence Check: determine whether abs or rel difference
READ(7,*) err, bitmx, itmx, iall, tolr

! Read the Solution Check Control:
!   ichk    => Frequency of check, 0 implies skip check
!   tchk    => Tolerence of solution check
READ(7,*) ichk, tchk

! Read the optional editing parameters
!   momp   => highest moment to be printed, [0->LAMBDA]
!   pmoaf  => flag to print moments of angular flux at boundaries: 0 for none, 1 for zero moment, 2 for all moments
READ(7,*) momp
READ(7,*) pmoaf

! Set up the extra needed info from the read input
apo = (qdord*(qdord+2))/8
order = lambda+1
ordsq = order*order

! Material Map
CALL readmt(mtfile)
 
   ! Angular quadrature
ALLOCATE(ang(apo,2), w(apo))
IF (qdtyp == 2) THEN
   INQUIRE(FILE=qdfile, EXIST=ex3)
   IF (qdfile == '        ' .OR. ex3 == .FALSE.) THEN
      WRITE(8,'(/,3x,A)') "ERROR: illegal entry for the qdfile name."
      iexit = iexit + 1
   END IF
   OPEN(UNIT=10, FILE=qdfile)
   READ(10,*)
   READ(10,*) (ang(n,1),ang(n,2),w(n),n=1,apo)
   ! Renormalize all the weights
   wtsum = SUM(w)
   DO n = 1, apo
      w(n) = w(n) * 0.25/wtsum
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

! Prepare initial information to go to all the processes in the system
! Integer package
ipak(1) = lambda
ipak(2) = meth
ipak(3) = nxt
ipak(4) = nyt
ipak(5) = npx
ipak(6) = npy
ipak(7) = ng
ipak(8) = nm
ipak(9) = bitmx
ipak(10) = itmx
ipak(11) = iall
ipak(12) = ichk
ipak(13) = apo
ipak(14) = order
ipak(15) = ordsq
ipak(16) = iexit
! Prepare the real package
rpak(1) = err
rpak(2) = tolr
rpak(3) = tchk
! Broadcast the information
CALL MPI_BCAST(ipak,16,MPI_INTEGER,irank,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(rpak,3,MPI_DOUBLE_PRECISION,irank,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(ang,(2*apo),MPI_DOUBLE_PRECISION,irank,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(w,apo,MPI_DOUBLE_PRECISION,irank,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(sigt,(ng*nm),MPI_DOUBLE_PRECISION,irank,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(sigs,(ng*ng*nm),MPI_DOUBLE_PRECISION,irank,MPI_COMM_WORLD,ierr)

RETURN
END SUBROUTINE input
