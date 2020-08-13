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
INTEGER, DIMENSION(25) :: ipak
! File Names
CHARACTER(8), INTENT(OUT) :: qdfile, xsfile, srcfile, mtfile, bcfile
LOGICAL :: ex1, ex2, ex3, ex4, ex5
REAL*8 :: wtsum
REAL*8, DIMENSION(2) :: rpak

INCLUDE 'mpif.h'

! Read the title of the case
READ(7,103) title
103 FORMAT(A80)

! Read Problem Size Specification:
!   qdord => Angular quadrature order
!   qdtyp => Angular quadrature type = 0/1/2 = TWOTRAN/EQN/Read-in
!   nxt   => Number of 'x' cells
!   nyt   => Number of 'y' cells
!   nzt   => Number of 'z' cells
!   npx   => Number of processes in x-dimension
!   npy   => Number of processes in y-dimension
!   npz   => Number of processes in z-dimension
!   sdnx  => #i cells per sub-domain
!   sdny  => #j cells per sub-domain
!   sdnz  => #k cells per sub-domain
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
READ(7,*) nxt, nyt, nzt
READ(7,*) npx, npy, npz
READ(7,*) sdnx, sdny, sdnz
READ(7,*) nm

! Set sizes for spatial mesh, read mesh sizes
!   dxt  => cell x-dimension
!   dyt  => cell y-dimension
!   dzt  => cell z-dimension
ALLOCATE(dxt(nxt), dyt(nyt), dzt(nzt))

READ(7,*) (dxt(i), i = 1, nxt)
READ(7,*) (dyt(j), j = 1, nyt)
READ(7,*) (dzt(k), k = 1, nzt)

! Read the boundary conditions
!   *bc = 0/1 = Vacuum/Reflective
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
IF (ex1 .EQV. .FALSE. .OR. ex2 .EQV. .FALSE. .OR. ex4 .EQV. .FALSE.) THEN
   WRITE(8,'(/,3X,A)') "ERROR: File does not exist for reading."
   iexit = iexit + 1
END IF
IF (MAXVAL(bc) == 2) THEN
   INQUIRE(FILE=bcfile, EXIST = ex5)
   IF (ex5 .EQV. .FALSE.) THEN
      WRITE(8,'(/,3X,A)') "ERROR: File does not exist for reading."
      iexit = iexit + 1
   END IF
END IF

! Read the iteration parameters
!   err    => Pointwise relative convergence criterion
!   bitmx => Parallel GS maximum number of iterations, pre-correction
!   tolr   => Tolerence for Convergence Check: determine whether abs or rel difference
READ(7,*) err, bitmx, tolr

! Read the optional editing parameters
!   itp   => iteration data printing  0=no, 1=yes
!   sfp   => scalar flux print flag   0=no, 1=yes, 2=root process only
READ(7,*) itp
READ(7,*) sfp

! Set up the extra needed info from the read input
! Check ensures size restrictions are met
apo = (qdord*(qdord+2))/8
nx = nxt/npx
ny = nyt/npy
nz = nzt/npz
nxsd = nx/sdnx
nysd = ny/sdny
nzsd = nz/sdnz
nsdp = nxsd*nysd*nzsd
neq = sdnx*sdny*sdnz
xys = sdnx*sdny
xzs = sdnx*sdnz
yzs = sdny*sdnz
bcs = apo*(xys+xzs+yzs)
bcs2 = xys+xzs+yzs

! Material map
CALL readmt(mtfile,iexit)

! Angular quadrature
ALLOCATE(ang(apo,3), w(apo))
IF (qdtyp == 2) THEN
   INQUIRE(FILE=qdfile, EXIST=ex3)
   IF (qdfile == '        ' .OR. ex3 .EQV. .FALSE.) THEN
      WRITE(8,'(/,3X,A)') "ERROR: illegal entry for the qdfile name."
      iexit = iexit + 1
   END IF
   OPEN(UNIT=10, FILE=qdfile)
   READ(10,*)
   READ(10,*) (ang(n,1),ang(n,2),ang(n,3),w(n),n=1,apo)
   ! Renormalize all the weights
   wtsum = SUM(w)
   DO n = 1, apo
      w(n) = w(n) * 0.125d0/wtsum
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
!CALL readbc(bcfile)

! Prepare the initial information to go to all the processes in the system
! Integer package
ipak(1) = nx
ipak(2) = ny
ipak(3) = nz
ipak(4) = npx
ipak(5) = npy
ipak(6) = npz
ipak(7) = sdnx
ipak(8) = sdny
ipak(9) = sdnz
ipak(10) = nxsd
ipak(11) = nysd
ipak(12) = nzsd
ipak(13) = nsdp
ipak(14) = nm
ipak(15) = bitmx
ipak(16) = apo
ipak(17) = neq
ipak(18) = xys
ipak(19) = xzs
ipak(20) = yzs
ipak(21) = bcs
ipak(22) = bcs2
ipak(23) = itp
ipak(24) = sfp
ipak(25) = iexit
! Prepare the real package
rpak(1) = err
rpak(2) = tolr
! Broadcast the information
CALL MPI_BCAST(ipak,25,MPI_INTEGER,irank,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(bc,6,MPI_INTEGER,irank,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(rpak,2,MPI_DOUBLE_PRECISION,irank,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(ang,(3*apo),MPI_DOUBLE_PRECISION,irank,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(w,apo,MPI_DOUBLE_PRECISION,irank,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(sigt,nm,MPI_DOUBLE_PRECISION,irank,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(sigs,nm,MPI_DOUBLE_PRECISION,irank,MPI_COMM_WORLD,ierr)

RETURN
END SUBROUTINE input
