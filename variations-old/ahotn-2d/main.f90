PROGRAM pidosits

!-------------------------------  PIDOSITS.1.2D------------------------------
!
!  Parallel Integral Discrete Ordinates & Source Iteration Transport Solver  
!  By RJZerr, March 2009
!
!  Using AHOT-N/WDD equations, solves neutron transport with
!   1. Integral Discrete Ordinates Approach from Azmy
!   2. Source Iteration Approach
!
!  Solves in parallel with the Parallel Block Jacobi approach. Each 
!   sub-domain solves for scalar flux and outgoing angular flux. Angular fluxes
!   are passed to neighbors to solve for new scalar flux. Iterate until
!   convergence is reached.
!
!  Features:  1. Double Precision
!             2. Module files for input, solution, and timing variables
!             3. General Order Nodal
!             4. Multigroup with Isotropic Downscattering
!             5. Printing and Solution Editing Options
!             6. Multi-file Input/Output
!
!  Use makefile in same directory to compile and link.
!
!----------------------------------------------------------------------------

USE invar
USE totvar
USE solvar
USE timevar
IMPLICIT NONE
CHARACTER(8) :: infile, outfile, qdfile, xsfile, srcfile, mtfile, bcfile
INTEGER :: statin, statout, ferr, ierr, root, iexit
INTEGER, DIMENSION(16) :: ipak
REAL*8, DIMENSION(3) :: rpak
REAL :: t1, t2, t3, t4, td1, td2, td3, tdiff1, tdiff2, tdiff3
LOGICAL :: existence

INCLUDE 'mpif.h'

! Initialize MPI
CALL MPI_INIT(ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD,irank,ierr)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD,isize,ierr)

! Set the process that does the reading
root = 0

! Have the root process above take care of input and output
IF (irank == root) THEN
   CALL GETARG(1, infile, statin)
   CALL GETARG(2, outfile, statout)
   IF (statin == -1 .OR. statout == -1) THEN
      PRINT *, "Failed reading input and/or output file names."
      ferr = 1
   ELSE
      OPEN (UNIT = 7, FILE = infile, STATUS = "OLD", ACTION = "READ")
      ! Check if the output file exists or not, then open appropriately
      INQUIRE (FILE = outfile, EXIST = existence)
      IF (existence == .TRUE.) THEN
         OPEN (UNIT = 8, FILE = outfile, STATUS = "OLD", ACTION = "WRITE")
      ELSE
         OPEN (UNIT = 8, FILE = outfile, STATUS = "NEW", ACTION = "WRITE")
      END IF
      ferr = 0
   END IF
END IF

! Broadcast the ferr to all process and stop program if necessary
CALL MPI_BCAST(ferr,1,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
IF (ferr == 1) THEN
   CALL MPI_FINALIZE(ierr)
   STOP
END IF

! Set up the introductory info
IF (irank == root) THEN
   CALL version
END IF

! Initialize arrays with data
ipak = 0
rpak = 0.0

! Read input data:
! Input will call dependency algorithms, namely the input check
IF (irank == root) THEN
   CALL input(qdfile, xsfile, srcfile, mtfile, bcfile, iexit)
ELSE
   ! Other processes get the integer package first, then unpacks
   CALL MPI_BCAST(ipak,16,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
   lambda = ipak(1)
   meth   = ipak(2)
   nxt    = ipak(3)
   nyt    = ipak(4)
   npx    = ipak(5)
   npy    = ipak(6)
   ng     = ipak(7)
   nm     = ipak(8)
   bitmx  = ipak(9)
   itmx   = ipak(10)
   iall   = ipak(11)
   ichk   = ipak(12)
   apo    = ipak(13)
   order  = ipak(14)
   ordsq  = ipak(15)
   iexit  = ipak(16)
   ! Other processes then receive the real package
   CALL MPI_BCAST(rpak,3,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
   ! Allocate data before continuing with rest of broadcasts
   ALLOCATE(ang(apo,2), w(apo), sigt(nm,ng), sigs(nm,ng,ng))
   CALL MPI_BCAST(ang,(2*apo),MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
   CALL MPI_BCAST(w,apo,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
   CALL MPI_BCAST(sigt,(ng*nm),MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
   CALL MPI_BCAST(sigs,(ng*ng*nm),MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
   ! Unpack rpak
   err  = rpak(1)
   tolr = rpak(2)
   tchk = rpak(3)
END IF

! Stop all processes if there was an error detected. 0 process writes out error messages
IF (iexit > 0) THEN
   CALL MPI_FINALIZE(ierr)
   STOP
END IF

! Echo the input data: by process 0 only
IF (irank == root) THEN
   CALL echo(infile, outfile, qdfile, xsfile, srcfile, mtfile, bcfile)
END IF

! All processes call for the parallel setup
CALL psds(root)

IF (irank == root) THEN
   ! Each process now has their information, can deallocate
   DEALLOCATE(dxt,dyt,matt,st,psi1t,psi2t,psi3t,psi4t)
END IF

! Solve the transport problem:
! Solve will call dependency algorithms: inner, weight, sweep
CALL solve(root)

! Print the output
IF (irank == root) THEN
   CALL output
END IF

! Get the times
CALL MPI_REDUCE(ttosolve,t1,1,MPI_REAL,MPI_MAX,root,MPI_COMM_WORLD,ierr)
IF (meth == 1) CALL MPI_REDUCE(tjmat,t2,1,MPI_REAL,MPI_MAX,root,MPI_COMM_WORLD,ierr)
CALL MPI_REDUCE(tsolve,t3,1,MPI_REAL,MPI_MAX,root,MPI_COMM_WORLD,ierr)
! Compute the differences and get max of them
IF (meth == 0) THEN
   td2 = tsolve - ttosolve
ELSE
   td1 = tjmat - ttosolve
   td2 = tsolve - tjmat
END IF
IF (meth == 1) CALL MPI_REDUCE(td1,tdiff1,1,MPI_REAL,MPI_MAX,root,MPI_COMM_WORLD,ierr)
CALL MPI_REDUCE(td2,tdiff2,1,MPI_REAL,MPI_MAX,root,MPI_COMM_WORLD,ierr)

! Time the end of the job
CALL CPU_TIME(tend)
CALL MPI_REDUCE(tend,t4,1,MPI_REAL,MPI_MAX,root,MPI_COMM_WORLD,ierr)
td3 = tend - tsolve
CALL MPI_REDUCE(td3,tdiff3,1,MPI_REAL,MPI_MAX,root,MPI_COMM_WORLD,ierr)

! Print the relevant times of the execution
IF (irank == root) THEN
   WRITE(8,'(/,2X,A,/)') "Fortran95 Timing Estimates with CPU_TIME..."
   WRITE(8,100)
   WRITE(8,101) "InputWrk", t1, t1
   IF (meth == 1) WRITE(8,101) "IDOmatrx", t2, tdiff1
   WRITE(8,101) "SolveTot", t3, tdiff2
   WRITE(8,101) "PrintOut", t4, tdiff3
   WRITE(8,102)
   100 FORMAT(5X,'WorkDone',3X,'Absolute(s)',6X,'Difference(s)')
   101 FORMAT(5X,A8,3X,F9.3,5X,F9.3)
   102 FORMAT(//,'*********************   END PROGRAM  ************************')
END IF

! Deallocate the allocated arrays
IF (irank == root) DEALLOCATE(ssum,flux,nxvec,nyvec)
DEALLOCATE(dx,dy,mat,ang,w,sigt,sigs,s,bcnvf)
DEALLOCATE(psi1,psi2,psi3,psi4,psio)

! End the PIDOSITS program
CALL MPI_FINALIZE(ierr)

STOP
END PROGRAM pidosits
