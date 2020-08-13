PROGRAM pidosits

!-------------------------------  PIDOSITS.1.3D  -----------------------------------
!
!  Parallel Integral Discrete Ordinates & Source Iteration Transport Solver  
!  By RJZerr, December 2010
!
!  Using DD equations, solves neutron transport with
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
!             3. DD with NO flux fixup and NO acceleration
!             4. Multigroup with Isotropic Downscattering
!             5. Printing and Solution Editing Options
!             6. Multi-file Input/Output
!             7. Full domain BCs include Vac/Refl/Inp
!
!  Use makefile in same directory to compile and link.
! 
!  UPDATES: add in latest improvements to PBJ: transpose, CG initial guess,
!             pre-conditioned CG, storing more variables in solvar, MPI_WTIME,
!             index ordering, etc.
!     updates taken from PBJ code and AHOT-N formalism added back in
!
!
!-----------------------------------------------------------------------------------

USE invar
USE totvar
USE solvar
IMPLICIT NONE
CHARACTER(8) :: infile, outfile, qdfile, xsfile, srcfile, mtfile, bcfile
INTEGER :: statin, statout, ferr, ierr, iexit
INTEGER, DIMENSION(26) :: ipak
REAL*8, DIMENSION(4) :: rpak
REAL*8 :: tstart, ttosolve, tjmat, tsolve, tend, tts, tj, ts, te, tres
REAL*8 :: t1, t2, t3, t4, td1, td2, td3, td4, tdiff1, tdiff2, tdiff3, tdiff4
LOGICAL :: existence

INCLUDE 'mpif.h'

! Initialize MPI
CALL MPI_INIT(ierr)

! Get the start time
tstart = MPI_WTIME()

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

! Broadcast the fer to all processes and stop program if necessary
CALL MPI_BCAST(ferr,1,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
IF (ferr == 1) THEN
   CALL MPI_FINALIZE(ierr)
   STOP
END IF

! Set up the introductory info
IF (irank == root) CALL version

! Read input data:
! Input will call dependency algorithms, namely the input check
IF (irank == root) THEN
   CALL input(qdfile, xsfile, srcfile, mtfile, bcfile, iexit)
ELSE
   ! Other processes get the integer package first, then unpacks
   CALL MPI_BCAST(ipak,26,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
   lambda = ipak(1)
   meth   = ipak(2)
   tpose  = ipak(3)
   sym    = ipak(4)
   idos   = ipak(5)
   invf   = ipak(6)
   pcf    = ipak(7)
   pcty   = ipak(8)
   nxt    = ipak(9)
   nyt    = ipak(10)
   nzt    = ipak(11)
   npx    = ipak(12)
   npy    = ipak(13)
   npz    = ipak(14)
   ng     = ipak(15)
   nm     = ipak(16)
   bitmx  = ipak(17)
   itmx   = ipak(18)
   iall   = ipak(19)
   apo    = ipak(20)
   order  = ipak(21)
   ordsq  = ipak(22)
   ordcb  = ipak(23)
   itp    = ipak(24)
   sfp    = ipak(25)
   iexit  = ipak(26)
   ! Other processes receive the bc array
   CALL MPI_BCAST(bc,6,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
   ! Other processes then receive the real package
   CALL MPI_BCAST(rpak,4,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
   ! Allocate data before continuing
   ALLOCATE(ang(apo,3), w(apo), sigt(nm,ng), sigs(nm,ng,ng))
   CALL MPI_BCAST(ang,(3*apo),MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
   CALL MPI_BCAST(w,apo,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
   CALL MPI_BCAST(sigt,(ng*nm),MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
   CALL MPI_BCAST(sigs,(ng*ng*nm),MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
   ! Unpack rpak
   sorw = rpak(1)
   err  = rpak(2)
   cgr  = rpak(3)
   tolr = rpak(4)
END IF

! Stop all processes if there was an error detected. 0 process writes out error messages
IF (iexit > 0) THEN
   CALL MPI_FINALIZE(ierr)
   STOP
END IF

! Echo the input data:
IF (irank == root) CALL echo(infile, outfile, qdfile, xsfile, srcfile, mtfile, bcfile)

! All processes call for the parallel setup
CALL psds

IF (irank == root) THEN
   ! Now each process has their information, can deallocate
   DEALLOCATE(dxt,dyt,dzt,matt,st)
   DEALLOCATE(posinz,neginz,posiny,neginy,posinx,neginx)
END IF

! Solve the transport problem:
! Solve will call dependency algorithms: inner, weight, sweep
CALL solve(ttosolve,tjmat,tsolve)

! Print the output
IF (irank == root) CALL output

! ------------------ Final Actions, Timing Data ----------------------------------
! Get the times, relative to 'tstart'
tts = ttosolve - tstart
IF (meth == 1) tj = tjmat - tstart
ts = tsolve - tstart

! Reduce to Max values for reporting
CALL MPI_REDUCE(tts,t1,1,MPI_DOUBLE_PRECISION,MPI_MAX,root,MPI_COMM_WORLD,ierr)
IF (meth == 1) CALL MPI_REDUCE(tj,t2,1,MPI_DOUBLE_PRECISION,MPI_MAX,root,MPI_COMM_WORLD,ierr)
CALL MPI_REDUCE(ts,t3,1,MPI_DOUBLE_PRECISION,MPI_MAX,root,MPI_COMM_WORLD,ierr)

! Compute the differences and get max of them
IF (meth == 1) THEN
   td1 = tjmat - ttosolve
   td2 = tsolve - tjmat
   CALL MPI_REDUCE(td1,tdiff1,1,MPI_DOUBLE_PRECISION,MPI_MAX,root,MPI_COMM_WORLD,ierr)
   CALL MPI_REDUCE(td2,tdiff2,1,MPI_DOUBLE_PRECISION,MPI_MAX,root,MPI_COMM_WORLD,ierr)
END IF
td3 = tsolve - ttosolve
CALL MPI_REDUCE(td3,tdiff3,1,MPI_DOUBLE_PRECISION,MPI_MAX,root,MPI_COMM_WORLD,ierr)

! Time the end of the job
tend = MPI_WTIME()
te = tend - tstart
CALL MPI_REDUCE(te,t4,1,MPI_DOUBLE_PRECISION,MPI_MAX,root,MPI_COMM_WORLD,ierr)
td4 = tend - tsolve
CALL MPI_REDUCE(td4,tdiff4,1,MPI_DOUBLE_PRECISION,MPI_MAX,root,MPI_COMM_WORLD,ierr)

! Get timing resolution
tres = MPI_WTICK()

! Print the relevant times of the execution
IF (irank == root) THEN
   WRITE(8,'(/,2X,A,/)') "Fortran95 Timing Estimates with MPI_WTIME..."
   WRITE(8,100)
   WRITE(8,101) "InputWrk", t1, t1
   IF (meth == 1) WRITE(8,101) "IDOmatrx", t2, tdiff1
   IF (meth == 1) WRITE(8,101) "IDOitrns", t3, tdiff2
   WRITE(8,101) "SolveTot", t3, tdiff3
   WRITE(8,101) "PrintOut", t4, tdiff4
   WRITE(8,'(/,5X,A,ES9.2)') "Timing Resolution (s): ", tres
   WRITE(8,102)
   100 FORMAT(5X,'WorkDone',3X,'Absolute(s)',6X,'Difference(s)')
   101 FORMAT(5X,A8,3X,F9.3,5X,F9.3)
   102 FORMAT(//,'*********************   END PROGRAM   ************************')
END IF

! Deallocate the allocated arrays
IF (irank == root) DEALLOCATE(ssum,nxvec,nyvec,nzvec)
IF (irank == root .AND. sfp == 1) DEALLOCATE(flux)
DEALLOCATE(dx,dy,dz,mat,ang,w,sigt,sigs,s,bcnvf)
DEALLOCATE(psii,psio)
IF (meth == 0) DEALLOCATE(f)
IF (meth == 1) DEALLOCATE(phi)

! End the PIDOSITS program
CALL MPI_FINALIZE(ierr)

STOP
END PROGRAM pidosits
