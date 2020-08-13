PROGRAM pidosits

!-------------------------------  PIDOSITS.1.1D  -----------------------------------
!
!  Parallel Integral Discrete Ordinates & Source Iteration Transport Solver  
!  By RJZerr, November 2009 (Special version to use DD: AHOT-N0 with weights=0)
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
!  Update from previous version:
!   Reduce the 3D version to 1D
!   ** Add anisotropic scattering approximations to both SI and ITMM
!   computations
!
!-----------------------------------------------------------------------------------

USE invar
USE totvar
USE solvar
USE timevar
IMPLICIT NONE
CHARACTER(8) :: infile, outfile, qdfile, xsfile, srcfile, mtfile, bcfile
INTEGER :: statin, statout, ferr, ierr, iexit, sta
INTEGER, DIMENSION(16) :: ipak
REAL*8, DIMENSION(4) :: rpak
REAL :: t1, t2, t3, t4, td1, td2, td3, td4, tdiff1, tdiff2, tdiff3, tdiff4
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
   CALL MPI_BCAST(ipak,16,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
   meth   = ipak(1)
   tpose  = ipak(2)
   sym    = ipak(3)
   idos   = ipak(4)
   anord  = ipak(5)
   nxt    = ipak(6)
   npx    = ipak(7)
   ng     = ipak(8)
   nm     = ipak(9)
   bitmx  = ipak(10)
   itmx   = ipak(11)
   ichk   = ipak(12)
   apo    = ipak(13)
   iexit  = ipak(14)
   itp    = ipak(15)
   sfp    = ipak(16)
   ! Other processes receive the bc array
   CALL MPI_BCAST(bc,2,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
   ! Other processes then receive the real package
   CALL MPI_BCAST(rpak,4,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
   ! Allocate data before continuing
   ALLOCATE(ang(apo), w(apo)) !sigt(nm,ng), sigs(nm,0:anord,ng,ng))
   ALLOCATE(sigt(nm,ng),STAT=sta)
!if (sta /= 0) print *, irank, sta
!print *, irank, sta
!   sta = sta + 0
   CALL MPI_BCAST(ang,(3*apo),MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
   CALL MPI_BCAST(w,apo,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
!   CALL MPI_BCAST(sigt,(ng*nm),MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
!   CALL MPI_BCAST(sigs,(ng*ng*(anord+1)*nm),MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
   ! Unpack rpak
   err = rpak(1)
   cgr = rpak(2)
   tolr = rpak(3)
   tchk = rpak(4)
END IF

! Stop all processes if there was an error detected. 0 process writes out error messages
IF (iexit > 0) THEN
   CALL MPI_FINALIZE(ierr)
   STOP
END IF

! Echo the input data:
IF (irank == root) CALL echo(infile, outfile, qdfile, xsfile, srcfile, mtfile, bcfile)

! All processes call for the parallel setup
!CALL psds

IF (irank == root) THEN
   ! Now each process has their information, can deallocate
   DEALLOCATE(dxt,matt,st)
   DEALLOCATE(posinx,neginx)
END IF

! Deallocate the allocated arrays
IF (irank == root) DEALLOCATE(ssum,sigs,sigt)
!IF (irank == root .AND. sfp == 1) DEALLOCATE(flux)
!DEALLOCATE(dx) !,mat,ang,w,sigt,sigs,s,bcnvf)
!DEALLOCATE(mat)
DEALLOCATE(ang)
DEALLOCATE(w)
!DEALLOCATE(sigs)
!DEALLOCATE(s)
!DEALLOCATE(bcnvf)
!DEALLOCATE(psii) !,psio)
!DEALLOCATE(f)
IF (irank /= root) DEALLOCATE(sigt,STAT=sta)
!if (irank /= root) print *, sta
! End the PIDOSITS program
CALL MPI_FINALIZE(ierr)

STOP
END PROGRAM pidosits
