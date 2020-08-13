PROGRAM pidosits

!-------------------------------  PIDOSITS.1.3D  -----------------------------------
!
!  Parallel Integral Discrete Ordinates & Source Iteration Transport Solver  
!  By RJZerr, April 2010
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
!    Special version that constructs Jmat with only nearest neighbor elements,
!    not full, global coupling
!
!
!-----------------------------------------------------------------------------------

USE invar
USE solvar
USE timevar
IMPLICIT NONE
CHARACTER(8) :: infile, outfile, qdfile, xsfile, srcfile, mtfile, bcfile
INTEGER :: statin, statout, ferr, ierr, iexit
LOGICAL :: existence

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

IF (ferr == 1) THEN
   STOP
END IF

! Set up the introductory info
CALL version

! Read input data:
CALL input(qdfile, xsfile, srcfile, mtfile, bcfile, iexit)

! Stop all processes if there was an error detected. 0 process writes out error messages
IF (iexit > 0) THEN
   STOP
END IF

! Echo the input data:
CALL echo(infile, outfile, qdfile, xsfile, srcfile, mtfile, bcfile)

! Solve the transport problem:
CALL solve

! Print the output
CALL output

CALL CPU_TIME(tend)

! Print the relevant times of the execution
WRITE(8,'(/,2X,A,/)') "Fortran95 Timing Estimates with CPU_TIME..."
WRITE(8,100)
WRITE(8,101) "InputWrk", ttosolve, ttosolve
WRITE(8,101) "IDOmatrx", tjmat, tjmat-ttosolve
WRITE(8,101) "IDOitrns", tsolve, tsolve-tjmat
WRITE(8,101) "SolveTot", tsolve, tsolve-ttosolve
WRITE(8,101) "PrintOut", tend, tend-tsolve
WRITE(8,102)
100 FORMAT(5X,'WorkDone',3X,'Absolute(s)',6X,'Difference(s)')
101 FORMAT(5X,A8,3X,F9.3,5X,F9.3)
102 FORMAT(//,'*********************   END PROGRAM   ************************')

! Deallocate the allocated arrays
DEALLOCATE(dx,dy,dz,mat,ang,w,sigt,sigs,s)
DEALLOCATE(phi)

STOP
END PROGRAM pidosits
