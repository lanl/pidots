SUBROUTINE version

!-------------------------------------------------------------
!
!    Defines the code version, date and time of execution
!
!-------------------------------------------------------------

IMPLICIT NONE
CHARACTER(12) :: Real_clock (3)
INTEGER, DIMENSION (8) :: datetime

WRITE (8,*) "CODE: Parallel Integral Discrete Ordinates ", ACHAR(38), " Source Iteration Transport Solver - 3D"
WRITE (8,*) "VERSION: 2.9.11-DD-1D-AnisotropicScattering"
WRITE (8,*) "AUTHOR: Robert Joseph Zerr"
CALL DATE_AND_TIME (Real_clock(1), Real_clock(2), Real_clock(3), datetime)
WRITE (8, '(1X "Ran on " I2 "-" I2 "-" I4)', ADVANCE="NO") datetime(2), datetime(3), datetime(1)
WRITE (8, '(" at time " I2 ":" I2 ":" I2)') datetime(5), datetime(6), datetime(7)
WRITE (8,*)
WRITE (8,*) "-----------------------------------------------------------------------------"

131 FORMAT(A,A1,A)

RETURN
END SUBROUTINE version
