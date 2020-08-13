SUBROUTINE readbc(bcfile)

!-------------------------------------------------------------
!
! Reads the incoming angular flux from the input file
!
!-------------------------------------------------------------

USE invar
USE totvar
IMPLICIT NONE
CHARACTER(8), INTENT(IN) :: bcfile
INTEGER :: n

! Allocate the size of the vectors for each quadrant from the dimensions and
! angular order of the problem
ALLOCATE(posinx(apo), neginx(apo))

posinx = 0.0
neginx = 0.0

! Open the BC data file to be read
! Only need to open file if one of the quadrants has angular flux inward
IF (MAXVAL(bc) == 2) THEN
   OPEN(UNIT=14, FILE=bcfile)
   ! Read first dummy line of file
   READ(14,*)       ! Dummy line
END IF

! Read the BC data from the file
! Eventually these matrices will be put into single vectors.
! Within each set of those it goes by angle, then by cell.
! There will be a vector for each direction but held in single array.
! For now, keep it as a vector for each end of slab: +X, -X
IF (bc(1) == 2) THEN
   DO n = 1, apo
      READ(14,*)                ! Dummy line
      READ(14,*) posinx(n)
   END DO
END IF
IF (bc(2) == 2) THEN
   DO n = 1, apo
      READ(14,*)                ! Dummy line
      READ(14,*) neginx(n)
   END DO
END IF

! Check that all the input angular fluxes are non-negative
IF (MINVAL(posinx) < 0.0 .OR. MINVAL(neginx) < 0.0) THEN
   warn = warn + 1
   WRITE(8,'(/,3X,A,/)') "WARNING: XBC angular flux from input should be zero or greater"
END IF

IF (MAXVAL(bc) == 2) CLOSE(14)

RETURN
END SUBROUTINE readbc
