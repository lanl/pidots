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
INTEGER :: n, i, o

! Allocate the size of the vectors for each quadrant from the dimensions and
! angular order of the problem
ALLOCATE(posiny(nxt,apo,2), neginy(nxt,apo,2))
ALLOCATE(posinx(nyt,apo,2), neginx(nyt,apo,2))

posiny = 0.0
neginy = 0.0
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
! Eventually these matrices will be put into single vectors that go ybc,
! xbc. Within each set of those it goes by angle, then by cell.
! There will be a vector for each quadrant.
! For now, keep it as a matrix for each face: +Y, -Y, +X, -X
! Read values for the four quadrants of each face in increasing order of 1-4. 1-4 are positive
! z direction + trig setup of quadrants.
IF (bc(1) == 2) THEN
   DO o = 1, 2     ! Actual octants 1, 2
      DO n = 1, apo
         READ(14,*)
         DO i = 1, nxt
            READ(14,*) posiny(i,n,o)
         END DO
      END DO
   END DO
END IF
IF (bc(2) == 2) THEN
   DO o = 1, 2     ! Actual octants 3, 4
      DO n = 1, apo
         READ(14,*)
         DO i = 1, nxt
            READ(14,*) neginy(i,n,o)
         END DO
      END DO
   END DO
END IF
IF (bc(3) == 2) THEN
   DO o = 1, 2     ! Actual octants 1, 4
      DO n = 1, apo
         READ(14,*)
         DO i = 1, nyt
            READ(14,*) posinx(i,n,o)
         END DO
      END DO
   END DO
END IF
IF (bc(4) == 2) THEN
   DO o = 1, 2     ! Actual octants 2, 3
      DO n = 1, apo
         READ(14,*)
         DO i = 1, nyt
            READ(14,*) neginx(i,n,o)
         END DO
      END DO
   END DO
END IF

! Check that all the input angular fluxes are non-negative
IF (MINVAL(posiny) < 0.0 .OR. MINVAL(neginy) < 0.0) THEN
   warn = warn + 1
   WRITE(8,'(/,3X,A,/)') "WARNING: YBC angular flux from input should be zero or greater"
END IF
IF (MINVAL(posinx) < 0.0 .OR. MINVAL(neginx) < 0.0) THEN
   warn = warn + 1
   WRITE(8,'(/,3X,A,/)') "WARNING: XBC angular flux from input should be zero or greater"
END IF

CLOSE(14)

RETURN
END SUBROUTINE readbc
