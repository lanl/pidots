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
INTEGER :: n, i, j, o

! Allocate the size of the vectors for each quadrant from the dimensions and
! angular order of the problem
ALLOCATE(posinz(nxt,nyt,apo,4), neginz(nxt,nyt,apo,4))
ALLOCATE(posiny(nxt,nzt,apo,4), neginy(nxt,nzt,apo,4))
ALLOCATE(posinx(nyt,nzt,apo,4), neginx(nyt,nzt,apo,4))

posinz = 0.0d0
neginz = 0.0d0
posiny = 0.0d0
neginy = 0.0d0
posinx = 0.0d0
neginx = 0.0d0

! Open the BC data file to be read
! Only need to open file if one of the quadrants has angular flux inward
IF (MAXVAL(bc) == 2) THEN
   OPEN(UNIT=14, FILE=bcfile)
   ! Read first dummy line of file
   READ(14,*)       ! Dummy line
END IF

! Read the BC data from the file
! Eventually these matrices will be put into single vectors that go zbc, ybc,
! xbc. Within each set of those it goes by angle, then by cell.
! There will be a vector for each octant.
! For now, keep it as a matrix for each face: +Z, -Z, +Y, -Y, +X, -X
! Read values for the four octants of each face in increasing order of 1-8. 1-4 are positive
! z direction + trig setup of quadrants. 5-8 are negative.
IF (bc(1) == 2) THEN     ! Have incoming flux on that face, vacuum remains zero
   DO o = 1, 4     ! Actual octants 1, 2, 3, 4
      DO n = 1, apo
         READ(14,*)      ! Dummy line with face, octant, angle
         DO j = 1, nyt
            DO i = 1, nxt
               READ(14,*) posinz(i,j,n,o)
            END DO
         END DO
      END DO
   END DO
END IF
IF (bc(2) == 2) THEN
   DO o = 1, 4     ! Actual octants 5, 6, 7, 8
      DO n = 1, apo
         READ(14,*)
         DO j = 1, nyt
            DO i = 1, nxt
               READ(14,*) neginz(i,j,n,o)
            END DO
         END DO
      END DO
   END DO
END IF
IF (bc(3) == 2) THEN
   DO o = 1, 4     ! Actual octants 1, 2, 5, 6
      DO n = 1, apo
         READ(14,*)
         DO j = 1, nzt
            DO i = 1, nxt
               READ(14,*) posiny(i,j,n,o)
            END DO
         END DO
      END DO
   END DO
END IF
IF (bc(4) == 2) THEN
   DO o = 1, 4     ! Actual octants 3, 4, 7, 8
      DO n = 1, apo
         READ(14,*)
         DO j = 1, nzt
            DO i = 1, nxt
               READ(14,*) neginy(i,j,n,o)
            END DO
         END DO
      END DO
   END DO
END IF
IF (bc(5) == 2) THEN
   DO o = 1, 4     ! Actual octants 1, 4, 5, 8
      DO n = 1, apo
         READ(14,*)
         DO j = 1, nzt
            DO i = 1, nyt
               READ(14,*) posinx(i,j,n,o)
            END DO
         END DO
      END DO
   END DO
END IF
IF (bc(6) == 2) THEN
   DO o = 1, 4     ! Actual octants 2, 3 , 6, 7
      DO n = 1, apo
         READ(14,*)
         DO j = 1, nzt
            DO i = 1, nyt
               READ(14,*) neginx(i,j,n,o)
            END DO
         END DO
      END DO
   END DO
END IF

! Check that all the input angular fluxes are non-negative
IF (MINVAL(posinz) < 0.0 .OR. MINVAL(neginz) < 0.0) THEN
   warn = warn + 1
   WRITE(8,'(/,3X,A,/)') "WARNING: ZBC angular flux from input should be zero or greater"
END IF
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
