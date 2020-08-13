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
INTEGER :: bcs, n, indx0, i

! Allocate the size of the vectors for each quadrant from the dimensions and
! angular order of the problem
bcs = apo*(nxt+nyt)*order
ALLOCATE(psi1t(bcs),psi2t(bcs),psi3t(bcs),psi4t(bcs))

! Initialize all elements of the vectors
psi1t = 0.0
psi2t = 0.0
psi3t = 0.0
psi4t = 0.0

! Open the BC data file to be read
! Only need to open file if one of the quadrants has angular flux inward
IF (q1bc == 1 .OR. q2bc == 1 .OR. q3bc == 1 .OR. q4bc == 1) THEN
   OPEN(UNIT=14, FILE=bcfile)
   ! Read first dummy line of file
   READ(14,*)       ! Dummy line
END IF

! Read the BC data from the file
! The order goes by quadrant, angle/ybc then angle/xbc, increasing order
! So for each quadrant, read the ybc's (1:nx) for angle 1, repeat for angle 2, ...
! Then do the same for the xbc's (*Changed from before for message passing ease*)
IF (q1bc == 1) THEN        ! Vacuum BCs don't need to be read, just set to zero
   DO n = 1, apo
      indx0 = (n-1)*nxt*order
      READ(14,*)         ! Dummy line
      READ(14,*) (psi1t(i), i = (indx0+1), (indx0+nxt*order))
   END DO
   DO n = 1, apo
      indx0 = apo*nxt*order + (n-1)*nyt*order
      READ(14,*)         ! Dummy line
      READ(14,*) (psi1t(i), i = (indx0+1), (indx0+nyt*order))
   END DO
END IF
IF (q2bc == 1) THEN
   DO n = 1, apo
      indx0 = (n-1)*nxt*order
      READ(14,*)
      READ(14,*) (psi2t(i), i = (indx0+1), (indx0+nxt*order))
   END DO
   DO n = 1, apo
      indx0 = apo*nxt*order + (n-1)*nyt*order
      READ(14,*)
      READ(14,*) (psi2t(i), i = (indx0+1), (indx0+nyt*order))
   END DO
END IF
IF (q3bc == 1) THEN
   DO n = 1, apo
      indx0 = (n-1)*nxt*order
      READ(14,*)
      READ(14,*) (psi3t(i), i = (indx0+1), (indx0+nxt*order))
   END DO
   DO n = 1, apo
      indx0 = apo*nxt*order + (n-1)*nyt*order
      READ(14,*)
      READ(14,*) (psi3t(i), i = (indx0+1), (indx0+nyt*order))
   END DO
END IF
IF (q4bc == 1) THEN
   DO n = 1, apo
      indx0 = (n-1)*nxt*order
      READ(14,*)
      READ(14,*) (psi4t(i), i = (indx0+1), (indx0+nxt*order))
   END DO
   DO n = 1, apo
      indx0 = apo*nxt*order + (n-1)*nyt*order
      READ(14,*)
      READ(14,*) (psi4t(i), i = (indx0+1), (indx0+nyt*order))
   END DO
END IF

! Check that all the input angular fluxes are non-negative
IF (MINVAL(psi1t) < 0.0 .OR. MINVAL(psi2t) < 0.0 .OR. MINVAL(psi3t) < 0.0 .OR. MINVAL(psi4t) < 0.0) THEN
   WRITE(8,'(/,3X,A,/)') "WARNING: BC angular flux from input should be zero or greater"
END IF

CLOSE(14)

RETURN
END SUBROUTINE readbc
