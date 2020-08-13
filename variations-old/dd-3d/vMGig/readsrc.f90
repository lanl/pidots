SUBROUTINE readsrc(srcfile,iexit)

!-------------------------------------------------------------
!
! Reads the source maps based on the format = 0/1
!
!-------------------------------------------------------------

USE invar
USE totvar
IMPLICIT NONE
INTEGER, INTENT(INOUT) :: iexit
INTEGER :: form, a, ni, i, j, k, g, xs, xe, ys, ye, zs, ze
REAL*8 :: q
CHARACTER(8), INTENT(IN) :: srcfile

! Set up the size of the arrays for the cross sections
ALLOCATE(st(nxt,nyt,nzt,ng))

! Initialize all elements of the source matrix to zero
st = 0.

! Open the source file for use
OPEN(UNIT=12, FILE=srcfile)

! First read the flag for the format of the source file
!   form = 0 or 1 only
!    0 => specification of source based on instructions
!         that relate to specific groups, moments, cells
!    1 => specification given for every group, moment, cell
READ(12,*)
READ(12,*) form

! Do the reading based on result
IF (form == 0) THEN
   ! Read the number of instructions, the group, the orders, the cells
   !  ni => number of instructions
   !  g  => the group
   !  xs => the starting x cell of the instruction
   !  ys => the starting y cell of the instruction
   !  zs => the starting z cell of the instruction
   !  xe => the ending x cell
   !  ye => the ending y cell
   !  ze => the ending z cell
   !  q  => temporary placeholder for source magnitude
   READ(12,*) ni
   IF (ni < 0) THEN
      WRITE(8,'(/,3X,A)') "ERROR: Illegal value # source instructions. Must be non-negative."
      iexit = iexit + 1
   END IF
   DO a = 1, ni
      READ(12,*) g
      READ(12,*) xs, ys, zs
      READ(12,*) xe, ye, ze
      READ(12,*) q
      DO k = zs, ze
         DO j = ys, ye
            DO i = xs, xe
               st(i,j,k,g) = q
            END DO
         END DO
      END DO
   END DO
   
! Or the reading is complete for all groups, moments, cells
ELSE IF (form == 1) THEN
   DO g = 1, ng
      READ(12,*)
      DO k = 1, nzt
         DO j = 1, nyt
            READ(12,*) (st(i,j,k,g), i = 1, nxt)
         END DO
      END DO
   END DO
! Only have two formats for the source file now
ELSE
   WRITE(8,*) "ERROR: Illegal value for the source file format, must be 0 or 1"
   iexit = iexit + 1
END IF

! Check the sources are all zero or greater
IF (MINVAL(st) < 0.) THEN
   WRITE(8,'(/,3X,A)') "ERROR: Illegal value for the source -- must be zero or greater"
   iexit = iexit + 1
END IF


RETURN
END SUBROUTINE readsrc
