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
INTEGER :: form, a, ni, i, j, k, l, g, xs, xe, ys, ye
REAL*8 :: q
CHARACTER(8), INTENT(IN) :: srcfile

! Set up the size of the arrays for the cross sections
ALLOCATE(st(nxt,nyt,0:lambda,0:lambda,ng))

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
   !  k  => x-moment (spatial) of the instruction
   !  l  => y-moment of the instruction
   !  xs => the starting x cell of the instruction
   !  ys => the starting y cell of the instruction
   !  xe => the ending x cell
   !  ye => the ending y cell
   !  q  => temporary placeholder for source magnitude
   READ(12,*) ni
   IF (ni <= 0) THEN
      WRITE(8,'(/,3x,A)') "ERROR: Illegal value # source instructions. Must be positive."
      iexit = iexit + 1
   END IF
   DO a = 1, ni
      READ(12,*) g, k, l
      READ(12,*) xs, ys
      READ(12,*) xe, ye
      READ(12,*) q
      DO j = ys, ye
         DO i = xs, xe
            st(i,j,k,l,g) = q
         END DO
      END DO
   END DO
   
! Or the reading is complete for all groups, moments, cells
ELSE IF (form == 1) THEN
   DO g = 1, ng
      DO k = 0, lambda
         DO l = 0, lambda
            READ(12,*)
            DO j = 1, nyt
               READ(12,*) (st(i,j,k,l,g), i = 1, nxt)
            END DO
         END DO
      END DO
   END DO
! Only have two formats for the source file now
ELSE
   WRITE(8,*) "ERROR: Illegal value for the source file format, must be 1 or 2"
   iexit = iexit + 1
END IF

! Check the sources are all zero or greater
IF (MINVAL(st) < 0.) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal value for the source -- must be zero or greater"
   iexit = iexit + 1
END IF

RETURN
END SUBROUTINE readsrc
