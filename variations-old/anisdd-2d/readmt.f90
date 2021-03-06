SUBROUTINE readmt(mtfile,iexit)

!-------------------------------------------------------------
!
! Reads the material map from file as either all cells or
!  by instructions like the source
!
!-------------------------------------------------------------

USE invar
USE totvar
IMPLICIT NONE
INTEGER, INTENT(INOUT) :: iexit
INTEGER :: form, ni, a, i, j, m, xs, xe, ys, ye
CHARACTER(8), INTENT(IN) :: mtfile

! Allocate the array size
ALLOCATE(matt(nxt,nyt))

! Open the file
OPEN(UNIT=13, FILE=mtfile)

! Read the dummy line
READ (13,*)

! First read the flag for the format of the material map file
!  form = 0 or 1 only
!  0 => specification of material map based on instructions
!       first set everything to base material '1' then use instructions to
!       overwrite
!  1 => specification of material map given for every cell
READ(13,*) form

! Do the reading based on the result
IF (form == 0) THEN
   ! Read the number of instructions then the starting cells then the ending cells 
   ! ni = number of instructions
   ! x/ys = starting x/y cell
   ! x/ye = ending x/y cell
   
   ! Set the entire material map to 1
   matt = 1

   ! Read the number of instructions to define other material locations
   READ(13,*) ni
   IF (ni < 0) THEN
      WRITE(8,'(/,3X,A)') "ERROR: Illegal value for # material instructions. Must be non-negative."
      iexit = iexit + 1
   END IF
   DO a = 1, ni
      READ(13,*) xs, ys
      READ(13,*) xe, ye
      READ(13,*) m
      DO j = ys, ye
         DO i = xs, xe
            matt(i,j) = m
         END DO
      END DO
   END DO

ELSE IF (form == 1) THEN
   ! Read by the y-row, then by the x-node
   DO j = 1, nyt
      READ (13,*) (matt(i,j), i = 1, nxt)
   END DO

ELSE
   WRITE(8,*) "ERROR: Illegal value for the material file format, must be 0 or 1"
   iexit = iexit + 1
END IF

RETURN
END SUBROUTINE readmt
