SUBROUTINE readmt(mtfile,iexit)

!-------------------------------------------------------------
!
! Reads the material map from file as either all cells or
!  by instructions like the source
!
!-------------------------------------------------------------

USE invar
IMPLICIT NONE
INTEGER, INTENT(INOUT) :: iexit
INTEGER :: form, ni, a, i, j, k, m, xs, xe, ys, ye, zs, ze
CHARACTER(8), INTENT(IN) :: mtfile

! Allocate the array size
ALLOCATE(mat(nx,ny,nz))

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
   ! x/y/zs = starting x/y/z cell
   ! x/y/ze = ending x/y/z cell
   
   ! Set the entire material map to 1
   mat = 1

   ! Read the number of instructions to define other material locations
   READ(13,*) ni
   IF (ni < 0) THEN
      WRITE(8,'(/,3X,A)') "ERROR: Illegal value for # material instructions. Must be non-negative."
      iexit = iexit + 1
   END IF
   DO a = 1, ni
      READ(13,*) xs, ys, zs
      READ(13,*) xe, ye, ze
      READ(13,*) m
      DO k = zs, ze
         DO j = ys, ye
            DO i = xs, xe
               mat(i,j,k) = m
            END DO
         END DO
      END DO
   END DO

ELSE IF (form == 1) THEN
   ! Read by the z-plane, then by the y-row, then by the x-node
   DO k = 1, nz
      READ (13,*)
      DO j = 1, ny
         READ (13,*) (mat(i,j,k), i = 1, nx)
      END DO
   END DO

ELSE
   WRITE(8,*) "ERROR: Illegal value for the material file format, must be 0 or 1"
   iexit = iexit + 1
END IF

RETURN
END SUBROUTINE readmt
