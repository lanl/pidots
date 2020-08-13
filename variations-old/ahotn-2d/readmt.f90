SUBROUTINE readmt(mtfile)

!-------------------------------------------------------------
!
! Reads the material map from file
!
!-------------------------------------------------------------

USE invar
USE totvar
IMPLICIT NONE
INTEGER :: i, j
CHARACTER(8), INTENT(IN) :: mtfile

! Allocate the array size
ALLOCATE(matt(nxt,nyt))

! Open the file
OPEN(UNIT=13, FILE=mtfile)


! Read the dummy line
READ (13,*)

! Read by the y-row, then by the x-node
DO j = 1, nyt
   READ (13,*) (matt(i,j), i = 1, nxt)
END DO

RETURN
END SUBROUTINE readmt
