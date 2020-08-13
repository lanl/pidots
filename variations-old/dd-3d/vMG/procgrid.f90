SUBROUTINE procgrid

!-------------------------------------------------------------
!
! Establishes the 3D cartesian process topology
!
!-------------------------------------------------------------

USE invar
IMPLICIT NONE
INTEGER, DIMENSION(3) :: dims
INTEGER :: ierr
LOGICAL :: reorder
LOGICAL, DIMENSION(3) :: periodic, remain

INCLUDE 'mpif.h'

! Set up the number of dimensions and create topology
dims(1) = npz
dims(2) = npy
dims(3) = npx
periodic(1) = .FALSE.
periodic(2) = .FALSE.
periodic(3) = .FALSE.
reorder = .TRUE.
CALL MPI_CART_CREATE(MPI_COMM_WORLD,3,dims,periodic,reorder,allcomm,ierr)

! Once created, suppress two indices at a time to create row/col/stack communicators
remain(1) = .FALSE.
remain(2) = .FALSE.
remain(3) = .TRUE.
CALL MPI_CART_SUB(allcomm,remain,xcomm,ierr)

remain(1) = .FALSE.
remain(2) = .TRUE.
remain(3) = .FALSE.
CALL MPI_CART_SUB(allcomm,remain,ycomm,ierr)

remain(1) = .TRUE.
remain(2) = .FALSE.
remain(3) = .FALSE.
CALL MPI_CART_SUB(allcomm,remain,zcomm,ierr)

! Get the ranks for the new communicators
CALL MPI_COMM_RANK(allcomm,nrank,ierr)
CALL MPI_COMM_RANK(xcomm,xrank,ierr)
CALL MPI_COMM_RANK(ycomm,yrank,ierr)
CALL MPI_COMM_RANK(zcomm,zrank,ierr)

RETURN
END SUBROUTINE procgrid 
