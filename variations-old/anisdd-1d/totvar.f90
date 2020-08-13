MODULE totvar

! Module for the variables for the entire sub-domian

IMPLICIT NONE
SAVE

! Problem size
REAL*8, DIMENSION(:), ALLOCATABLE :: dxt

! Cell materials
INTEGER, DIMENSION(:), ALLOCATABLE :: matt

! Boundary conditions
REAL*8, DIMENSION(:), ALLOCATABLE :: posinx, neginx

! Source data
REAL*8, DIMENSION(:,:), ALLOCATABLE :: st

END MODULE
