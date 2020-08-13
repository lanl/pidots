MODULE timevar

! Module to store the timing variables

IMPLICIT NONE
SAVE

! Called in MAIN
! Final time of the program
REAL :: tend

! Called in SOLVE
! Time to reach solution phase, time after construction of matrices, after system solved
REAL :: ttosolve, tjmat, tsolve

END MODULE
