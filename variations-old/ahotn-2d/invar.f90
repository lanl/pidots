MODULE invar

! Module to store the AHOT input variables

IMPLICIT NONE
SAVE

! Parallel initialization
INTEGER :: isize, irank

! Title
CHARACTER(80) :: title

! Problem Size Specifications
INTEGER :: lambda, meth, qdord, qdtyp, nxt, nyt, npx, npy, ng, nm
REAL*8, DIMENSION(:), ALLOCATABLE :: dx, dy

! Problem size per block
INTEGER :: nx, ny

! Cell materials
INTEGER, DIMENSION(:,:), ALLOCATABLE :: mat

! Iteration Controls
REAL*8 :: err, tolr
INTEGER :: bitmx, itmx, iall

! Solution check frequency
REAL*8 :: tchk
INTEGER :: ichk

! Extra variables derived from input
INTEGER :: apo, order, ordsq
REAL*8, DIMENSION(:,:), ALLOCATABLE :: ssum

! Angular quadrature input
REAL*8, DIMENSION(:,:), ALLOCATABLE :: ang
REAL*8, DIMENSION(:), ALLOCATABLE :: w

! Boundary Conditions by quadrant (quadrants 1-4 same as basic trig. quadrants)
INTEGER :: q1bc, q2bc, q3bc, q4bc
REAL*8, DIMENSION(:), ALLOCATABLE :: psi1, psi2, psi3, psi4

! Cross section input
REAL*8, DIMENSION(:,:), ALLOCATABLE :: sigt
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: sigs

! Source data
REAL*8, DIMENSION(:,:,:,:,:), ALLOCATABLE :: s

! Editing data
INTEGER :: momp, pmoaf

! Parallel Topology setup
INTEGER :: allcomm, xcomm, ycomm, nrank, xrank, yrank
INTEGER, DIMENSION(:), ALLOCATABLE :: nxvec, nyvec

! Solution methodology
INTEGER :: itmflag

! Print the matrix
INTEGER :: matrix

END MODULE
