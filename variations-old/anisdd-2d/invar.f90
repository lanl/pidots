MODULE invar

! Module to store the AHOT input variables

IMPLICIT NONE
SAVE

! Parallel initialization
INTEGER :: isize, irank, root

! Title
CHARACTER(80) :: title

! Problem Size Specifications
INTEGER :: meth, tpose, sym, idos, qdord, qdtyp, anord, nxt, nyt, npx, npy, ng, nm

! Problem size per block
INTEGER :: nx, ny
REAL*8, DIMENSION(:), ALLOCATABLE :: dx, dy

! Cell materials
INTEGER, DIMENSION(:,:), ALLOCATABLE :: mat

! Iteration Controls
REAL*8 :: err, cgr, tolr
INTEGER :: bitmx, itmx

! Solution check frequency
REAL*8 :: tchk
INTEGER :: ichk

! Extra variables derived from input
INTEGER :: apo, warn, sord, nmom
REAL*8, DIMENSION(:,:), ALLOCATABLE :: ssum

! Angular quadrature input
REAL*8, DIMENSION(:,:), ALLOCATABLE :: ang
REAL*8, DIMENSION(:), ALLOCATABLE :: omg, w

! Boundary conditions by octant
! (Octants 1-4 like trig quadrants with positive z)
INTEGER, DIMENSION(4) :: bc
REAL*8, DIMENSION(:,:), ALLOCATABLE :: psii

! Cross section input
REAL*8, DIMENSION(:,:), ALLOCATABLE :: sigt
REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE :: sigs

! Source data
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: s

! Print iteration data
INTEGER :: itp

! Print fluxes option
INTEGER :: sfp

! Parallel Topology setup
INTEGER :: allcomm, xcomm, ycomm
INTEGER :: nrank, xrank, yrank
INTEGER, DIMENSION(:), ALLOCATABLE :: nxvec, nyvec

END MODULE
