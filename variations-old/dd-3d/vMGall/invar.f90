MODULE invar

! Module to store the AHOT input variables

IMPLICIT NONE
SAVE

! Parallel initialization
INTEGER :: isize, irank, root

! Title
CHARACTER(80) :: title

! Problem Size Specifications
INTEGER :: tpose, sym, idos, invf, pcf, pcty
INTEGER :: qdord, qdtyp, nxt, nyt, nzt, npx, npy, npz, nm
REAL*8 :: sorw

! Problem size per block
INTEGER :: nx, ny, nz
REAL*8, DIMENSION(:), ALLOCATABLE :: dx, dy, dz

! Cell materials
INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: mat

! Iteration Controls
REAL*8 :: err, cgr, tolr
INTEGER :: bitmxn, itmx, vitmx, bitmx, bitmx2

! Extra variables derived from input
INTEGER :: apo, warn

! Angular quadrature input
REAL*8, DIMENSION(:,:), ALLOCATABLE :: ang
REAL*8, DIMENSION(:), ALLOCATABLE :: w

! Boundary conditions by octant
! (Octants 1-4 like trig quadrants with positive z, 5-8 is negative z)
INTEGER, DIMENSION(6) :: bc
REAL*8, DIMENSION(:,:), ALLOCATABLE :: psii

! Cross section input
REAL*8, DIMENSION(:), ALLOCATABLE :: sigt, sigs

! Source data
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: s

! Print iteration data
INTEGER :: itp

! Print fluxes option
INTEGER :: sfp

! Parallel Topology setup
INTEGER :: allcomm, xcomm, ycomm, zcomm
INTEGER :: nrank, xrank, yrank, zrank
INTEGER, DIMENSION(:), ALLOCATABLE :: nxvec, nyvec, nzvec

END MODULE
