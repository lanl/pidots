MODULE invar

! Module to store the AHOT input variables

IMPLICIT NONE
SAVE

! Parallel initialization
INTEGER :: isize, irank, root

! Title
CHARACTER(80) :: title

! Problem Size Specifications
INTEGER :: qdord, qdtyp, nxt, nyt, nzt, npx, npy, npz, nm

! Problem size per block
INTEGER :: nx, ny, nz, sdnx, sdny, sdnz, nxsd, nysd, nzsd, nsdp
INTEGER :: neq, xys, xzs, yzs, bcs, bcs2
REAL*8, DIMENSION(:), ALLOCATABLE :: dx, dy, dz

! Cell materials
INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: mat

! Iteration Controls
REAL*8 :: err, tolr
INTEGER :: bitmx

! Extra variables derived from input
INTEGER :: apo, warn

! Angular quadrature input
REAL*8, DIMENSION(:,:), ALLOCATABLE :: ang
REAL*8, DIMENSION(:), ALLOCATABLE :: w

! Boundary conditions by octant
! (Octants 1-4 like trig quadrants with positive z, 5-8 is negative z)
INTEGER, DIMENSION(6) :: bc
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: psii

! Cross section input
REAL*8, DIMENSION(:), ALLOCATABLE :: sigt
REAL*8, DIMENSION(:), ALLOCATABLE :: sigs

! Source data
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: s

! Print iteration data
INTEGER :: itp

! Print fluxes option
INTEGER :: sfp

! Parallel Topology setup
INTEGER :: allcomm, xcomm, ycomm, zcomm
INTEGER :: nrank, xrank, yrank, zrank

END MODULE
