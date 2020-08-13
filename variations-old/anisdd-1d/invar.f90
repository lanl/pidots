MODULE invar

! Module to store the AHOT input variables

IMPLICIT NONE
SAVE

! Parallel initialization
INTEGER :: isize, irank, root

! Title
CHARACTER(80) :: title

! Problem Size Specifications
INTEGER :: meth, tpose, sym, idos, qdord, qdtyp, anord, nxt, npx, ng, nm

! Problem size per block
INTEGER :: nx
REAL*8, DIMENSION(:), ALLOCATABLE :: dx

! Cell materials
INTEGER, DIMENSION(:), ALLOCATABLE :: mat

! Iteration Controls
REAL*8 :: err, cgr, tolr
INTEGER :: bitmx, itmx

! Solution check frequency
REAL*8 :: tchk
INTEGER :: ichk

! Extra variables derived from input
INTEGER :: apo, warn, sord
REAL*8, DIMENSION(:,:), ALLOCATABLE :: ssum

! Angular quadrature input
REAL*8, DIMENSION(:), ALLOCATABLE :: ang, w

! Boundary conditions: positive x, negative x
INTEGER, DIMENSION(2) :: bc
REAL*8, DIMENSION(:,:), ALLOCATABLE :: psii

! Cross section input
REAL*8, DIMENSION(:,:), ALLOCATABLE :: sigt
REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE :: sigs

! Source data
REAL*8, DIMENSION(:,:), ALLOCATABLE :: s

! Print iteration data
INTEGER :: itp

! Print fluxes option
INTEGER :: sfp

! Parallel Topology setup
INTEGER, DIMENSION(:), ALLOCATABLE :: nxvec

END MODULE
