MODULE invar

! Module to store the AHOT input variables

IMPLICIT NONE
SAVE

! Title
CHARACTER(80) :: title

! Problem Size Specifications
INTEGER :: qdord, qdtyp, nxt, nyt, nzt, ng, nm

! Problem size per block
INTEGER :: nx, ny, nz
REAL*8, DIMENSION(:), ALLOCATABLE :: dx, dy, dz

! Cell materials
INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: mat

! Iteration Controls
REAL*8 :: err, tolr
INTEGER :: itmx

! Extra variables derived from input
INTEGER :: apo, warn
REAL*8, DIMENSION(:,:), ALLOCATABLE :: ssum

! Angular quadrature input
REAL*8, DIMENSION(:,:), ALLOCATABLE :: ang
REAL*8, DIMENSION(:), ALLOCATABLE :: w

INTEGER, DIMENSION(6) :: bc

! Cross section input
REAL*8, DIMENSION(:,:), ALLOCATABLE :: sigt
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: sigs

! Source data
REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE :: s

END MODULE
