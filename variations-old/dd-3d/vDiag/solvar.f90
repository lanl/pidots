MODULE solvar

! Module to store the solution variables

IMPLICIT NONE
SAVE

! Variables for subroutine 'weight'
REAL*8 :: ex, ey, ez

! Scalar flux moments: ITM
REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE :: phi

! Convergence flag
INTEGER, DIMENSION(:), ALLOCATABLE :: cnvf

! Gamma elements
REAL*8 :: gaa,  gaxy,  gaxz,  gayz
REAL*8 :: gxya, gxyxy, gxyxz, gxyyz
REAL*8 :: gxza, gxzxy, gxzxz, gxzyz
REAL*8 :: gyza, gyzxy, gyzxz, gyzyz

! X and Y matrices
REAL*8 :: xmat, xold, yold, zold
REAL*8, DIMENSION(:), ALLOCATABLE   :: ymat
REAL*8, DIMENSION(:,:), ALLOCATABLE :: zmat

! Matrices for solving the system of equations
REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE :: jmat
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: src, sv

END MODULE
