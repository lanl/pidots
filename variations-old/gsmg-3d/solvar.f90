MODULE solvar

! Module to store the matsweep variables

IMPLICIT NONE
SAVE

! Variables for subroutine 'weight'
REAL*8 :: ex, ey, ez

! Gamma elements
REAL*8 :: gaa,  gaxy,  gaxz,  gayz
REAL*8 :: gxya, gxyxy, gxyxz, gxyyz
REAL*8 :: gxza, gxzxy, gxzxz, gxzyz
REAL*8 :: gyza, gyzxy, gyzxz, gyzyz

! X and Y matrices
REAL*8, DIMENSION(:,:,:), ALLOCATABLE     :: xmat, xold, yold, zold
REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE   :: ymat
REAL*8, DIMENSION(:,:,:,:,:), ALLOCATABLE :: zmat

! Matrices for composing the kmats
REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE :: zbcxy, ybcxy, xbcxy
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: zbcxz, ybcxz, xbcxz
REAL*8, DIMENSION(:,:), ALLOCATABLE :: zbcyz, zbcxyo, zbcxzo, zbcyzo
REAL*8, DIMENSION(:,:), ALLOCATABLE :: ybcyz, ybcxyo, ybcxzo, ybcyzo
REAL*8, DIMENSION(:,:), ALLOCATABLE :: xbcyz, xbcxyo, xbcxzo, xbcyzo

END MODULE
