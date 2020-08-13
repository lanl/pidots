MODULE solvar

! Module to store the solution variables

IMPLICIT NONE
SAVE

! Variables for allocating operators' sizes
INTEGER :: neq, bcs, bcs2, xys, xzs, yzs

! Variables for subroutine 'weight'
REAL*8 :: ex, ey, ez

! Scalar flux final solution
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: flux

! Scalar flux moments and Angular flux out
REAL*8, DIMENSION(:), ALLOCATABLE :: phiold
REAL*8, DIMENSION(:,:), ALLOCATABLE :: phi, phir, oldpsi
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: psio
REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE :: phir2

! Gamma elements
REAL*8 :: gaa,  gaxy,  gaxz,  gayz
REAL*8 :: gxya, gxyxy, gxyxz, gxyyz
REAL*8 :: gxza, gxzxy, gxzxz, gxzyz
REAL*8 :: gyza, gyzxy, gyzxz, gyzyz

! X and Y matrices
REAL*8, DIMENSION(:,:,:), ALLOCATABLE     :: xmat, xold, yold, zold
REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE   :: ymat
REAL*8, DIMENSION(:,:,:,:,:), ALLOCATABLE :: zmat

! Matrices for solving the system of equations
INTEGER, DIMENSION(:,:), ALLOCATABLE :: piv
REAL*8, DIMENSION(:,:), ALLOCATABLE :: sv
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: av, jmat
REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE :: kmat

! Matrices for composing the kmats
REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE :: zbcxy, ybcxy, xbcxy
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: zbcxz, ybcxz, xbcxz
REAL*8, DIMENSION(:,:), ALLOCATABLE :: zbcyz, zbcxyo, zbcxzo, zbcyzo
REAL*8, DIMENSION(:,:), ALLOCATABLE :: ybcyz, ybcxyo, ybcxzo, ybcyzo
REAL*8, DIMENSION(:,:), ALLOCATABLE :: xbcyz, xbcxyo, xbcxzo, xbcyzo

! Matrices for computing the angular flux out of the domain at the boundaries
REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE :: jpsi
REAL*8, DIMENSION(:,:,:,:,:), ALLOCATABLE :: kpsi

END MODULE
