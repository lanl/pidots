MODULE solvar

! Module to store the solution variables

IMPLICIT NONE
SAVE

! Variables for allocating operators' sizes
INTEGER :: neq, bcs, bcs2, xys, xzs, yzs

! Variables for subroutine 'weight'
REAL*8 :: ex, ey, ez

! Scalar flux final solution
REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE :: flux
REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE :: afx, afy, afz

! Scalar flux moments: SI
REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE :: f
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: fold

! Scalar flux moments: ITM
REAL*8, DIMENSION(:,:), ALLOCATABLE :: phi
REAL*8, DIMENSION(:), ALLOCATABLE :: phiold

! Convergence flag
INTEGER, DIMENSION(:), ALLOCATABLE :: cnvf, bcnvf

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
INTEGER, DIMENSION(:), ALLOCATABLE :: piv
REAL*8, DIMENSION(:,:), ALLOCATABLE :: jmat
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: kmat
REAL*8, DIMENSION(:), ALLOCATABLE :: src, sv, dv

! Matrices for composing the kmats
REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE :: zbcxy, ybcxy, xbcxy
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: zbcxz, ybcxz, xbcxz
REAL*8, DIMENSION(:,:), ALLOCATABLE :: zbcyz, zbcxyo, zbcxzo, zbcyzo
REAL*8, DIMENSION(:,:), ALLOCATABLE :: ybcyz, ybcxyo, ybcxzo, ybcyzo
REAL*8, DIMENSION(:,:), ALLOCATABLE :: xbcyz, xbcxyo, xbcxzo, xbcyzo

! Matrices for computing the angular flux out of the domain at the boundaries
REAL*8, DIMENSION(:,:), ALLOCATABLE :: rpsi
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: psio, jpsi
REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE :: kpsi

END MODULE
