MODULE solvar

! Module to store the solution variables

IMPLICIT NONE
SAVE

! Variables for allocating operators' sizes
INTEGER :: neq, bcs, bcs2, xys, xzs, yzs
INTEGER :: cneq, cbcs, cbcs2, cxys, cxzs, cyzs

! Variables for subroutine 'weight'
REAL*8 :: ex, ey, ez

! Scalar flux final solution
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: flux

! Scalar flux moments and Angular flux out
REAL*8, DIMENSION(:), ALLOCATABLE :: phi, phiold, f, fold, phir
REAL*8, DIMENSION(:,:), ALLOCATABLE :: oldpsi, psio, cpsii, cpsio, rrpsi

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
INTEGER, DIMENSION(:), ALLOCATABLE :: piv, cpiv
REAL*8, DIMENSION(:,:), ALLOCATABLE :: jmat, cjmat
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: kmat, ckmat
REAL*8, DIMENSION(:), ALLOCATABLE :: src, sv, dv

! Matrices for composing the kmats
REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE :: zbcxy, ybcxy, xbcxy
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: zbcxz, ybcxz, xbcxz
REAL*8, DIMENSION(:,:), ALLOCATABLE :: zbcyz, zbcxyo, zbcxzo, zbcyzo
REAL*8, DIMENSION(:,:), ALLOCATABLE :: ybcyz, ybcxyo, ybcxzo, ybcyzo
REAL*8, DIMENSION(:,:), ALLOCATABLE :: xbcyz, xbcxyo, xbcxzo, xbcyzo

! Matrices for computing the angular flux out of the domain at the boundaries
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: jpsi, cjpsi
REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE :: kpsi, ckpsi

END MODULE
