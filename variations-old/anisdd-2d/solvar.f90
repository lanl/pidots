MODULE solvar

! Module to store the solution variables

IMPLICIT NONE
SAVE

! Variables for subroutine 'weight'
REAL*8 :: ex, ey

! Scalar flux final solution
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: flux

! Scalar flux moments: SI
REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE :: f
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: e, fold

! Source moments
REAL*8, DIMENSION(:,:), ALLOCATABLE :: sm

! Scalar flux moments: ITM
REAL*8, DIMENSION(:,:), ALLOCATABLE :: phi
REAL*8, DIMENSION(:), ALLOCATABLE :: phiold

! Convergence flag
INTEGER, DIMENSION(:), ALLOCATABLE :: cnvf, bcnvf

! Gamma elements
REAL*8, DIMENSION(:), ALLOCATABLE :: gaa, gxa, gya
REAL*8 :: gax, gay, gxx, gxy, gyx, gyy

! X and Y matrices
REAL*8, DIMENSION(:,:,:), ALLOCATABLE     :: xmat, xold, yold
REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE   :: ymat

! Matrices for solving the system of equations
REAL*8, DIMENSION(:,:), ALLOCATABLE :: jmat
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: kmat
REAL*8, DIMENSION(:), ALLOCATABLE :: src, sv

! Matrices for composing the kmats
REAL*8, DIMENSION(:,:), ALLOCATABLE :: ybcx, xbcx
REAL*8, DIMENSION(:), ALLOCATABLE :: ybcy, ybcyo, ybcxo
REAL*8, DIMENSION(:), ALLOCATABLE :: xbcy, xbcyo, xbcxo

! Matrices for computing the angular flux out of the domain at the boundaries
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: psio, jpsi
REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE :: kpsi

END MODULE
