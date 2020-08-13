MODULE solvar

! Module to store the solution variables

IMPLICIT NONE
SAVE

! Variables for subroutine 'weight'
REAL*8 :: ex

! Scalar flux final solution
REAL*8, DIMENSION(:,:), ALLOCATABLE :: flux

! Scalar flux vectors
REAL*8, DIMENSION(:,:), ALLOCATABLE :: f
REAL*8, DIMENSION(:), ALLOCATABLE :: fold

! Source vector with sources for all moments
REAL*8, DIMENSION(:,:), ALLOCATABLE :: sm

! Convergence flag
INTEGER, DIMENSION(:), ALLOCATABLE :: cnvf, bcnvf

! Gamma elements
REAL*8, DIMENSION(:), ALLOCATABLE :: gaa, gxa
REAL*8 :: gax, gxx

! X and Y matrices
REAL*8, DIMENSION(:,:), ALLOCATABLE :: xmat, xold

! Matrices for solving the system of equations
REAL*8, DIMENSION(:,:), ALLOCATABLE :: jmat
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: kmat
REAL*8, DIMENSION(:), ALLOCATABLE :: src, sv

! Variables for composing the kmats
REAL*8 :: xbc, xbco

! Matrices for computing the angular flux out of the domain at the boundaries
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: psio, jpsi
REAL*8, DIMENSION(:,:), ALLOCATABLE :: kpsi

END MODULE
