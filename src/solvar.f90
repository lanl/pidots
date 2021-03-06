MODULE solvar

  ! Module to store the solution variables

  IMPLICIT NONE
  SAVE

  ! Variables for subroutine 'weight'
  REAL*8 :: ex, ey, ez

  REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: flux, fluxt

  ! Scalar flux moments: ITM
  REAL*8, DIMENSION(:,:), ALLOCATABLE :: phi, phiold

  ! Convergence flag
  INTEGER :: bcnvf

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
  REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: jmat
  REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE :: kmat
  REAL*8, DIMENSION(:,:), ALLOCATABLE :: src, sv

  ! Matrices for composing the kmats
  REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE :: zbcxy, ybcxy, xbcxy
  REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: zbcxz, ybcxz, xbcxz
  REAL*8, DIMENSION(:,:), ALLOCATABLE :: zbcyz, zbcxyo, zbcxzo, zbcyzo
  REAL*8, DIMENSION(:,:), ALLOCATABLE :: ybcyz, ybcxyo, ybcxzo, ybcyzo
  REAL*8, DIMENSION(:,:), ALLOCATABLE :: xbcyz, xbcxyo, xbcxzo, xbcyzo

  ! Matrices for computing the angular flux out of the domain at the boundaries
  REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: psio
  REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE :: jpsi
  REAL*8, DIMENSION(:,:,:,:,:), ALLOCATABLE :: kpsi

  !TODO REMOVE TIMING STUFF BELOW
  DOUBLE PRECISION:: com_max(1000,1000)=0, dtime, stime, etime
  DOUBLE PRECISION:: com_std(1000,1000)=0, com_avg(1000,1000)=0
  DOUBLE PRECISION:: com_temp_R(1000, 1000)=0
  INTEGER:: com_temp_I(1000,1000)=0, com_count(1000,1000)=0
  INTEGER:: s_indx, r_indx
  !END TODO

END MODULE solvar
