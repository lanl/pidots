MODULE itmvar

! Module to store the ITM variables

IMPLICIT NONE
SAVE

INTEGER, DIMENSION(:,:), ALLOCATABLE :: piv, fpiv

REAL*8, DIMENSION(:), ALLOCATABLE :: frphi
REAL*8, DIMENSION(:,:), ALLOCATABLE :: frpsi, foldpsi, phi, phiold, phir, sv, rphi, f, fr, fsv
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: psii, psio, av, oldpsi, rpsi, fpi, fpo, fav, jmat, fjmt
REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE :: kmat, fkmt, jpsi, fjps
REAL*8, DIMENSION(:,:,:,:,:), ALLOCATABLE :: kpsi, fkps

END MODULE
