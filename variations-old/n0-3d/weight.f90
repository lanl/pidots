SUBROUTINE weight

!-------------------------------------------------------------
!
!  Computes the spatial weights
!  From Azmy 1988 in NSE, equations 20a and 20b used
!  First find the epsilons
!  Input: order, dx, dy, dz, total XS, mu, eta, xi
!  Returns: weight for cell (i,j,k) in x (alpha) and 
!            y (beta) directions and z (gamma) directions
!
!-------------------------------------------------------------

USE solvar
IMPLICIT NONE
REAL*8 :: ext, eyt, ezt, shx, shy, shz, chx, chy, chz, smxe, smye, smze
REAL*8 :: alpha, beta, gamma

! Compute the coefficients and the hyperbolic functions
ext = 0.5/ex
eyt = 0.5/ey
ezt = 0.5/ez
shx = SINH(ext)
shy = SINH(eyt)
shz = SINH(ezt)
chx = COSH(ext)
chy = COSH(eyt)
chz = COSH(ezt)

! Compute the even summation term
smxe = shx/ext
smye = shy/eyt
smze = shz/ezt

! Set the weights
alpha = (chx - smxe)/shx
beta  = (chy - smye)/shy
gamma = (chz - smze)/shz

! Instead of using the weights, store the weighted factors for in/out faces
ai =  0.5*(1.0-alpha)
ao = -0.5*(1.0+alpha)
bi =  0.5*(1.0-beta)
bo = -0.5*(1.0+beta)
gi =  0.5*(1.0-gamma)
go = -0.5*(1.0+gamma)

RETURN
END SUBROUTINE weight
