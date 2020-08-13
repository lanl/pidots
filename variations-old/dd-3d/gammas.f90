SUBROUTINE gammas

!-------------------------------------------------------------
!
! Compute gamma elements from known analytic form of gamma
! 
!-------------------------------------------------------------

USE solvar
IMPLICIT NONE
REAL*8 :: alf

! Set the scalar alpha
alf = 1.0/(-0.125-0.25*ez-0.25*ey-0.25*ex)

! Now need 16 formulas to make 16 gamma elements
! gaa
gaa = alf*(-0.125)

! gax, gaxz, gayz
gaxy = alf*(-0.25)*ez
gaxz = alf*(-0.25)*ey
gayz = alf*(-0.25)*ex

! gxya, gxza, gyza
gxya = alf*(-0.25) 
gxza = gxya
gyza = gxya

! gxyxy, gxyxz, gxyyz, gxzxy, gxzxz, gxzyz, gyzxy, gyzxz, gyzyz
gxyxy = alf*(0.5*(0.25+0.5*(ey+ex))-0.25*ez)
gxyxz = 2.0*gaxz
gxyyz = 2.0*gayz

gxzxy = 2.0*gaxy
gxzxz = alf*(0.5*(0.25+0.5*(ez+ex))-0.25*ey)
gxzyz = 2.0*gayz

gyzxy = 2.0*gaxy
gyzxz = 2.0*gaxz
gyzyz = alf*(0.5*(0.25+0.5*(ez+ey))-0.25*ex)

RETURN
END SUBROUTINE gammas
