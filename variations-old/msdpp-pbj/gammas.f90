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
alf = 1.0d0/(-0.125d0-0.25d0*ez-0.25d0*ey-0.25d0*ex)

! Now need 16 formulas to make 16 gamma elements
! gaa
gaa = alf*(-0.125d0)

! gax, gaxz, gayz
gaxy = alf*(-0.25d0)*ez
gaxz = alf*(-0.25d0)*ey
gayz = alf*(-0.25d0)*ex

! gxya, gxza, gyza
gxya = alf*(-0.25d0)
gxza = gxya
gyza = gxya

! gxyxy, gxyxz, gxyyz, gxzxy, gxzxz, gxzyz, gyzxy, gyzxz, gyzyz
gxyxy = alf*(0.5d0*(0.25d0+0.5d0*(ey+ex))-0.25d0*ez)
gxyxz = 2.0d0*gaxz
gxyyz = 2.0d0*gayz

gxzxy = 2.0d0*gaxy
gxzxz = alf*(0.5d0*(0.25d0+0.5d0*(ez+ex))-0.25d0*ey)
gxzyz = 2.0d0*gayz

gyzxy = 2.0d0*gaxy
gyzxz = 2.0d0*gaxz
gyzyz = alf*(0.5d0*(0.25d0+0.5d0*(ez+ey))-0.25d0*ex)

RETURN
END SUBROUTINE gammas
