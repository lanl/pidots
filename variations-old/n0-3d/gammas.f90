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
alf = 1.0/(ao*bo*go-ao*bo*ez-ao*go*ey-bo*go*ex)

! Now need 16 formulas to make 16 gamma elements
! gaa
gaa = alf*ao*bo*go

! gax, gaxz, gayz
gaxy = alf*ez*ao*bo*(go-gi)
gaxz = alf*ey*ao*go*(bo-bi)
gayz = alf*ex*bo*go*(ao-ai)

! gxya, gxza, gyza
gxya = -alf*ao*bo
gxza = -alf*ao*go
gyza = -alf*bo*go

! gxyxy, gxyxz, gxyyz, gxzxy, gxzxz, gxzyz, gyzxy, gyzxz, gyzyz
gxyxy = alf*(-ao*bo*ez+gi*(ao*bo-bo*ex-ao*ey))
gxyxz = alf*ey*ao*(bi-bo)
gxyyz = alf*ex*bo*(ai-ao)

gxzxy = alf*ez*ao*(gi-go)
gxzxz = alf*(-ao*go*ey+bi*(ao*go-go*ex-ao*ez))
gxzyz = alf*ex*go*(ai-ao)

gyzxy = alf*ez*bo*(gi-go)
gyzxz = alf*ey*go*(bi-bo)
gyzyz = alf*(-bo*go*ex+ai*(bo*go-go*ey-bo*ez))

RETURN
END SUBROUTINE gammas
