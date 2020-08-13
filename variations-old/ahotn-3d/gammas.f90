SUBROUTINE gammas

!-------------------------------------------------------------
!
! Construct the gamma sub-matrices from gmat
! 
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER :: i, j

! Now make the submatrices from gmat
! gaa
DO j = 1, ordcb
   DO i = 1, ordcb
      gaa(i,j) = gmat(i,j)
   END DO
END DO
! gaxy, gaxz, & gayz
DO j = 1, ordsq
   DO i = 1, ordcb
      gaxy(i,j) = gmat(i,(ordcb+j))
      gaxz(i,j) = gmat(i,(ordcb+ordsq+j))
      gayz(i,j) = gmat(i,(ordcb+2*ordsq+j))
   END DO
END DO

! gxya, gxza, & gyza
DO j = 1, ordcb
   DO i = 1, ordsq
      gxya(i,j) = gmat((ordcb+i),j)
      gxza(i,j) = gmat((ordcb+ordsq+i),j)
      gyza(i,j) = gmat((ordcb+2*ordsq+i),j)
   END DO
END DO

! gxyxy, gxyxz, gxyyz, gxzxy, gxzxz, gxzyz, gyzxy, gyzxz, gyzyz
DO j = 1, ordsq
   DO i = 1, ordsq
      gxyxy(i,j) = gmat((ordcb+i),(ordcb+j))
      gxyxz(i,j) = gmat((ordcb+i),(ordcb+ordsq+j))
      gxyyz(i,j) = gmat((ordcb+i),(ordcb+2*ordsq+j))
      gxzxy(i,j) = gmat((ordcb+ordsq+i),(ordcb+j))
      gxzxz(i,j) = gmat((ordcb+ordsq+i),(ordcb+ordsq+j))
      gxzyz(i,j) = gmat((ordcb+ordsq+i),(ordcb+2*ordsq+j))
      gyzxy(i,j) = gmat((ordcb+2*ordsq+i),(ordcb+j))
      gyzxz(i,j) = gmat((ordcb+2*ordsq+i),(ordcb+ordsq+j))
      gyzyz(i,j) = gmat((ordcb+2*ordsq+i),(ordcb+2*ordsq+j))
   END DO
END DO

RETURN
END SUBROUTINE gammas
