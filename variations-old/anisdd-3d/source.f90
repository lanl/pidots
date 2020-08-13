SUBROUTINE source(g)

!--------------------------------------------------------
!
! Set the source up from the read-in external 00 source
!  and the downscattering events
!
!--------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER, INTENT(IN) :: g
INTEGER :: gp, i, j, k, m, l, ll, t1, t2, ieq
REAL*8 :: xsct

! Initialize sm(:,g)
sm(:,g) = 0.0

! Place external 00 moment source in sm
DO k = 1, nz
   DO j = 1, ny
      DO i = 1, nx
         ieq = ((k-1)*ny*nx + (j-1)*nx + (i-1))*nmom + 1
         sm(ieq,g) = s(i,j,k,g)
      END DO
   END DO
END DO

! Reset the source as external + scattering
IF (g > 1) THEN ! Downscattering only considered
   DO gp = 1, (g-1)
      DO k = 1, nz
         DO j = 1, ny
            DO i = 1, nx
               m = mat(i,j,k)
               t1 = ((k-1)*ny*nx + (j-1)*nx + (i-1))*nmom
               DO l = 0, anord
                  xsct = sigs(l,m,g,gp)
                  ! Update the l,0th moment
                  t2 = l**2 + 1
                  ieq = t1 + t2
                  IF (meth == 1) THEN
                     sm(ieq,g) = sm(ieq,g) + xsct*phi(ieq,g)
                  ELSE
                     sm(ieq,g) = sm(ieq,g) + xsct*f(t2,i,j,k,gp)
                  END IF
                  ! Now update the ll>0 moments
                  DO ll = 1, l
                     t2 = l**2 + 2*ll
                     ieq = t1 + t2
                     IF (meth == 1) THEN
                        sm(ieq,g) = sm(ieq,g) + xsct*phi(ieq,g)            ! Even function
                        sm(ieq+1,g) = sm(ieq+1,g) + xsct*phi(ieq+1,g)      ! Odd function
                     ELSE
                        sm(ieq,g) = sm(ieq,g) + xsct*f(t2,i,j,k,gp)        ! Even function
                        sm(ieq+1,g) = sm(ieq+1,g) + xsct*f(t2+1,i,j,k,gp)  ! Odd function
                     END IF
                  END DO
               END DO
            END DO
         END DO
      END DO
   END DO
END IF

RETURN
END SUBROUTINE source
