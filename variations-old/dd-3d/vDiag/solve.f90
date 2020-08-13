SUBROUTINE solve

!-------------------------------------------------------------
!
!  Directs the solution by either calling for a mesh
!   sweep (DD-SI) or for a ITM solution (DD-ITM)
! 
!-------------------------------------------------------------

USE invar
USE solvar
USE timevar
IMPLICIT NONE
INTEGER :: i, j, k, m, g, gp, ierr, temp, ieq, its
REAL*8 :: xsct

ALLOCATE(phi(nx,ny,nz,ng),src(nx,ny,nz),sv(nx,ny,nz))
ALLOCATE(cnvf(ng))

! Mark the beginning of the solution phase
WRITE (8,*)
WRITE (8,*) "-------------------------- THE SOLUTION ----------------------------------"
WRITE (8,*)

! Get the time to reach this point
CALL CPU_TIME(ttosolve)

! Start the loop over all energy groups
DO g = 1, ng
   ! Reset the source as external + scattering
   IF (g > 1) THEN ! Downscattering only considered
      DO gp = 1, (g-1)
         DO k = 1, nz
            DO j = 1, ny
               DO i = 1, nx
                  m = mat(i,j,k)
                  xsct = sigs(m,g,gp)
                  s(i,j,k,g) = s(i,j,k,g) + xsct*phi(i,j,k,gp)
               END DO
            END DO
         END DO
      END DO
   END IF

   ! Check method, if IDO, need to construct matrices for solving first
   CALL idomats(g)
   CALL CPU_TIME(tjmat)

   CALL idot(g,its)

   CALL CPU_TIME(tsolve)

   IF (cnvf(g) /= 1) THEN
      warn = warn + 1
      WRITE (8,'(/,1X,A,I2,A,I5,A,/)') "WARNING: Group ", g, " did not converge."
   END IF

   WRITE(8,'(1X,A,I4,A)') "Group", g, " Integral discrete ordinates..."

END DO

! Deallocate variables no longer needed
DEALLOCATE(src,sv,cnvf,jmat)

RETURN
END SUBROUTINE solve
