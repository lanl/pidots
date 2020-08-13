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
INTEGER :: i, j, k, m, g, gp, bit, ierr, temp, ieq, its
REAL*8 :: xsct

INCLUDE 'mpif.h'

neq = nx*ny*nz
xys = nx*ny
xzs = nx*nz
yzs = ny*nz
bcs = apo*(xys+xzs+yzs)
bcs2 = xys+xzs+yzs
rn1 = rn + 1
! Allocate the solution vectors/matrices
ALLOCATE(f(neq,ng))
ALLOCATE(e(neq),bs(neq),ba(bcs,8))
ALLOCATE(jmat(neq,neq))
IF (tpose == 1) THEN
   ALLOCATE(kmat(bcs,neq,8), jpsi(neq,bcs,8))
ELSE
   ALLOCATE(kmat(neq,bcs,8), jpsi(bcs,neq,8))
END IF
ALLOCATE(kpsi(bcs2,bcs2,apo,8))
ALLOCATE(psio(bcs,8))

! Mark the beginning of the solution phase
IF (irank == root) THEN
   WRITE (8,*)
   WRITE (8,*) "-------------------------- THE SOLUTION ----------------------------------"
   WRITE (8,*)
END IF

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
                  ieq = i + (j-1)*nx + (k-1)*xys
                  s(i,j,k,g) = s(i,j,k,g) + xsct*f(ieq,gp)
               END DO
            END DO
         END DO
      END DO
   END IF

   CALL idomats(g)
   CALL CPU_TIME(tjmat)

   ! Call for the GMRES solver: 0=MGS, 1=CGS
   IF (om == 1) THEN
      CALL gmrescgs(g)
   ELSE
      CALL gmresmgs(g)
   END IF
   CALL CPU_TIME(tsolve)
END DO

! Deallocate variables no longer needed
DEALLOCATE(e,jmat,kmat,jpsi,kpsi)
IF (pcf == 1) THEN
   DEALLOCATE(prec)
   IF (pinv == 0) DEALLOCATE(piv)
END IF

! Check the scalar flux printing flag
IF (sfp == 1) THEN
   ! Root allocates the primary flux solution holder
   IF (irank == root) ALLOCATE(flux(nxt,nyt,nzt,ng))

   ! Bring all the blocks' solutions for all the groups back to root for output
   CALL idogthr
END IF

RETURN
END SUBROUTINE solve
