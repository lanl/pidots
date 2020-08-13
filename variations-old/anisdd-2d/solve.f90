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
INTEGER :: i, j, t, m, g, gp, neq, bcs, bcs2, bit, ierr, temp, ieq, its, l, ll
REAL*8 :: xsct

INCLUDE 'mpif.h'

! Set the number of scattering moments for allocating dimensions
sord = anord + 1
nmom = sord*(sord+1)/2

neq = nx*ny*nmom
bcs = apo*(nx+ny)
! Allocate the solution vectors/matrices depending on method
IF (meth == 0) THEN
   ALLOCATE(f(nmom,nx,ny,ng), sm(neq,ng))
   ALLOCATE(e(nmom,nx,ny))
   ALLOCATE(fold(nmom,nx,ny))
ELSE
   bcs2 = nx+ny
   ALLOCATE(phi(neq,ng), sm(neq,ng))
   ALLOCATE(phiold(neq), src(neq), sv(neq))
   IF (tpose == 1) THEN
      ALLOCATE(kmat(bcs,neq,4), jpsi(neq,bcs,4))
   ELSE
      ALLOCATE(kmat(neq,bcs,4), jpsi(bcs,neq,4))
   END IF
   ALLOCATE(kpsi(bcs2,bcs2,apo,4))
END IF

! Allocate the outward angular flux value (both SI and ITM use this)
ALLOCATE(psio(bcs,4,ng))

ALLOCATE(cnvf(ng), bcnvf(ng))

! Mark the beginning of the solution phase
IF (irank == root) THEN
   WRITE (8,*)
   WRITE (8,*) "-------------------------- THE SOLUTION ----------------------------------"
   WRITE (8,*)
END IF

! Get the time to reach this point
CALL CPU_TIME(ttosolve)

! Initialize the source moment vector, then place external zero moment source in it
sm = 0.0
DO g = 1, ng
   DO j = 1, ny
      DO i = 1, nx
         ieq = ((j-1)*nx + (i-1))*nmom + 1
         sm(ieq,g) = s(i,j,g)
      END DO
   END DO
END DO

! Start the loop over all energy groups
DO g = 1, ng
   ! Reset the source as external + scattering
   IF (g > 1) THEN ! Downscattering only considered
      DO gp = 1, (g-1)
         DO j = 1, ny
            DO i = 1, nx
               m = mat(i,j)
               DO l = 0, anord
                  xsct = sigs(l,m,g,gp)
                  DO ll = 0, l
                     t = l*(l+1)/2 + ll + 1
                     ieq = ((j-1)*nx + (i-1))*nmom + t
                     IF (meth == 1) THEN
                        sm(ieq,g) = sm(ieq,g) + xsct*phi(ieq,gp)
                     ELSE
                        sm(ieq,g) = sm(ieq,g) + xsct*f(t,i,j,gp)
                     END IF
                  END DO
               END DO
            END DO
         END DO
      END DO
   END IF

   ! Check method, if IDO, need to construct matrices for solving first
   IF (meth == 1) THEN
      CALL idomats(g,neq,bcs)
      CALL CPU_TIME(tjmat)
   END IF

   ! Iterations of the parallel blocks
   DO bit = 1, bitmx

      ! Check which solution scheme will be employed
      IF (meth == 0) THEN
         CALL inner(g,its)
      ELSE IF (meth == 1) THEN
         CALL idot(g,its,neq,bcs,bcs2)
      END IF

      CALL MPI_REDUCE(cnvf(g),temp,1,MPI_INTEGER,MPI_MIN,root,MPI_COMM_WORLD,ierr)
      IF (irank == root) THEN
         IF (temp /= 1) THEN
            warn = warn + 1
            WRITE (8,'(/,1X,A,I2,A,I5,A,/)') "WARNING: Group ", g, " does not have a converged &
                   solution in block iteration ", bit, " for at least one of the blocks."
         END IF
      END IF

      ! Now have a scalar flux and angular flux moments either from SI or ITM
      ! Call for a convergence check and send angular flux moments if necessary
      CALL pbj(g,neq,bit,its)

      IF (bcnvf(g) == 1) EXIT
      IF (bcnvf(g) == 0 .AND. bit == bitmx) THEN
         CALL MPI_FINALIZE(ierr)
      END IF
   END DO

   IF (irank == root) THEN
      IF (meth == 0) WRITE(8,'(1X,A,I4,A)') "Group", g, " source iterations..."
      IF (meth == 1) THEN
         WRITE(8,'(1X,A,I4,A)') "Group", g, " Integral discrete ordinates..."
         IF (idos == 0) WRITE(8,'(3X,A)') "...with direct solution..."
         IF (idos == 1) WRITE(8,'(3X,A)') "...with CG iterations..."
      END IF
   END IF

END DO

! Deallocate variables no longer needed
IF (meth == 0) THEN
   DEALLOCATE(e,fold,sm,cnvf)
ELSE
   DEALLOCATE(phiold,sm,src,sv,kmat,jpsi,kpsi,cnvf)
END IF

! Check the scalar flux printing flag
IF (sfp == 1) THEN
   ! Root allocates the primary flux solution holder
   ! Don't need all moments, will only be using the 00 moment.
   IF (irank == root) ALLOCATE(flux(nxt,nyt,ng))

   ! Bring all the blocks' solutions for all the groups back to root for output
   ! Use different routines for easier reading
   IF (meth == 0) THEN
      CALL sigthr(neq)
   ELSE
      CALL idogthr(neq)
   END IF
END IF

RETURN
END SUBROUTINE solve
