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
INTEGER :: i, l, ieq, m, g, gp, bit, ierr, temp, its, nsm
REAL*8 :: xsct

INCLUDE 'mpif.h'

! Set the number of scattering moments for allocating dimensions
sord = anord + 1
nsm = nx*(sord)

! Allocate the solution vectors/matrices depending on method
IF (meth == 0) THEN
   ALLOCATE(f(nsm,ng), sm(nsm,ng))
   ALLOCATE(fold(nsm))
ELSE
   ALLOCATE(f(nsm,ng), sm(nsm,ng))
   ALLOCATE(fold(nsm),src(nsm),sv(nsm))
   IF (tpose == 1) THEN
      ALLOCATE(kmat(apo,nsm,2), jpsi(nsm,apo,2))
   ELSE
      ALLOCATE(kmat(nsm,apo,2), jpsi(apo,nsm,2))
   END IF
   ALLOCATE(kpsi(apo,2))
END IF

! Allocate the outward angular flux value (both SI and ITM use this)
ALLOCATE(psio(apo,2,ng))

ALLOCATE(cnvf(ng), bcnvf(ng))

! Mark the beginning of the solution phase
IF (irank == root) THEN
   WRITE (8,*)
   WRITE (8,*) "-------------------------- THE SOLUTION ----------------------------------"
   WRITE (8,*)
END IF

! Get the time to reach this point
CALL CPU_TIME(ttosolve)

! Initialize source moment vector, then place external zero moment source in it
sm = 0.0
DO g = 1, ng
   DO i = 1, nx
      ieq = (i-1)*sord + 1
      sm(ieq,g) = s(i,g)
   END DO
END DO

! Start the loop over all energy groups
DO g = 1, ng
   ! Reset the source for scattering (isotropic and anisotropic)
   IF (g > 1) THEN ! Downscattering only considered, no fission
      DO gp = 1, (g-1)
         DO i = 1, nx
            m = mat(i)
            DO l = 0, anord
               xsct = sigs(m,l,g,gp)
               ieq = (i-1)*sord + 1 + l
               sm(ieq,g) = sm(ieq,g) + xsct*f(ieq,gp)
            END DO
         END DO
      END DO
   END IF

   ! Set fold
   fold = 0.0
   ! Check method, if IDO, need to construct matrices for solving first
   IF (meth == 1) THEN
      CALL idomats(g,nsm)
      CALL CPU_TIME(tjmat)
   END IF

   ! Iterations of the parallel blocks
   DO bit = 1, bitmx

      ! Check which solution scheme will be employed
      IF (meth == 0) THEN
         CALL inner(g,its,nsm)
      ELSE IF (meth == 1) THEN
         CALL idot(g,nsm,its)
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
      CALL pbj(g,nsm,bit,its)

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

! End the loop over all groups
END DO

! Deallocate variables no longer needed
IF (meth == 0) THEN
   DEALLOCATE(fold,sm,cnvf)
ELSE
   DEALLOCATE(fold,sm,src,sv,kmat,jpsi,kpsi,cnvf)
   IF (idos == 1) DEALLOCATE(jmat)
END IF

! Check the scalar flux printing flag
IF (sfp == 1) THEN
   ! Root allocates the primary flux solution holder
   IF (irank == root) ALLOCATE(flux(nxt*sord,ng))

   ! Bring all the blocks' solutions for all the groups back to root for output
   CALL gthr
END IF

RETURN
END SUBROUTINE solve
