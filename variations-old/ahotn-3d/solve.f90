SUBROUTINE solve(ttosolve,tjmat,tsolve)

!-------------------------------------------------------------
!
!  Directs the solution by either calling for a mesh
!   sweep (DD-SI) or for a ITM solution (DD-ITM)
! 
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER :: i, j, k, t, u, v, m, g, gp, bit, ierr, temp, ieq, its
REAL*8 :: xsct
REAL*8, INTENT(OUT) :: ttosolve, tjmat, tsolve

INCLUDE 'mpif.h'

neq = nx*ny*nz*ordcb
xys = nx*ny*ordsq
xzs = nx*nz*ordsq
yzs = ny*nz*ordsq
bcs = apo*(xys+xzs+yzs)
! Allocate the solution vectors/matrices depending on method
IF (meth == 0) THEN
   ALLOCATE(f(nx,ny,nz,0:lambda,0:lambda,0:lambda,ng))
   ALLOCATE(fold(nx,ny,nz,0:lambda,0:lambda,0:lambda))
ELSE
   bcs2 = xys+xzs+yzs
   ALLOCATE(phi(neq,ng))
   ALLOCATE(phiold(neq),src(neq),sv(neq))
   IF (tpose == 1) THEN
      ALLOCATE(kmat(bcs,neq,8),jpsi(neq,bcs,8))
   ELSE
      ALLOCATE(kmat(neq,bcs,8),jpsi(bcs,neq,8))
   END IF
   ALLOCATE(kpsi(bcs2,bcs2,apo,8))
END IF

! Allocate the outward angular flux value (both SI and ITM use this)
ALLOCATE(psio(bcs,8,ng))

ALLOCATE(cnvf(ng), bcnvf(ng))

! Mark the beginning of the solution phase
IF (irank == root) THEN
   WRITE (8,*)
   WRITE (8,*) "-------------------------- THE SOLUTION ----------------------------------"
   WRITE (8,*)
END IF

! Get the time to reach this point
ttosolve = MPI_WTIME()

DO g = 1, ng
   ! Reset the source as external + scattering
   IF (g > 1) THEN ! Downscattering only considered
      DO gp = 1, (g-1)
         DO k = 1, nz
            DO j = 1, ny
               DO i = 1, nx
                  m = mat(i,j,k)
                  xsct = sigs(m,g,gp)
                  DO v = 0, lambda
                     DO u = 0, lambda
                        DO t = 0, lambda
                           IF (meth == 0) THEN
                              s(i,j,k,t,u,v,g) = s(i,j,k,t,u,v,g) + xsct*f(i,j,k,t,u,v,gp)
                           ELSE
                              ieq = ((i-1) + (j-1)*nx + (k-1)*nx*ny)*ordcb + t*ordsq + u*order + v + 1
                              s(i,j,k,t,u,v,g) = s(i,j,k,t,u,v,g) + xsct*phi(ieq,gp)
                           END IF
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
   END IF

   ! Check method, if IDO, need to construct matrices for solving first
   IF (meth == 1) THEN
      CALL idomats(g)
      phiold = 0.0
      tjmat = MPI_WTIME()
   END IF

   ! Iterations of the parallel blocks
   DO bit = 1, bitmx

      ! Check which solution scheme will be employed
      IF (meth == 0) THEN
         CALL inner(g,its)
      ELSE IF (meth == 1) THEN
         CALL idot(g,its,bit)
      END IF

!      CALL MPI_REDUCE(cnvf(g),temp,1,MPI_INTEGER,MPI_MIN,root,MPI_COMM_WORLD,ierr)
!      IF (irank == root) THEN
!         IF (temp /= 1) THEN
!            warn = warn + 1
!            WRITE (8,'(/,1X,A,I2,A,I5,A,/)') "WARNING: Group ", g, " does not have a converged &
!                   solution in block iteration ", bit, " for at least one of the blocks."
!         END IF
!      END IF

      ! Now have a scalar flux and angular flux moments either from SI or ITM
      ! Call for a convergence check and send angular flux moments if necessary
      CALL pbj(g,bit,its,tsolve)

      IF (bcnvf(g) == 1) EXIT
      IF (bcnvf(g) == 0 .AND. bit == bitmx) THEN
         CALL MPI_FINALIZE(ierr)
         STOP
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
   DEALLOCATE(fold,cnvf)
ELSE
   DEALLOCATE(phiold,src,sv,kmat,jpsi,kpsi,cnvf)
   IF (idos == 0 .AND. invf == 0) THEN
      DEALLOCATE(jmat)
      IF (sym == 0) DEALLOCATE(piv)
   END IF
   IF (idos == 1) THEN
      DEALLOCATE(jmat)
      IF (pcf == 1) DEALLOCATE(dv)
   END IF
END IF

! Check the flux printing flag
IF (sfp == 1) THEN
   ! Root allocates the primary flux solution holder
   IF (irank == root) ALLOCATE(flux(nxt,nyt,nzt,0:lambda,0:lambda,0:lambda,ng))

   ! Bring all the blocks' solutions for all the groups back to root for output
   ! Use different routines for easier reading
   IF (meth == 0) THEN
      CALL sigthr
   ELSE
      CALL idogthr
   END IF
END IF

RETURN
END SUBROUTINE solve
