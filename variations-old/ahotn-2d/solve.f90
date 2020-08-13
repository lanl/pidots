SUBROUTINE solve(root)

!-------------------------------------------------------------
!
!  Directs the solution by either calling for a mesh
!   sweep (AHOT-N) or for a ITM solution (AHOT-N-NS)
! 
!-------------------------------------------------------------

USE invar
USE solvar
USE timevar
IMPLICIT NONE
INTEGER, INTENT(IN) :: root
INTEGER :: i, j, k, l, m, g, gp, neq, bcs, bcs2, bit, temp
INTEGER :: pi, pj, root2, eqs, ieq, tag, indx, jndx, trank, ierr
INTEGER, DIMENSION(2) :: coord
INTEGER, DIMENSION(3) :: istat
REAL*8 :: xsct
REAL*8, DIMENSION(:,:), ALLOCATABLE :: tmphi
REAL*8, DIMENSION(:,:,:,:,:), ALLOCATABLE :: tmpf

INCLUDE 'mpif.h'

! Root allocates the primary flux solution holder
IF (irank == root) ALLOCATE(flux(nxt,nyt,0:lambda,0:lambda,ng))

neq = nx*ny*ordsq
bcs = apo*(nx+ny)*order
! Allocate the solution vectors/matrices depending on method
IF (meth == 0) THEN
   ALLOCATE(f(nx,ny,0:lambda,0:lambda,ng))
   ALLOCATE(e(nx,ny,0:lambda,0:lambda))
   ALLOCATE(fold(nx,ny,0:lambda,0:lambda))
ELSE
   bcs2 = (nx+ny)*order
   ALLOCATE(phi(neq,ng))
   ALLOCATE(phiold(neq), src(neq), sv(neq))
   ALLOCATE(kmat(neq,bcs,4), jpsi(bcs,neq,4), kpsi(bcs2,bcs2,apo,4))   
END IF

! Allocate the outward angular flux value (both SI and ITM use this)
ALLOCATE(psio(bcs,4,ng))

ALLOCATE(cnvf(ng), bcnvf(ng))

! Intitialize warn to indicate where warnings may occur
warn = 0

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
         DO j = 1, ny
            DO i = 1, nx
               m = mat(i,j)
               xsct = sigs(m,g,gp)
               DO l = 0, lambda
                  DO k = 0, lambda
                     IF (meth == 0) THEN
                        s(i,j,k,l,g) = s(i,j,k,l,g) + xsct*f(i,j,k,l,gp)
                     ELSE
                        ieq = ((i-1) + (j-1)*nx)*ordsq + k*order + l + 1
                        s(i,j,k,l,g) = s(i,j,k,l,g) + xsct*phi(ieq,gp)
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
  
   IF (meth == 0) THEN
      IF (irank == root) WRITE(8,'(1X,A,I4,A)') "Group", g, " source iterations..."
   ELSE
      IF (irank == root) WRITE(8,'(1X,A,I4,A)') "Group", g, " Integral discrete ordinates..."
   END IF

   ! Iterations of the parallel blocks
   DO bit = 1, bitmx

      ! Check which solution scheme will be employed
      IF (meth == 0) THEN
         ! Call for the inner iteration: SI
         CALL inner(g)
      ELSE IF (meth == 1) THEN
         CALL idot(g)
      END IF
      CALL MPI_REDUCE(cnvf(g),temp,1,MPI_INTEGER,MPI_MIN,root,MPI_COMM_WORLD,ierr)
      IF (irank == root) THEN
         IF (temp /= 1) THEN
            warn = warn + 1
            WRITE (8,'(/,1X,A,I2,A,I5,A,/)') "WARNING: Group ", g, " does not have a converged &
                   solution in block iteration ", bit, " for at least one of the blocks"
         END IF
      END IF
      ! Now have a scalar flux and angular flux moments either from SI or ITM
      ! Call for a convergence check and send angular flux moments if necessary
      CALL pbj(g,bit,root)

      IF (bcnvf(g) == 1) EXIT
      IF (bcnvf(g) == 0 .AND. bit == bitmx) THEN
         CALL MPI_FINALIZE(ierr)
         STOP
      END IF
   END DO
END DO

! Get the time out of the solution
CALL CPU_TIME(tsolve)

! Setup for parallel
IF (irank == root) root2 = nrank
CALL MPI_BCAST(root2,1,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)

! Bring all the blocks' solutions for all the groups back to root for output
IF (meth == 0) THEN
   IF (nrank /= root2) THEN
      tag = 100 + nrank
      CALL MPI_SEND(f,neq*ng,MPI_DOUBLE_PRECISION,root2,tag,allcomm,ierr)
   ELSE
      jndx = 0
      DO pj = 0, npy-1
         coord(2) = pj
         indx = 0
         DO pi = 0, npx-1
            coord(1) = pi
            eqs = nxvec(pi)*nyvec(pj)*ordsq*ng
            CALL MPI_CART_RANK(allcomm,coord,trank,ierr)
            IF (trank /= nrank) THEN
               ALLOCATE(tmpf(nxvec(pi),nyvec(pj),0:lambda,0:lambda,ng))
               tmpf = 0.0
               tag = 100 + trank
               CALL MPI_RECV(tmpf,eqs,MPI_DOUBLE_PRECISION,trank,tag,allcomm,istat,ierr)
               DO g = 1, ng
                  DO k = 0, lambda
                     DO l = 0, lambda
                        DO j = 1, nyvec(pj)
                           DO i = 1, nxvec(pi)
                              flux(i+indx,j+jndx,k,l,g) = tmpf(i,j,k,l,g)
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
               DEALLOCATE(tmpf)
            ELSE IF (trank == nrank) THEN
               ! Root copies its own solution into flux
               DO g = 1, ng
                  DO k = 0, lambda
                     DO l = 0, lambda
                        DO j = 1, ny
                           DO i = 1, nx
                              flux(i+indx,j+jndx,k,l,g) = f(i,j,k,l,g)
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END IF
            indx = indx + nxvec(pi)
         END DO
         jndx = jndx + nyvec(pj)
      END DO
   END IF
ELSE
   IF (nrank /= root2) THEN
      tag = 100 + nrank
      CALL MPI_SEND(phi,neq*ng,MPI_DOUBLE_PRECISION,root2,tag,allcomm,ierr)
   ELSE
      jndx = 0
      DO pj = 0, npy-1
         coord(2) = pj
         indx = 0
         DO pi = 0, npx-1
            coord(1) = pi
            eqs = nxvec(pi)*nyvec(pj)*ordsq
            CALL MPI_CART_RANK(allcomm,coord,trank,ierr)
            IF (trank /= nrank) THEN
               ALLOCATE(tmphi(eqs,ng))
               tmphi = 0.0
               tag = 100 + trank
               CALL MPI_RECV(tmphi,eqs*ng,MPI_DOUBLE_PRECISION,trank,tag,allcomm,istat,ierr)
               DO g = 1, ng
                  DO k = 0, lambda
                     DO l = 0, lambda
                        DO j = 1, nyvec(pj)
                           DO i = 1, nxvec(pi)
                              ieq = ((i-1) + (j-1)*nxvec(pi))*ordsq + k*order + l + 1
                              flux(i+indx,j+jndx,k,l,g) = tmphi(ieq,g)
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
               DEALLOCATE(tmphi)
            ELSE IF (trank == nrank) THEN
               ! Root copies its own solution into flux
               DO g = 1, ng
                  DO k = 0, lambda
                     DO l = 0, lambda
                        DO j = 1, ny
                           DO i = 1, nx
                              ieq = ((i-1) + (j-1)*nx)*ordsq + k*order + l + 1
                              flux(i+indx,j+jndx,k,l,g) = phi(ieq,g)
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END IF
            indx = indx + nxvec(pi)
         END DO
         jndx = jndx + nyvec(pj)
      END DO
   END IF
END IF

IF (meth == 0) THEN
   DEALLOCATE(f,e,fold,cnvf)
ELSE
   DEALLOCATE(phi,phiold,src,sv,kmat,jpsi,kpsi,cnvf)
END IF

RETURN
END SUBROUTINE solve
