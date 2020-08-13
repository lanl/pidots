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
CHARACTER(8) :: fl
INTEGER :: i, j, k, m, g, gp, bit, ierr, temp, ieq, its, o, n, jeq
REAL*8 :: xsct, tmp
REAL*8, DIMENSION(:), ALLOCATABLE :: rphi, rsv
REAL*8, DIMENSION(:,:), ALLOCATABLE :: rpsi, rpave

INCLUDE 'mpif.h'

neq = nx*ny*nz
xys = nx*ny
xzs = nx*nz
yzs = ny*nz
bcs = apo*(xys+xzs+yzs)
! Allocate the solution vectors/matrices depending on method
IF (meth == 0) THEN
   ALLOCATE(f(nx,ny,nz,ng))
   ALLOCATE(fold(nx,ny,nz))
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
                  IF (meth == 0) THEN
                     s(i,j,k,g) = s(i,j,k,g) + xsct*f(i,j,k,gp)
                  ELSE
                     ieq = i + (j-1)*nx + (k-1)*xys
                     s(i,j,k,g) = s(i,j,k,g) + xsct*phi(ieq,gp)
                  END IF
               END DO
            END DO
         END DO
      END DO
   END IF

   ! Check method, if IDO, need to construct matrices for solving first
   IF (meth == 1) THEN
      CALL idomats(g)
      phiold = 0.0
      CALL CPU_TIME(tjmat)
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
      CALL pbj(g,bit,its)

      IF (bcnvf(g) == 1) EXIT
      IF (bcnvf(g) == 0 .AND. bit == bitmx) THEN
!         CALL MPI_FINALIZE(ierr)
!         STOP
         EXIT
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

! Check the scalar flux printing flag
!IF (sfp == 1) THEN
   ! Root allocates the primary flux solution holder
   IF (irank == root) ALLOCATE(flux(nxt,nyt,nzt,ng))
   IF (irank == root) THEN
      ALLOCATE(afx(nyt,nzt,apo,8),afy(nxt,nzt,apo,8),afz(nxt,nyt,apo,8))
      afz = 0.0
      afy = 0.0
      afx = 0.0
   END IF

   ! Bring all the blocks' solutions for all the groups back to root for output
   ! Use different routines for easier reading
   IF (meth == 0) THEN
      CALL sigthr
   ELSE
      CALL idogthr
      CALL psigthr
   END IF

! Print flux to special file
IF (irank == root) THEN
   OPEN (UNIT=20, FILE="finephi")
   DO k = 1, nzt
      DO j = 1, nyt
         DO i = 1, nxt
            WRITE(20,*) flux(i,j,k,1)
         END DO
      END DO
   END DO
! z
   DO o = 1, 8
      DO n = 1, apo
         DO j = 1, nyt
            DO i = 1, nxt
               WRITE(20,*) afz(i,j,n,o)
            END DO
         END DO
      END DO
   END DO
! y
   DO o = 1, 8
      DO n = 1, apo
         DO k = 1, nzt
            DO i = 1, nxt
               WRITE(20,*) afy(i,k,n,o)
            END DO
         END DO
      END DO
   END DO
! x
   DO o = 1, 8
      DO n = 1, apo
         DO k = 1, nzt
            DO j = 1, nyt
               WRITE(20,*) afx(j,k,n,o)
            END DO
         END DO
      END DO
   END DO
END IF

!END IF

! Compute the residual -- first all fine grid contributions:
ALLOCATE(rphi(neq),rsv(neq))
ALLOCATE(rpsi(bcs,8))
rphi = 0.0
rpsi = 0.0
CALL residual(1,rphi,rpsi)

! Restrict
! phi
! Get average of rphi
rsv = 0.0
tmp = SUM(rphi)/neq
CALL MPI_GATHER(tmp,1,MPI_DOUBLE_PRECISION,rsv,1,MPI_DOUBLE_PRECISION,root,allcomm,ierr)
! psi
ALLOCATE(rpave(3*apo,8))
rpave = 0.0
DO o = 1, 8
   DO n = 1, apo
      ieq = n
      DO j = 1, ny
         DO i = 1, nx
            jeq = (n-1)*xys + (j-1)*nx + i
            rpave(ieq,o) = rpave(ieq,o) + rpsi(jeq,o)
         END DO
      END DO
      rpave(ieq,o) = rpave(ieq,o)/xys
   END DO
   DO n = 1, apo
      ieq = apo + n
      DO k = 1, nz
         DO i = 1, nx
            jeq = apo*xys + (n-1)*xzs + (k-1)*nx + i
            rpave(ieq,o) = rpave(ieq,o) + rpsi(jeq,o)
         END DO
      END DO
      rpave(ieq,o) = rpave(ieq,o)/xzs
   END DO
   DO n = 1, apo
      ieq = 2*apo + n
      DO k = 1, nz
         DO j = 1, ny
            jeq = apo*(xys+xzs) + (n-1)*yzs + (k-1)*ny + j
            rpave(ieq,o) = rpave(ieq,o) + rpsi(jeq,o)
         END DO
      END DO
      rpave(ieq,o) = rpave(ieq,o)/yzs
   END DO
END DO
CALL rpsigthr(rpsi,rpave)

! Print out to file
IF (irank == root) THEN
   OPEN (UNIT=21, FILE="residuals")
   DO i = 1, neq
      WRITE(21,*) rsv(i)
   END DO

   DO j = 1, 8
      DO i = 1, bcs
         WRITE(21,*) rpsi(i,j)
      END DO
   END DO
END IF

WRITE(fl,'(''psi.'',I1.1)') nrank
OPEN(UNIT=25, FILE=TRIM(fl))
DO j = 1, 8
   DO i = 1, bcs
      WRITE(25,*) psii(i,j)
   END DO
END DO







! Deallocate variables no longer needed
IF (meth == 0) THEN
   DEALLOCATE(fold,cnvf)
ELSE
   DEALLOCATE(phiold,sv,src,kmat,jpsi,kpsi,cnvf)
   IF (idos == 0 .AND. invf == 0) THEN
      DEALLOCATE(jmat,jmat2)
      IF (sym == 0) DEALLOCATE(piv)
   END IF
   IF (idos == 1) THEN
      DEALLOCATE(jmat)
      IF (pcf == 1) DEALLOCATE(dv)
   END IF
END IF

RETURN
END SUBROUTINE solve
