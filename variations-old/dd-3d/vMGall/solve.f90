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
INTEGER :: i, j, k, m, bit, ierr, temp, ieq, its, o, n, jeq, vit, bcnvf, mgcnvf, cnt
REAL*8 :: xsct, phie
REAL*8, DIMENSION(:), ALLOCATABLE :: rphi, rsv
REAL*8, DIMENSION(:,:), ALLOCATABLE :: rpsi, rpave

INCLUDE 'mpif.h'

neq = nx*ny*nz
xys = nx*ny
xzs = nx*nz
yzs = ny*nz
bcs = apo*(xys+xzs+yzs)
bcs2 = xys+xzs+yzs

! Allocate ITMM operators
ALLOCATE(phi(neq),phiold(neq),src(neq),sv(neq),phir(neq))
IF (tpose == 1) THEN
   ALLOCATE(kmat(bcs,neq,8),jpsi(neq,bcs,8))
ELSE
   ALLOCATE(kmat(neq,bcs,8),jpsi(bcs,neq,8))
END IF
ALLOCATE(kpsi(bcs2,bcs2,apo,8))


cneq = npx*npy*npz
cxys = npx*npy
cxzs = npx*npz
cyzs = npy*npz
cbcs = apo*(cxys+cxzs+cyzs)
cbcs2 = cxys+cxzs+cyzs

! Coarse grid operators
IF (irank == root) THEN
   ALLOCATE(f(cneq),fold(cneq))
   IF (tpose == 1) THEN
      ALLOCATE(ckmat(cbcs,cneq,8),cjpsi(cneq,cbcs,8))
   ELSE
      ALLOCATE(ckmat(cneq,cbcs,8),cjpsi(cbcs,cneq,8))
   END IF
   ALLOCATE(ckpsi(cbcs2,cbcs2,apo,8))

   ALLOCATE(flux(nxt,nyt,nzt))

   ALLOCATE(cpsii(cbcs,8),cpsio(cbcs,8))
END IF

! Allocate the outward angular flux value (both SI and ITM use this)
ALLOCATE(psio(bcs,8))
! Allocate the copy of the old inward flux for use in computing residuals
ALLOCATE(oldpsi(bcs,8))

! Allocate and initialize residual vectors
ALLOCATE(rphi(neq),rsv(cneq))
ALLOCATE(rpsi(bcs,8),rpave(3*apo,8))
IF (irank == root) ALLOCATE(rrpsi(cbcs,8))
rphi = 0.0
rpsi = 0.0


! Mark the beginning of the solution phase
IF (irank == root) THEN
   WRITE (8,*)
   WRITE (8,*) "-------------------------- THE SOLUTION ----------------------------------"
   WRITE (8,*)
END IF

! Get the time to reach this point
CALL CPU_TIME(ttosolve)

! Need to construct matrices for solving first
CALL idomats
phiold = 0.0
! Call for the coarse idomats
IF (irank == root) CALL cidomats
CALL CPU_TIME(tjmat)

cnt = 0

! Start the multigrid iterations
DO vit = 1, vitmx


!if (vit == 2) then
!   print *, cnt
!   exit
!end if



   !------------------------------------------------------------------------------------------------
   ! Do the fine grid iterations of the parallel blocks
   DO bit = 1, bitmx

      cnt = cnt + 1

      ! Call to compute phi and psio
      CALL idot(its,bit)

      !CALL MPI_REDUCE(cnvf,temp,1,MPI_INTEGER,MPI_MIN,root,MPI_COMM_WORLD,ierr)
      !IF (irank == root) THEN
      !   IF (temp /= 1) THEN
      !      warn = warn + 1
      !      WRITE (8,'(/,1X,A,I2,A,I5,A,/)') "WARNING: Group ", g, " does not have a converged &
      !             solution in block iteration ", bit, " for at least one of the blocks."
      !   END IF
      !END IF

      ! Now have a scalar flux and angular flux moments either from SI or ITM
      ! Call for a convergence check and send angular flux moments if necessary
      CALL pbj(bit,its,bcnvf)

      IF (bcnvf == 1) THEN
         mgcnvf = 1
         EXIT
      END IF
      IF (bcnvf == 0 .AND. bit == bitmx) THEN
         mgcnvf = 0
         EXIT
      END IF
   END DO

!mgcnvf = 1
!print *, cnt

   IF (mgcnvf == 1) THEN
      IF (irank == root) WRITE(8,'(1X,A,I5)') "A) Multigrid converged in this many v-cycles: ", vit
      IF (irank == root) WRITE(8,'(1X,A,I5)') "A) Total iterations: ", cnt
      CALL CPU_TIME(tsolve)
      EXIT
   END IF

   ! Bring all the blocks' solutions for all the groups back to root for output
   CALL idogthr
   CALL psigthr

   ! Call for the residual computation
   CALL residual(rphi,rpsi)

   ! Call to restrict
   CALL restrict(rphi,rsv,rpsi,rpave)

   !--------------------------------------------------------------------------------------------------
   ! Now want root to solve the coarse grid problem (VAC ONLY FOR NOW)
   IF (irank == root) THEN

      cnt = cnt + 1

      cpsii = 0.0
      CALL cidot(rsv)
   END IF

   !--------------------------------------------------------------------------------------------------
   ! Return to fine grid
   ! Now have the error and can scatter it back
   CALL MPI_SCATTER(f,1,MPI_DOUBLE_PRECISION,phie,1,MPI_DOUBLE_PRECISION,root,allcomm,ierr)

   ! Add error onto phi and perform more PBJ iterations
   phi = phi + phie*phir

!if (nrank == 0) print*, rsv

!if (nrank == 0) print *, f
!print *, irank, phie, phie*phir(1)

   bcnvf = 0
   DO bit = 1, bitmx2

      cnt = cnt + 1

      CALL igidot(its,bit)
      CALL pbj(bit,its,bcnvf)
      IF (bcnvf == 1) THEN
         mgcnvf = 1
         EXIT
      END IF
      IF (bcnvf == 0 .AND. bit == bitmx) THEN
         mgcnvf = 0
         EXIT
      END IF
   END DO

   IF (mgcnvf == 1) THEN
      IF (irank == root) WRITE(8,'(1X,A,I5)') "B) Multigrid converged in this many v-cycles: ", vit
      IF (irank == root) WRITE(8,'(1X,A,I5)') "B) Total iterations: ", cnt
      CALL CPU_TIME(tsolve)
      EXIT
   END IF
   IF (mgcnvf == 0 .AND. vit == vitmx) THEN
      IF (irank == root) THEN
         PRINT *, "Did not converge"
      END IF
      CALL CPU_TIME(tsolve)
      EXIT
   END IF
 
END DO

IF (irank == root) THEN
   WRITE(8,'(1X,A)') "Integral discrete ordinates..."
   IF (idos == 0) WRITE(8,'(3X,A)') "...with direct solution..."
   IF (idos == 1) WRITE(8,'(3X,A)') "...with CG iterations..."
END IF


CALL idogthr




! Deallocate variables no longer needed
IF (irank==root) DEALLOCATE(fold)
DEALLOCATE(phiold,sv,src,kmat,jpsi,kpsi)
IF (idos == 0 .AND. invf == 0) THEN
   DEALLOCATE(jmat)
   IF (sym == 0) DEALLOCATE(piv)
END IF
IF (idos == 1) THEN
   DEALLOCATE(jmat)
   IF (pcf == 1) DEALLOCATE(dv)
END IF

RETURN
END SUBROUTINE solve
