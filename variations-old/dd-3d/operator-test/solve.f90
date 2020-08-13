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
INTEGER :: i, j, k, m, bit, ierr, temp, ieq
REAL*8 :: xsct
REAL*8, INTENT(OUT) :: ttosolve, tjmat, tsolve

INCLUDE 'mpif.h'

neq = nx*ny*nz
xys = nx*ny
xzs = nx*nz
yzs = ny*nz
bcs = apo*(xys+xzs+yzs)
bcs2 = xys+xzs+yzs

! Allocate operators/vectors
ALLOCATE(phi(neq))
ALLOCATE(phiold(neq),src(neq),sv(neq))
ALLOCATE(kmatz(xys,neq,apo,8),kmaty(xzs,neq,apo,8),kmatx(yzs,neq,apo,8))
ALLOCATE(jpsiz(neq,xys,apo,8),jpsiy(neq,xzs,apo,8),jpsix(neq,yzs,apo,8))

ALLOCATE(kpsizz(xys,xys,apo,8),kpsizy(xzs,xys,apo,8),kpsizx(yzs,xys,apo,8))
ALLOCATE(kpsiyz(xys,xzs,apo,8),kpsiyy(xzs,xzs,apo,8),kpsiyx(yzs,xzs,apo,8))
ALLOCATE(kpsixz(xys,yzs,apo,8),kpsixy(xzs,yzs,apo,8),kpsixx(yzs,yzs,apo,8))

! Allocate the outward angular flux value (both SI and ITM use this)
ALLOCATE(psioz(xys,apo,8),psioy(xzs,apo,8),psiox(yzs,apo,8))

! Mark the beginning of the solution phase
IF (irank == root) THEN
   WRITE (8,*)
   WRITE (8,*) "-------------------------- THE SOLUTION ----------------------------------"
   WRITE (8,*)
END IF

! Get the time to reach this point
ttosolve = MPI_WTIME()

! Need to construct matrices for solving first
CALL idomats
phiold = 0.0
tjmat = MPI_WTIME()

! Iterations of the parallel blocks
DO bit = 1, bitmx

   CALL idot(bit)

   ! Now have a scalar flux and angular flux moments
   ! Call for a convergence check and send angular flux moments if necessary
   CALL pbj(bit,tsolve)

   IF (bcnvf == 1) EXIT
   IF (bcnvf == 0 .AND. bit == bitmx) THEN
      CALL MPI_FINALIZE(ierr)
      STOP
   END IF
END DO

IF (irank == root) THEN
   WRITE(8,'(1X,A)') "Integral discrete ordinates..."
   WRITE(8,'(3X,A)') "...with direct solution..."
END IF

! Deallocate variables no longer needed
DEALLOCATE(phiold,src,sv)
DEALLOCATE(jmat)
DEALLOCATE(piv)
DEALLOCATE(kmatz,kmaty,kmatx,jpsiz,jpsiy,jpsix)
DEALLOCATE(kpsizz,kpsizy,kpsizx)
DEALLOCATE(kpsiyz,kpsiyy,kpsiyx)
DEALLOCATE(kpsixz,kpsixy,kpsixx)

! Check the scalar flux printing flag
IF (sfp == 1) THEN
   ! Root allocates the primary flux solution holder
   IF (irank == root) ALLOCATE(flux(nxt,nyt,nzt))

   ! Bring all the blocks' solutions for all the groups back to root for output
   CALL idogthr
END IF

RETURN
END SUBROUTINE solve
