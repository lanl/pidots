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
INTEGER :: i, j, k, bit, ierr, temp, ieq, i1, sdi
REAL*8, INTENT(OUT) :: ttosolve, tjmat, tsolve

INCLUDE 'mpif.h'

! Allocate the solution vectors/matrices
ALLOCATE(phi(neq,nsdp))
ALLOCATE(phiold(neq,nsdp),src(neq,nsdp),sv(neq,nsdp),piv(neq,nsdp))
ALLOCATE(jmat(neq,neq,nsdp))
ALLOCATE(kmat(bcs,neq,8,nsdp), jpsi(neq,bcs,8,nsdp))
ALLOCATE(kpsi(bcs2,bcs2,apo,8,nsdp))

! Allocate the outward angular flux value
ALLOCATE(psio(bcs,8,nsdp))

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
phiold = 0.0d0
tjmat = MPI_WTIME()

! Iterations of the parallel blocks
DO bit = 1, bitmx

   ! Red iterations
   sdi = 0
   DO k = 1, nzsd
      DO j = 1, nysd
         i1 = MOD(j+k,2) + 1
         DO i = i1, nxsd, 2
            sdi = (k-1)*nysd*nxsd + (j-1)*nxsd + i
            CALL idot(sdi,bit)
         END DO
      END DO
   END DO

   ! Pass the angular flux of red sub-domains to black sub-domains
   CALL pgsred

   ! Black iterations
   sdi = 0
   DO k = 1, nzsd
      DO j = 1, nysd
         i1 = MOD(j+k+1,2) + 1
         DO i = i1, nxsd, 2
            sdi = (k-1)*nysd*nxsd + (j-1)*nxsd + i
            CALL idot(sdi,bit)
         END DO
      END DO
   END DO

   ! Now have a scalar flux and angular flux for all red and black sub-domains. 
   ! Call for convergence check and black to red passing
   CALL pgsblk(bit,tsolve)

   IF (bcnvf == 1) EXIT
!   IF (bcnvf == 0 .AND. bit == bitmx) THEN
!      CALL MPI_FINALIZE(ierr)
!   END IF
END DO

IF (irank == root) THEN
   WRITE(8,'(1X,A)') "Integral discrete ordinates..."
   WRITE(8,'(3X,A)') "...with direct solution..."
END IF

! Combine all the scalar fluxes across sub-domains on each processor individually
ALLOCATE(flux(nx,ny,nz))
CALL idocomb

! Deallocate variables no longer needed
DEALLOCATE(phi,phiold,src,sv,piv,jmat,kmat,jpsi,kpsi,psio)

! Check the scalar flux printing flag
IF (sfp == 1) THEN
   ! Root allocates the primary flux solution holder
   IF (irank == root) ALLOCATE(fluxt(nxt,nyt,nzt))
   CALL flxgthr
END IF

RETURN
END SUBROUTINE solve
