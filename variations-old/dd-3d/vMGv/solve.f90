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
INTEGER :: i, j, k, m, bit, ierr, temp, ieq, o, n, jeq
INTEGER :: v, sp, vit, cnvf, cnt, xt, yt, zt, bitmxs
REAL*8, DIMENSION(:), ALLOCATABLE :: rphi
REAL*8, DIMENSION(:,:), ALLOCATABLE :: rpsi

INCLUDE 'mpif.h'

neq  = nx*ny*nz
cneq = neq/8
xys  = nx*ny
xzs  = nx*nz
yzs  = ny*nz
bcs  = apo*(xys+xzs+yzs)
cbcs = apo*(xys+xzs+yzs)/4
bcs2 = xys+xzs+yzs

! Allocate ITMM operators and vectors
ALLOCATE(phi(neq,ns),phiold(neq),sv(neq,ns),av(bcs,8,ns),phir(neq,ns),piv(neq,ns))
ALLOCATE(oldphi(neq))
ALLOCATE(psio(bcs,8,ns),oldpsi(bcs,8))
ALLOCATE(jmat(neq,neq,ns))
ALLOCATE(kmat(bcs,neq,8,ns),jpsi(neq,bcs,8,ns))
ALLOCATE(kpsi(bcs2,bcs2,apo,8,ns))

! Allocate temporary residual vectors
ALLOCATE(rphi(neq))
ALLOCATE(rpsi(bcs,8))

IF (irank == root) ALLOCATE(flux(nxt,nyt,nzt))

! Mark the beginning of the solution phase
IF (irank == root) THEN
   WRITE (8,*)
   WRITE (8,*) "-------------------------- THE SOLUTION ----------------------------------"
   WRITE (8,*)
END IF

! Get the time to reach this point
CALL CPU_TIME(ttosolve)

! Need to construct matrices for solving first
! Progress down V-Cycle stages to construct operators once
DO v = 1, ns
   sp = 2**(v-1)
   xt = MOD(xrank,sp)
   yt = MOD(yrank,sp)
   zt = MOD(zrank,sp)
   IF (xt == 0 .AND. yt == 0 .AND. zt == 0) CALL idomats(v,sp)
END DO
CALL CPU_TIME(tjmat)

! Initialize phiold, the iteration counter, the convergence flag
phiold = 0.0
cnt = 0
cnvf  = 0

! Start the multigrid iterations
DO vit = 1, vitmx

   !--------------------------------  Pre-Correction  -------------------------------------------
   ! Work down with bitmx1 iterations the V-Cycle stages until the 
   ! last stage where it is solved in one iteration exactly
   DO v = 1, ns-1
      sp = 2**(v-1)
      xt = MOD(xrank,sp)
      yt = MOD(yrank,sp)
      zt = MOD(zrank,sp)
      ! Call for the solution loop only if your rank participates in the current stage

if (v /= 1) then 
   bitmxs = bitmx1 !/3
else
   bitmxs = bitmx1
end if


      IF (xt == 0 .AND. yt == 0 .AND. zt == 0) THEN

         ! Reset the psii for v>1 because the RHS vector has changed
         IF (v > 1) psii(:,:,v) = 0.0
oldpsi = 0.0
         DO bit = 1, bitmxs
            cnt = cnt + 1
            CALL idot(v,bit)

            ! Now have scalar flux and angular flux
            ! Call for convergence check and send angular flux moments to neighbors
            CALL pbj(v,sp,bit,bitmxs,cnvf)

            ! Use bcnvf to exit V-cycles if converged
            IF (cnvf == 1) EXIT
         END DO

         IF (cnvf == 1) EXIT

         ! Call for the residual computation
         CALL residual(v,rphi,rpsi)

         ! Call to restrict
         CALL restrict(v,sp,rphi,rpsi, vit)
      END IF
   END DO

   IF (cnvf == 1) THEN
      IF (irank == root) WRITE(8,'(1X,A,I5)') "A) Multigrid converged in this many v-cycles: ", vit
      IF (irank == root) WRITE(8,'(1X,A,I5)') "A) Total iterations: ", cnt
      CALL CPU_TIME(tsolve)
      EXIT
   END IF

   !--------------------------------  Coarsest Grid ---------------------------------------------
   ! Only need to solve the coarsest grid in one iteration since it is a 1x1x1 problem
   IF (nrank == 0) CALL idot(ns,1)

!   IF (nrank == 0) print*, vit, phi(:,ns)

!if (nrank ==0 .or. nrank == 42) print*, "pre", vit, nrank, sv(:,2)
!if (nrank ==0) print*, "pre", vit, nrank, phi(:,3)

!if (nrank == 0 .or. nrank == 63) print*, "pre", vit, nrank, phi(:,1)
!if (nrank == 0 .or. nrank == 63) print*, "pre", vit, nrank, phir(:,1)
!if (nrank == 0 .or. nrank == 42) print*, "pre", vit, nrank, phir(:,2)



   !--------------------------------  Post-Correction  ------------------------------------------
   ! Work up with bitmx2 iterations the V-Cycle stages until the 
   ! top stage (fine grid) where convergence is tested. Start at next to coarsest grid.
   DO v = ns-1, 1, -1
      sp = 2**(v-1)
      xt = MOD(xrank,sp)
      yt = MOD(yrank,sp)
      zt = MOD(zrank,sp)
      ! Call for the solution loop only if your rank participates in the current stage
      IF (xt == 0 .AND. yt == 0 .AND. zt == 0) THEN


!if (vit == 2 .and. v==1 .and. nrank==10) print*, nrank, phi(:,v)
!if (nrank == 0 .and. v == 1) print*, "a", nrank, phi(:,2)

if (vit==5 .and. v==1) then
   oldphi = phi(:,v)
end if

         ! Interpolation step: send out correction info to finer grid, add it to phi(:,v)
         CALL interp(v,sp,vit)

!if ((nrank == 0 .or. nrank == 63) .and. v==1) print*, "post", vit, nrank, phi(:,1)


if (vit==5 .and. v==1) then
   oldphi = phi(:,v) - oldphi
!   print *, nrank, maxval(abs(oldphi))
end if


!if (v==1 .and. nrank == 0) print*, nrank, vit, phi(:,2)

!IF (vit == 2 .and. v==1 .and. nrank ==10) print*, nrank, phi(:,v)

         ! Now perform the iterations
         DO bit = 1, bitmx2
            cnt = cnt + 1
            CALL igidot(v,bit, vit,sp)
!IF (vit == 2 .and. v==1 .and. nrank ==10 .and. bit == 1) print*, nrank, phi(:,v)
!IF (vit == 2 .and. v==1 .and. nrank ==10 .and. bit == 2) print*, nrank, phi(:,v)

            CALL pbj(v,sp,bit,bitmx2,cnvf)

            ! Use bcnvf to exit V-cycles if converged
            IF (cnvf == 1) EXIT
         END DO
         IF (cnvf == 1) EXIT
      END IF
   END DO

   ! Final checks on convergence
   IF (cnvf == 1) THEN
      IF (irank == root) WRITE(8,'(1X,A,I5)') "B) Multigrid converged in this many v-cycles: ", vit
      IF (irank == root) WRITE(8,'(1X,A,I5)') "B) Total iterations: ", cnt
      CALL CPU_TIME(tsolve)
      EXIT
   END IF
   IF (cnvf == 0 .AND. vit == vitmx) THEN
      IF (irank == root) THEN
         PRINT *, "Did not converge"
      END IF
      CALL CPU_TIME(tsolve)
      EXIT
   END IF

! Ends the V-cycle iterative sequence
END DO

!if (nrank == 10) print *, nrank, phi(:,1)

!------------------------------------------------------------------------------------------------

IF (irank == root) THEN
   WRITE(8,'(1X,A)') "Integral discrete ordinates..."
   WRITE(8,'(3X,A)') "...with direct solution..."
END IF

CALL idogthr

! Deallocate variables no longer needed
DEALLOCATE(phiold,sv,av,phir,piv,psio,oldpsi,psii)
DEALLOCATE(jmat,kmat,jpsi,kpsi)
DEALLOCATE(rphi,rpsi)

RETURN
END SUBROUTINE solve
