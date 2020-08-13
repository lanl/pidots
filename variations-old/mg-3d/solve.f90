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
CHARACTER(8) :: fl
INTEGER :: i, j, k, m, bit, ierr, temp, ieq, o, n, jeq
INTEGER :: v, spx, spy, spz, vit, cnvf, cnt, xt, yt, zt
INTEGER :: kc, jc, ic, cneq, cbcs, bitsp, crsf, tx, ty, tz
REAL, INTENT(OUT) :: ttosolve, tjmat, tsolve
REAL*8 :: sdcf
REAL*8, DIMENSION(:), ALLOCATABLE :: rphi
REAL*8, DIMENSION(:,:), ALLOCATABLE :: rpsi

INCLUDE 'mpif.h'

neq  = nx*ny*nz
xys  = nx*ny
xzs  = nx*nz
yzs  = ny*nz
bcs  = apo*(xys+xzs+yzs)
bcs2 = xys+xzs+yzs

! Allocate ITMM operators and vectors
ALLOCATE(phi(neq,ns),phiold(neq),sv(neq,ns),av(bcs,8,ns),phir(neq,ns),piv(neq,ns))
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
spx = 1
spy = 1
spz = 1
DO v = 1, ns
   xt = MOD(xrank,spx)
   yt = MOD(yrank,spy)
   zt = MOD(zrank,spz)
   IF (xt == 0 .AND. yt == 0 .AND. zt == 0) CALL idomats(v,spx,spy,spz)
   IF (v < ns) THEN
      spx = spx*xcf(v)
      spy = spy*ycf(v)
      spz = spz*zcf(v)
   END IF
END DO
CALL CPU_TIME(tjmat)

! Initialize phiold, the iteration counter, the convergence flag
phiold = 0.0
cnt = 0
cnvf = 0

! Start the multigrid iterations
DO vit = 1, vitmx

   !--------------------------------  Pre-Correction  -------------------------------------------
   ! Work down with bitmx1 iterations the V-Cycle stages until the 
   ! last stage where it is solved in one iteration exactly
   spx = 1
   spy = 1
   spz = 1
   DO v = 1, ns-1
      xt = MOD(xrank,spx)
      yt = MOD(yrank,spy)
      zt = MOD(zrank,spz)
      kc = zcf(v)
      jc = ycf(v)
      ic = xcf(v)
      tz = nz
      ty = ny
      tx = nx
      sdcf = 1.0
      crsf = 0
      IF (kc > nz) THEN
         crsf = 1
         sdcf = sdcf*REAL(nz)/REAL(kc)
         tz = kc
         kc = nz
      END IF
      IF (jc > ny) THEN
         crsf = 1
         sdcf = sdcf*REAL(ny)/REAL(jc)
         ty = jc
         jc = ny
      END IF
      IF (ic > nx) THEN
         crsf = 1
         sdcf = sdcf*REAL(nx)/REAL(ic)
         tx = ic
         ic = nx
      END IF
      IF (crsf == 1) THEN
         IF (irank == root .AND. vit == 1) THEN
            warn = warn + 1
            WRITE(8,'(/,3X,A,I2,/)') "WARNING: Sub-domain coarsening factor is greater than &
                                      number of cells in one or more dimensions for stage ", v
         END IF
      END IF
      cneq = neq/(kc*jc*ic)
      cbcs = apo*(xys/(ic*jc) + xzs/(ic*kc) + yzs/(jc*kc))
      ! IF (vit == 1 .AND. v == 1) ALLOCATE(phir2(tx,ty,tz,ns))

      ! Call for the solution loop only if your rank participates in the current stage
      IF (xt == 0 .AND. yt == 0 .AND. zt == 0) THEN
         ! Reset the psii for v>1 because the RHS vector has changed
         IF (v > 1) psii(:,:,v) = 0.0
         oldpsi = 0.0

         IF (v > 1) THEN
            bitsp = 1
         ELSE
            bitsp = bitmx1
         END IF

         DO bit = 1, bitsp
            cnt = cnt + 1
            IF (v > 1) THEN
               CALL idotv(v,bit)
            ELSE
               CALL idot(v,bit)
            END IF

            ! Now have scalar flux and angular flux
            ! Call for convergence check and send angular flux moments to neighbors
            CALL pbj(v,spx,spy,spz,bit,bitsp,cnvf)

            ! Use cnvf to exit V-cycles if converged
            IF (cnvf == 1) EXIT
         END DO

         IF (cnvf == 1) EXIT

         ! Call for the residual computation
         CALL residual(v,rphi,rpsi)
         ! Restrict the residual
         CALL restrict(v,spx,spy,spz,kc,jc,ic,cneq,cbcs,crsf,sdcf,rphi,rpsi)
      END IF

      spx = spx*xcf(v)
      spy = spy*ycf(v)
      spz = spz*zcf(v)
   END DO

   IF (cnvf == 1) THEN
      IF (irank == root) WRITE(8,'(1X,A,I5)') "A) Multigrid converged in this many v-cycles: ", vit
      IF (irank == root) WRITE(8,'(1X,A,I5)') "A) Total iterations: ", cnt
      CALL CPU_TIME(tsolve)
      EXIT
   END IF

   !--------------------------------  Coarsest Grid ---------------------------------------------
   ! Only need to solve the coarsest grid in one iteration since it is a 1x1x1 problem
   IF (nrank == 0) THEN
      CALL cgidot(ns)
      cnt = cnt + 1
   END IF

   !--------------------------------  Post-Correction  ------------------------------------------
   ! Work up with bitmx2 iterations the V-Cycle stages until the 
   ! top stage (fine grid) where convergence is tested. Start at next to coarsest grid.
   DO v = ns-1, 1, -1
      spx = spx/xcf(v)
      spy = spy/ycf(v)
      spz = spz/zcf(v)
      xt = MOD(xrank,spx)
      yt = MOD(yrank,spy)
      zt = MOD(zrank,spz)
      kc = zcf(v)
      jc = ycf(v)
      ic = xcf(v)
      crsf = 0
      IF (kc > nz) THEN
         crsf = 1
         kc = nz
      END IF
      IF (jc > ny) THEN
         crsf = 1
         jc = ny
      END IF
      IF (ic > nx) THEN
         crsf = 1
         ic = nx
      END IF
      cneq = neq/(kc*jc*ic)

      ! Call for the solution loop only if your rank participates in the current stage
      IF (xt == 0 .AND. yt == 0 .AND. zt == 0) THEN
         ! Interpolation step: send out correction info to finer grid, add it to phi(:,v)
         CALL interp(v,spx,spy,spz,kc,jc,ic,cneq)

         IF (v > 1) THEN
            ! Already interpolated the flux; add to the counter
            cnt = cnt + 1
         ELSE
            ! Now perform the iterations if v > 1
            DO bit = 1, bitmx2
               cnt = cnt + 1
               CALL igidot(v,bit)

               CALL pbj(v,spx,spy,spz,bit,bitmx2,cnvf)

               ! Use bcnvf to exit V-cycles if converged
               IF (cnvf == 1) EXIT
            END DO
         END IF
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
         WRITE(8,*) "   DID NOT CONVEGE   "
      END IF
      CALL CPU_TIME(tsolve)
      EXIT
   END IF

! Ends the V-cycle iterative sequence
END DO

!------------------------------------------------------------------------------------------------

IF (irank == root) THEN
   WRITE(8,'(1X,A)') "Integral discrete ordinates..."
   WRITE(8,'(3X,A)') "...with direct solution..."
END IF

CALL idogthr

! Deallocate variables no longer needed
DEALLOCATE(phiold,sv,av,phir,piv,psio,oldpsi,psii)              ! ,phir2)
DEALLOCATE(jmat,kmat,jpsi,kpsi)
DEALLOCATE(rphi,rpsi)

RETURN
END SUBROUTINE solve
