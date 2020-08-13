SUBROUTINE solve(ttosolve,tjmat,tsolve)

!-------------------------------------------------------------
!
!  Directs the solution by either calling for a mesh
!   sweep (DD-SI) or for a ITM solution (DD-ITM)
! 
!-------------------------------------------------------------

USE invar
USE itmvar
USE totvar
IMPLICIT NONE
INTEGER :: i, j, k, bit, ierr, temp, ieq, i1, sdi, cnvf, cnt, v, vit, bitsp, bitd
INTEGER :: spx, spy, spz, kc, jc, ic, crsf, cneq, cbcs, xt, yt, zt
!INTEGER, DIMENSION(neq,nsdp) :: piv
!INTEGER, DIMENSION(neq,ns1) :: fpiv
REAL*8 :: sdcf
REAL*8, INTENT(OUT) :: ttosolve, tjmat, tsolve
!REAL*8, DIMENSION(neq) :: frphi
!REAL*8, DIMENSION(bcs,8) :: frpsi, foldpsi
!REAL*8, DIMENSION(neq,nsdp) :: phi, phiold, phir, sv, rphi
!REAL*8, DIMENSION(neq,ns1) :: f, fr, fsv
!REAL*8, DIMENSION(bcs,8,nsdp) :: psii, psio, av, oldpsi, rpsi
!REAL*8, DIMENSION(bcs,8,ns1) :: fpi, fpo, fav
!REAL*8, DIMENSION(neq,neq,nsdp) :: jmat
!REAL*8, DIMENSION(neq,neq,ns1) :: fjmt
!REAL*8, DIMENSION(bcs,neq,8,nsdp) :: kmat
!REAL*8, DIMENSION(bcs,neq,8,ns1) :: fkmt
!REAL*8, DIMENSION(neq,bcs,8,nsdp) :: jpsi
!REAL*8, DIMENSION(neq,bcs,8,ns1) :: fjps
!REAL*8, DIMENSION(bcs2,bcs2,apo,8,nsdp) :: kpsi
!REAL*8, DIMENSION(bcs2,bcs2,apo,8,ns1) :: fkps

INCLUDE 'mpif.h'

! Allocate
ALLOCATE(piv(neq,nsdp),fpiv(neq,ns1))
ALLOCATE(frphi(neq),frpsi(bcs,8),foldpsi(bcs,8))
ALLOCATE(phi(neq,nsdp),phiold(neq,nsdp),phir(neq,nsdp),sv(neq,nsdp),rphi(neq,nsdp))
ALLOCATE(f(neq,ns1),fr(neq,ns1),fsv(neq,ns1))
ALLOCATE(psii(bcs,8,nsdp),psio(bcs,8,nsdp),av(bcs,8,nsdp),oldpsi(bcs,8,nsdp),rpsi(bcs,8,nsdp))
ALLOCATE(fpi(bcs,8,ns1),fpo(bcs,8,ns1),fav(bcs,8,ns1))
ALLOCATE(jmat(neq,neq,nsdp),fjmt(neq,neq,ns1))
ALLOCATE(kmat(bcs,neq,8,nsdp),fkmt(bcs,neq,8,ns1))
ALLOCATE(jpsi(neq,bcs,8,nsdp),fjps(neq,bcs,8,ns1))
ALLOCATE(kpsi(bcs2,bcs2,apo,8,nsdp),fkps(bcs2,bcs2,apo,8,ns1))

! Mark the beginning of the solution phase
IF (irank == root) THEN
   WRITE (8,*)
   WRITE (8,*) "-------------------------- THE SOLUTION ----------------------------------"
   WRITE (8,*)
END IF

! Get the time to reach this point
ttosolve = MPI_WTIME()

! Need to construct matrices for solving first
! Progress down V-cycle stages to construct operators once
spx = 1
spy = 1
spz = 1
DO v = 0, ns1
   IF (v == 0) THEN
      CALL idomats(spx,spy,spz,jmat,kmat,jpsi,kpsi,piv,sv,av)
   ELSE
      xt = MOD(xrank,spx)
      yt = MOD(yrank,spy)
      zt = MOD(zrank,spz)
      IF (xt == 0 .AND. yt == 0 .AND. zt == 0) THEN
         CALL idomatsf(spx,spy,spz,fjmt(:,:,v),fkmt(:,:,:,v),fjps(:,:,:,v),fkps(:,:,:,:,v),fpiv(:,v))
      END IF
      IF (v < ns1) THEN
         spx = spx*xcf(v)
         spy = spy*ycf(v)
         spz = spz*zcf(v)
      END IF
   END IF
END DO

! Initialize psii, phiold, the iteration counter, the convergence flag
psii = 0.0
phiold = 0.0
oldpsi = 0.0
cnt = 0
cnvf = 0

! Get the construction time
tjmat = MPI_WTIME()

! Hardwire max iterations for v>0
bitsp = 1

! Set dummy bitmx2
bitd = bitmx2 + 1

! Start the Multigrid iterations
DO vit = 1, vitmx

   !------------------------------------ Pre-Correction -------------------------------------------
   ! Work down with bitmx1 iterations the V-cycle stages until the last stage
   ! where it is solved in one iteration exactly
   ! Start with the coarsening step of multiple to one s-d per processor: STAGE 0
   spx = 1
   spy = 1
   spz = 1
   DO v = 0, ns2
      kc = zcf(v)
      jc = ycf(v)
      ic = xcf(v)
      sdcf = 1.0
      crsf = 0
      IF (kc > sdnz) THEN
         crsf = 1
         sdcf = sdcf*REAL(sdnz)/REAL(kc)
         kc = sdnz
      END IF
      IF (jc > sdny) THEN
         crsf = 1
         sdcf = sdcf*REAL(sdny)/REAL(jc)
         jc = sdny
      END IF
      IF (ic > sdnx) THEN
         crsf = 1
         sdcf = sdcf*REAL(sdnx)/REAL(ic)
         ic = sdnx
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

      ! All processes do PGS iterations and restrict to one s-d per processor
      IF (v == 0) THEN
         DO bit = 1, bitmx1
            cnt = cnt + 1
            ! Red iterations
            sdi = 0
            DO k = 1, nzsd
               DO j = 1, nysd
                  i1 = MOD(j+k,2) + 1
                  DO i = i1, nxsd, 2
                     sdi = (k-1)*nysd*nxsd + (j-1)*nxsd + i
                     CALL idot(phi(:,sdi),sv(:,sdi),psii(:,:,sdi),psio(:,:,sdi),av(:,:,sdi), &
                               jmat(:,:,sdi),kmat(:,:,:,sdi),jpsi(:,:,:,sdi),kpsi(:,:,:,:,sdi),piv(:,sdi))
                  END DO
               END DO
            END DO
            ! Copy/Send from red to black
            CALL pgsred(bit,bitmx1,psii,psio,oldpsi)
               
            ! Black iterations
            sdi = 0
            DO k = 1, nzsd
               DO j = 1, nysd
                  i1 = MOD(j+k+1,2) + 1
                  DO i = i1, nxsd, 2
                     sdi = (k-1)*nysd*nxsd + (j-1)*nxsd + i
                     CALL idot(phi(:,sdi),sv(:,sdi),psii(:,:,sdi),psio(:,:,sdi),av(:,:,sdi), &
                               jmat(:,:,sdi),kmat(:,:,:,sdi),jpsi(:,:,:,sdi),kpsi(:,:,:,:,sdi),piv(:,sdi))
                  END DO
               END DO
            END DO
            ! Check for convergence and pass from black to red
            CALL pgsblk(bit,bitmx1,cnvf,phi,phiold,psii,psio,oldpsi)

            ! Exit the bit do loop if converged
            IF (cnvf == 1) EXIT

!print *, bit
!print *, phi
         END DO

         ! Exit the v loop
         IF (cnvf == 1) EXIT

         ! Call for the residual computation
         CALL residual(psii,oldpsi,rphi,rpsi,kmat,kpsi)


!do bit = 1, nsdp
!   print *, bit
!   print *, rphi(:,bit)
!end do

         ! Call to restrict to one sub-domain per processor
         CALL restrict(kc,jc,ic,cneq,cbcs,sdcf,rphi,rpsi,phi,phir,fsv(:,1),fav(:,:,1))

      ! Otherwise only certain processes move on
      ELSE
         ! Determine if your rank says you'll move on in this stage
         xt = MOD(xrank,spx)
         yt = MOD(yrank,spy)
         zt = MOD(zrank,spz)

         ! Call for the solution loop only if your rank participates in the current stage
         IF (xt == 0 .AND. yt == 0 .AND. zt == 0) THEN
            ! Reinitialize
            foldpsi = 0.0
            DO bit = 1, bitsp
               cnt = cnt + 1
               ! Call for the residual equation iteration
               CALL idotf(f(:,v),fsv(:,v),fpiv(:,v),fpo(:,:,v),fav(:,:,v),fjmt(:,:,v),fjps(:,:,:,v))
               ! Now call for communication in a PBJ sense because v>0
               CALL pbj(spx,spy,spz,fpi(:,:,v),fpo(:,:,v))
            END DO

            ! Call for the residual computation and restriction
            CALL residualf(fpi(:,:,v),foldpsi,frphi,frpsi,fkmt(:,:,:,v),fkps(:,:,:,:,v))
            CALL restrictf(v,spx,spy,spz,kc,jc,ic,cneq,cbcs,crsf,sdcf,frphi,frpsi,f(:,v),fr(:,v),fsv(:,v+1),fav(:,:,v+1))
         END IF

         ! Update spx/y/z for next loop
         spx = spx*xcf(v)
         spy = spy*ycf(v)
         spz = spz*zcf(v)
     END IF
   END DO

   ! Break out of the vit loops
   IF (cnvf == 1) THEN
      IF (irank == root) THEN
         WRITE(8,'(1X,A,I5)') "A) Multigrid converged in this many v-cycles: ", vit
         WRITE(8,'(1X,A,I5)') "A) Total iterations: ", cnt
      END IF
      tsolve = MPI_WTIME()
      EXIT
   END IF

   !----------------------------------------- Coarsest Grid ----------------------------------------------
   IF (nrank == 0) THEN
      CALL idotcg(fpiv(:,ns1),f(:,ns1),fsv(:,ns1),fjmt(:,:,ns1))
      cnt = cnt + 1
   END IF

!if (nrank == root) print *, fsv(:,1)
!if (nrank == root) print *, f(:,1)

   !---------------------------------------- Post-Correction ---------------------------------------------
   ! Work back up the V-Cycle stages until the top stage (fine grid) where
   ! convergence is tested. Start at next coarsest grid
   DO v = ns2, 0, -1
      kc = zcf(v)
      jc = ycf(v)
      ic = xcf(v)
      crsf = 0
      IF (kc > sdnz) THEN
         crsf = 1
         kc = sdnz
      END IF
      IF (jc > sdny) THEN
         crsf = 1
         jc = sdny
      END IF
      IF (ic > sdnx) THEN
         crsf = 1
         ic = sdnx
      END IF
      cneq = neq/(kc*jc*ic)

      IF (v > 0) THEN
         ! Adjust parallel parameters
         spx = spx/xcf(v)
         spy = spy/ycf(v)
         spz = spz/zcf(v)
         xt = MOD(xrank,spx)
         yt = MOD(yrank,spy)
         zt = MOD(zrank,spz)
         ! Call for the solution loop only if your rank participates in the current stage
         IF (xt == 0 .AND. yt == 0 .AND. zt == 0) THEN
            ! Interpolation step: send out correction info to finer grid, add it to f(:,v)
            CALL interpf(v,spx,spy,spz,kc,jc,ic,cneq,f(:,v),f(:,v+1),fr(:,v))

            ! Already interpolated the flux; add to the counter
            cnt = cnt + 1
         END IF

      ! Otherwise interpolate back to multiple sub-domains per processor and
      ! iterate with bitmx2 iterations
      ELSE

!print *, phi
        ! All processes interpolate to many sub-domains per p
         CALL interp(kc,jc,ic,cneq,phi,f(:,v+1),phir)

!print *, phi

         ! Now perform the iterations since v=0
         DO bit = 1, bitmx2
            cnt = cnt + 1
            ! Red iterations
            sdi = 0
            DO k = 1, nzsd
               DO j = 1, nysd
                  i1 = MOD(j+k,2) + 1
                  DO i = i1, nxsd, 2
                     sdi = (k-1)*nysd*nxsd + (j-1)*nxsd + i
                     CALL idotig(bit,phi(:,sdi),sv(:,sdi),psii(:,:,sdi),psio(:,:,sdi),av(:,:,sdi), &
                                 jmat(:,:,sdi),kmat(:,:,:,sdi),jpsi(:,:,:,sdi),kpsi(:,:,:,:,sdi),piv(:,sdi))
                  END DO
               END DO
            END DO
            ! Copy/Send from red to black
            CALL pgsred(bit,bitd,psii,psio,oldpsi)
               
            ! Black iterations
            sdi = 0
            DO k = 1, nzsd
               DO j = 1, nysd
                  i1 = MOD(j+k+1,2) + 1
                  DO i = i1, nxsd, 2
                     sdi = (k-1)*nysd*nxsd + (j-1)*nxsd + i

!                     CALL idot(phi(:,sdi),sv(:,sdi),psii(:,:,sdi),psio(:,:,sdi),av(:,:,sdi), &
!                               jmat(:,:,sdi),kmat(:,:,:,sdi),jpsi(:,:,:,sdi),kpsi(:,:,:,:,sdi),piv(:,sdi))


                     CALL idotig(bit,phi(:,sdi),sv(:,sdi),psii(:,:,sdi),psio(:,:,sdi),av(:,:,sdi), &
                                  jmat(:,:,sdi),kmat(:,:,:,sdi),jpsi(:,:,:,sdi),kpsi(:,:,:,:,sdi),piv(:,sdi))
                  END DO
               END DO
            END DO
            ! Check for convergence and pass from black to red
            ! Don't need oldpsi, so just use bitd, dummy value
            CALL pgsblk(bit,bitd,cnvf,phi,phiold,psii,psio,oldpsi)

            ! Exit the bit do loop if converged
            IF (cnvf == 1) EXIT
         END DO

         ! Exit the v do loop if converged
         IF (cnvf == 1) EXIT
      END IF
   END DO

   ! Final checks on convergence
   ! Break out of the vit loops
   IF (cnvf == 1) THEN
      IF (irank == root) THEN
         WRITE(8,'(1X,A,I5)') "B) Multigrid converged in this many v-cycles: ", vit
         WRITE(8,'(1X,A,I5)') "B) Total iterations: ", cnt
      END IF
      tsolve = MPI_WTIME()
      EXIT
   END IF
   IF (cnvf == 0 .AND. vit == vitmx) THEN
      IF (irank == root) WRITE(8,*) "    DID NOT CONVERGE    "
      tsolve = MPI_WTIME()
      EXIT
   END IF

! Ejnd the V-cycle iterative sequence
END DO

!------------------------------------------------------------------------------------------

IF (irank == root) THEN
   WRITE(8,'(1X,A)') "Integral discrete ordinates..."
   WRITE(8,'(3X,A)') "...with direct solution..."
END IF

! Combine all the scalar fluxes across sub-domains on each processor individually
CALL idocomb(phi)

! Check the scalar flux printing flag
IF (sfp == 1) CALL flxgthr 

! Deallocate
DEALLOCATE(piv,fpiv)
DEALLOCATE(frphi,frpsi,foldpsi,phi,phiold,phir,sv,rphi,f,fr,fsv)
DEALLOCATE(psii,psio,av,oldpsi,rpsi,fpi,fpo,fav)
DEALLOCATE(jmat,kmat,jpsi,kpsi)
DEALLOCATE(fjmt,fkmt,fjps,fkps)

RETURN
END SUBROUTINE solve
