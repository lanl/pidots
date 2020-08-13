SUBROUTINE idomats(spx,spy,spz,jmat,kmat,jpsi,kpsi,piv,sv,av)

!-------------------------------------------------------------
!
!  Integral discrete ordinates matrices
!
!   Control the construction of jmat, kmat, jpsi, kpsi
!   Matrices nned to be made only once, then idot.f90 will
!    use them to solve the system of equations
!
!-------------------------------------------------------------

USE invar
IMPLICIT NONE
INTEGER, INTENT(IN) :: spx, spy, spz
INTEGER :: ieq, i, j, k, info, m, t, sdi, ii, jj, kk, rsx, rsy, rsz
INTEGER, DIMENSION(neq,nsdp), INTENT(OUT) :: piv
REAL*8, DIMENSION(neq) :: src
REAL*8, DIMENSION(neq,nsdp), INTENT(OUT) :: sv
REAL*8, DIMENSION(bcs,8,nsdp), INTENT(OUT) :: av
REAL*8, DIMENSION(neq,neq,nsdp), INTENT(OUT) :: jmat
REAL*8, DIMENSION(bcs,neq,8,nsdp), INTENT(OUT) :: kmat
REAL*8, DIMENSION(neq,bcs,8,nsdp), INTENT(OUT) :: jpsi
REAL*8, DIMENSION(bcs2,bcs2,apo,8,nsdp), INTENT(OUT) :: kpsi

! Initialize the matrices
jmat = 0.0
kmat = 0.0
jpsi = 0.0
kpsi = 0.0

!------------------------------------------------------------------
! Want to make all the operators for each sub-domain one at a time
sdi = 0
DO k = 1, nzsd
   rsz = (k-1)*sdnz
   DO j = 1, nysd
      rsy = (j-1)*sdny
      DO i = 1, nxsd
         sdi = sdi + 1
         rsx = (i-1)*sdnx

         ! Can now start the single sweep for constructing all matrices
         CALL matsweep(rsx,rsy,rsz,spx,spy,spz,jmat(:,:,sdi),kmat(:,:,:,sdi),jpsi(:,:,:,sdi),kpsi(:,:,:,:,sdi))

         ! Begin constructing the RHS
         DO kk = 1, sdnz
            DO jj = 1, sdny
               DO ii = 1, sdnx
                  ieq = ii + (jj-1)*sdnx + (kk-1)*xys
                  src(ieq) = s(ii+rsx,jj+rsy,kk+rsz)/sigs(mat(ii+rsx,jj+rsy,kk+rsz))
               END DO
            END DO
         END DO

         ! Multiply the source vector by the Jacobi iteration matrix
         DO jj = 1, neq
            sv(jj,sdi) = DOT_PRODUCT(jmat(:,jj,sdi),src)
         END DO

         ! Multiply the source vector by Jpsi to form the other saved vector
         DO jj = 1, 8
            DO ii = 1, bcs
               av(ii,jj,sdi) = DOT_PRODUCT(jpsi(:,ii,jj,sdi),src)
            END DO
         END DO

         ! Form the LHS matrix
         jmat(:,:,sdi) = -jmat(:,:,sdi)
         DO t = 1, neq
            jmat(t,t,sdi) = 1.0 + jmat(t,t,sdi)
         END DO

         !----------------------------------------------------------------
         ! Need to factor jmat

         CALL dgetrf(neq,neq,jmat(:,:,sdi),neq,piv(:,sdi),info)
         IF (info /= 0) THEN
            WRITE (8,'(//,1X,A)') "ERROR: jmat either has illegal value or is singular, cannot factor."
            STOP
         END IF

         ! Have needed values: src, sv, jmat, kmat, jpsi, kpsi, piv

      END DO        ! x sub-domain divisions
   END DO           ! y sub-domain divisions
END DO              ! z sub-domain divisions

!----------------------------------------------------------------------------------------

RETURN
END SUBROUTINE idomats
