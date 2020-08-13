SUBROUTINE idomats

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
USE solvar
IMPLICIT NONE
INTEGER :: ieq, i, j, k, info, m, t, sdi, ii, jj, kk, rsx, rsy, rsz

! Initialize the matrices
jmat = 0.0d0
kmat = 0.0d0
jpsi = 0.0d0
kpsi = 0.0d0

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
         CALL matsweep(sdi,rsx,rsy,rsz)

         ! Begin constructing the RHS
         DO kk = 1, sdnz
            DO jj = 1, sdny
               DO ii = 1, sdnx
                  ieq = ii + (jj-1)*sdnx + (kk-1)*xys
                  src(ieq,sdi) = s(ii+rsx,jj+rsy,kk+rsz)/sigs(mat(ii+rsx,jj+rsy,kk+rsz))
               END DO
            END DO
         END DO

         ! Multiply the source vector by the Jacobi iteration matrix
         DO jj = 1, neq
            sv(jj,sdi) = DOT_PRODUCT(jmat(:,jj,sdi),src(:,sdi))
         END DO

         ! Form the LHS matrix
         jmat(:,:,sdi) = -jmat(:,:,sdi)
         DO t = 1, neq
            jmat(t,t,sdi) = 1.0d0 + jmat(t,t,sdi)
         END DO


!print *, sdi
!print *, jmat

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
