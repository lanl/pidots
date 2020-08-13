SUBROUTINE idocomb

!-------------------------------------------------------------
!
! Gathers all the phi vectors and puts solution into flux 
! 
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER :: i, j, k, ii, jj, kk, rsx, rsy, rsz, ieq, sdi

DO k = 1, nzsd
   rsz = (k-1)*sdnz
   DO j = 1, nysd
      rsy = (j-1)*sdny
      DO i = 1, nxsd
         rsx = (i-1)*sdnx
         sdi = (k-1)*nxsd*nysd + (j-1)*nxsd + i
         DO kk = 1, sdnz
            DO jj = 1, sdny
               DO ii = 1, sdnx
                  ieq = (kk-1)*xys + (jj-1)*sdnx + ii
                  flux(ii+rsx,jj+rsy,kk+rsz) = phi(ieq,sdi)
               END DO
            END DO
         END DO
      END DO
   END DO
END DO

RETURN
END SUBROUTINE idocomb
