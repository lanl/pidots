SUBROUTINE output

!-------------------------------------------------------------
!
!    Echo the flux output
!
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER :: i, j, k, l, g, quad, n

! Start the echo of the output for each group
DO g = 1, ng   
   ! Check if the flux converged
      WRITE (8,112) "========== Group ", g, " Converged Scalar Flux Zero Moment =========="
      k = 0
      l = 0
      DO j = 1, nyt
         WRITE (8,*) " Row(j) : ", j
         WRITE (8,113) (flux(i,j,k,l,g), i = 1, nxt)
      END DO
      ! Print the optional flux moments as determined by momp
      IF (momp > 0) THEN
         DO l = 0, momp
            DO k = 0, momp
               IF (k == 0 .AND. l == 0) CYCLE
               IF (k > iall .OR. l > iall) THEN
                  warn = warn + 1
                  WRITE (8,'(/,1X,A,/)') "WARNING: the printed flux moment below is outside the converged orders"
               END IF
               WRITE (8,*)  
               WRITE (8,114) "----- Group: ", g, ", X-Moment: ", k, ", Y-Moment: ", l, " Scalar Flux -----"
               DO j = 1, nyt
                  WRITE (8,*) " Row(j) : ", j
                  WRITE (8,113) (flux(i,j,k,l,g), i = 1, nxt)
               END DO
            END DO
         END DO
      END IF
END DO

!! Call for the printed outward angular fluxes at the boundaries if desired
!IF (pmoaf /= 0) THEN
!   DO g = 1, ng
!      WRITE(8,*)
!      WRITE(8,112) "===== Group ", g, " Zero Moment Angular Flux Outward at Boundaries ====="
!      DO quad = 1, 4
!         DO n = 1, apo
!            WRITE(8,117) "---Quadrant:", quad, " Angle:", n, " y-dir. angular flux---"
!            WRITE(8,116) (psio(i,n,quad,ng), i = 1, nx*order, order)
!            WRITE(8,117) "---Quadrant:", quad, " Angle:", n, " x-dir. angular flux---"
!            WRITE(8,116) (psio(i,n,quad,ng), i = (nx*order+1), (nx+ny)*order, order)
!         END DO
!      END DO
!   END DO
!END IF
!! Print the non-zero moments if desired
!IF (pmoaf == 2) THEN
!   DO g = 1, ng
!      DO k = 2, order
!         WRITE(8,*)
!         WRITE(8,117) "=== Group ", g, " Moment ", k, " Angular Flux Outward at Boundaries ==="
!         DO quad = 1, 4
!            DO n = 1, apo
!               WRITE(8,117) "---Quadrant:", quad, " Angle:", n, " y-dir. angular flux---"
!               WRITE(8,116) (psio(i,n,quad,ng), i = k, nx*order, order)
!               WRITE(8,117) "---Quadrant:", quad, " Angle:", n, " x-dir. angular flux---"
!               WRITE(8,116) (psio(i,n,quad,ng), i = (nx*order+k), (nx+ny)*order, order)
!            END DO
!         END DO
!      END DO
!   END DO
!END IF

112 FORMAT(1X, A, I4, A)
113 FORMAT(2X, 8ES14.6)
114 FORMAT(1X, A, I3, A, I3, A, I3, A)
115 FORMAT(1X, A, I4, A)
116 FORMAT(2X, 8ES14.6)
117 FORMAT(2X, A, I3, A, I3, A)

! Report the number of warnings
WRITE (8,'(//,1X,I2,A,/)') warn, " warnings have been issued"
      
! End the input echo
WRITE (8,*)
WRITE (8,*) "-------------------------- END SOLUTION ----------------------------------"
WRITE (8,*)

RETURN
END SUBROUTINE output
