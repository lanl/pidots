SUBROUTINE inner(g,its,nsm)

!-------------------------------------------------------------
!
!  Directs the inner iterations
!   Calls for the mesh sweep in 'sweep'
!   Evaulates convergence
! 
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER, INTENT(IN) :: g, nsm
INTEGER :: i, it
INTEGER, INTENT(OUT) :: its
REAL*8 :: df, dfmx
REAL*8, DIMENSION(nx*sord) :: e

! Initialize the previous flux iterate
e = 0.0
! Start the iterations
DO it = 1, itmx
   ! Call for the mesh sweep
   CALL sweep(g,e)
   
   ! Compare new and old flux iterates for user chosen range of moments, iall
   dfmx = -1.0
   DO i = 1, nsm, sord     ! Set stride so only check 0 moment convergence
      ! Compute the difference depending on 'e' value
      IF (ABS(e(i)) >= tolr) THEN
         df = ABS((f(i,g) - e(i))/e(i))
      ELSE
         df = ABS(f(i,g) - e(i))
      END IF
      ! Find the largest value
      IF (df > dfmx) THEN
         dfmx = df
      END IF
   END DO

   ! Print whether or not convergence was reached
   IF (dfmx > err .AND. it < itmx) THEN
      ! Set previous iterate of flux equal to current iterate
      e = f(:,g)
   ELSE IF (dfmx < err) THEN
      cnvf(g) = 1
      its = it
      EXIT
   ELSE IF (it == itmx) THEN
      cnvf(g) = 0
      its = it
      EXIT
   END IF

! End the iterations
END DO

! 111 FORMAT(2X,'Gr',I3,' It ',I5,' Pos ',3I4,' DfMx ',ES11.3,' Flx ',ES11.3, '  Time(s) ', F9.3)

RETURN
END SUBROUTINE inner
