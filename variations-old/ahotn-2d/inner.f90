SUBROUTINE inner(g)

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
INTEGER, INTENT(IN) :: g
INTEGER :: i, j, k, l, it
REAL*8 :: df, dfmx

! Initialize the previous flux iterate
e = 0.0
! Start the iterations
DO it = 1, itmx
   ! Call for the mesh sweep
   CALL sweep(g)
   
   ! Compare new and old flux iterates for user chosen range of moments, iall
   dfmx = -1.0
   DO l = 0, iall
      DO k = 0, iall
         DO j = 1, ny
            DO i = 1, nx
               ! Compute the difference depending on 'e' value
               IF (e(i,j,k,l) >= tolr) THEN
                  df = ABS((f(i,j,k,l,g) - e(i,j,k,l))/e(i,j,k,l))
               ELSE
                  df = ABS(f(i,j,k,l,g) - e(i,j,k,l))
               END IF
               ! Find the largest value
               IF (df > dfmx) THEN
                  dfmx = df
               END IF
            END DO
         END DO
      END DO
   END DO
   
   ! Print whether or not convergence was reached
   IF (dfmx > err .AND. it < itmx) THEN
      ! Set previous iterate of flux equal to current iterate
      e = f(:,:,:,:,g)
   ELSE IF (dfmx < err) THEN
      cnvf(g) = 1
      EXIT
   ELSE IF (it == itmx) THEN
      cnvf(g) = 0
      EXIT
   END IF

! End the iterations
END DO

RETURN
END SUBROUTINE inner
