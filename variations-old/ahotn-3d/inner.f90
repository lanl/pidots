SUBROUTINE inner(g,its)

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
INTEGER :: i, j, k, t, u, v, it
INTEGER, INTENT(OUT) :: its
REAL*8 :: df, dfmx
REAL*8, DIMENSION(nx,ny,nz,0:lambda,0:lambda,0:lambda) :: e

! Initialize the previous flux iterate
e = 0.0
! Start the iterations
DO it = 1, itmx
   ! Call for the mesh sweep
   CALL sweep(g,e)
   
   ! Compare new and old flux iterates for user chosen range of moments, iall
   dfmx = -1.0
   DO v = 0, iall
      DO u = 0, iall
         DO t = 0, iall
            DO k = 1, nz
               DO j = 1, ny
                  DO i = 1, nx
                     ! Compute the difference depending on 'e' value
                     IF (e(i,j,k,t,u,v) >= tolr) THEN
                        df = ABS((f(i,j,k,t,u,v,g) - e(i,j,k,t,u,v))/e(i,j,k,t,u,v))
                     ELSE
                        df = ABS(f(i,j,k,t,u,v,g) - e(i,j,k,t,u,v))
                     END IF
                     ! Find the largest value
                     IF (df > dfmx) THEN
                        dfmx = df
                     END IF
                  END DO
               END DO
            END DO
         END DO
      END DO
   END DO

   ! Print whether or not convergence was reached
   IF (dfmx > err .AND. it < itmx) THEN
      ! Set previous iterate of flux equal to current iterate
      e = f(:,:,:,:,:,:,g)
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

!111 FORMAT(2X,'Gr',I3,' It ',I5,' Pos ',3I4,' Mom ',3I2,' DfMx ',ES11.3,' Flx ',ES11.3, '  Time(s) ', F9.3)

RETURN
END SUBROUTINE inner
