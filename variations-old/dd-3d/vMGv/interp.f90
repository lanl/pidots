SUBROUTINE interp(v,sp, vit)

!------------------------------------------------------------------
!
! Coarse grid processes send correction to fine grid processes.
! Then interpolate the correction using phir.
!
!------------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER, INTENT(IN) :: v, sp, vit
INTEGER :: kt, jt, it, k, j, i, ieq, jeq
REAL*8, DIMENSION(cneq) :: rsv

! First get the correction sent to the fine grid processes
CALL corrscat(v,sp,rsv, vit)

!if ((nrank == 0 .or. nrank==42) .and. v==2) print *, "rsv", vit, nrank, rsv


! Now interpolate
DO k = 1, nz
   kt = (k+1)/2
   DO j = 1, ny
      jt = (j+1)/2
      DO i = 1, nx
         it = (i+1)/2
         ieq = (k-1)*xys + (j-1)*nx + i
         jeq = (kt-1)*xys/4 + (jt-1)*nx/2 + it
         phi(ieq,v) = phi(ieq,v) + phir(ieq,v)*rsv(jeq)
      END DO
   END DO
END DO

!if ((nrank == 0 .or. nrank==42) .and. v==2) print *, vit, nrank, phi(:,2)


RETURN
END SUBROUTINE interp
