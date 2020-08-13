SUBROUTINE precond(inv,tran,vec)

!------------------------------------------------------------------
!
! Apply the preconditioner to the incoming vector 'vec' depending
!  on the preconditioner flag pinv
!
!------------------------------------------------------------------

USE solvar
IMPLICIT NONE
INTEGER, INTENT(IN) :: inv, tran
INTEGER :: i, info
REAL*8, DIMENSION(neq), INTENT(INOUT) :: vec
REAL*8, DIMENSION(neq) :: tmp

IF (inv == 0) THEN
   IF (tran == 1) THEN
      CALL dgetrs('T',neq,1,prec,neq,piv,vec,neq,info)
   ELSE
      CALL dgetrs('N',neq,1,prec,neq,piv,vec,neq,info)
   END IF
   IF (info /= 0) THEN
      WRITE (8,'(//,1X,A)') "ERROR: dgetrs cannot solve."
      STOP
   END IF
ELSE
   IF (tran == 1) THEN
      DO i = 1, neq
         tmp(i) = DOT_PRODUCT(prec(:,i),vec)
      END DO
      vec = tmp
   ELSE
      vec = MATMUL(prec,vec)
   END IF
END IF

RETURN
END SUBROUTINE precond
