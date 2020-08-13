SUBROUTINE afcm(i,xs,xe,incx,n,oct,mu)

!-------------------------------------------------------------
!
!  Angular Flux Coefficient Matrix
!  Directs the construction of the matrices for working with
!   the incoming angular flux. Constructs kmat's and kpsi's
!
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER, INTENT(IN) :: i, xs, xe, incx, n, oct
INTEGER :: ieq
REAL*8, INTENT(IN) :: mu
REAL*8, DIMENSION(sord) :: ktmp

! Save the BC matrices for updates
IF (i /= xs) THEN
   xbco = xbc
END IF

! Update BC matrices with xbco
IF (i /= xs) THEN
   xbc = gxx*xbco
END IF

! Adjust the values for BCs to their adjacent cells
IF (i == xs) THEN
   xbc = gxx
END IF

! Update the appropriate kmat given the octant of interest
! Call conk to construct the row of kmat, then set it in kmat#
! Set according to 'tpose' flag
IF (tpose == 1) THEN
   ieq = (i-1)*sord + 1
   CALL conk(n,i,xs,incx,mu,ktmp)
   IF (oct == 1) THEN
      kmat(n,ieq:ieq+anord,1) = kmat(n,ieq:ieq+anord,1) + ktmp
   ELSE IF (oct == 2) THEN
      kmat(n,ieq:ieq+anord,2) = kmat(n,ieq:ieq+anord,2) + ktmp
   END IF
   
   ! Use the data in the xbc variable at xe to compute kpsi
   IF (i == xe) THEN
      kpsi(n,oct) = xbc
   END IF
!-------------------------------------------------------------------------------
ELSE
   ieq = (i-1)*sord + 1
   CALL conk(n,i,xs,incx,mu,ktmp)
   IF (oct == 1) THEN
      kmat(ieq:ieq+anord,n,1) = kmat(ieq:ieq+anord,n,1) + ktmp
   ELSE IF (oct == 2) THEN
      kmat(ieq:ieq+anord,n,2) = kmat(ieq:ieq+anord,n,2) + ktmp
   END IF

   ! Use the data in the xbc variable at xe to compute kpsi
   IF (i == xe) THEN
      kpsi(n,oct) = xbc
   END IF
END IF

RETURN
END SUBROUTINE afcm
