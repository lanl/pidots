SUBROUTINE jima(i,xs,xe,incx,n,oct,mu,tsxs)

!-------------------------------------------------------------
!
!  Jacobian Iteration MAtrix
!  Directs the equations to set up the matrices that work
!   with the phi and source vectors. takes in parameters from
!   matsweep and solves for the updates to jmat and jpsi 
!
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER, INTENT(IN) :: i, xs, xe, incx, n, oct
INTEGER :: ii, l, ieq, jeq
REAL*8, INTENT(IN) :: mu
REAL*8, DIMENSION(0:anord), INTENT(IN) :: tsxs
REAL*8 :: wt, plm1, pl

! Set the quad weight
wt = w(n)

! Save the X matrix for updates
IF (i /= xs) THEN
   xold = xmat
END IF

! Update X matrix with old X matrix
IF (i /= xs) THEN
   DO ii = xs, i, incx
      xmat(:,ii) = gxx*xold(:,ii)
   END DO
END IF

! Append X matrix
DO l = 0, anord
   xmat(l,i) = gxa(l)*tsxs(l)
END DO

! Start to formulate Jacobian matrix
! Depends on the 'tpose' flag how the indices are incremented
! tpose = 1 --> make transposed Jphi and Jpsi
IF (tpose == 1) THEN
   jeq = (i-1)*sord + 1
   ! Diagonal Block
   plm1 = 0.0
   pl   = 1.0
   DO l = 0, anord
      jmat(jeq:jeq+anord,jeq+l) = jmat(jeq:jeq+anord,jeq+l) + wt*pl*gaa*tsxs
      CALL legpoly(l,incx,mu,plm1,pl)
   END DO

   ! Off-diagonal contribution from X
   IF (i /= xs) THEN
      DO ii = xs, i, incx
         ieq = (ii-1)*sord + 1
         plm1 = 0.0
         pl   = 1.0
         DO l = 0, anord
            jmat(ieq:ieq+anord,jeq+l) = jmat(ieq:ieq+anord,jeq+l) + wt*pl*gax*xold(:,ii)
            CALL legpoly(l,incx,mu,plm1,pl)
         END DO
      END DO
   END IF

   ! Now check to see if the computational cell is at end
   ! Place an xmat set of values if it is
   IF (i == xe) THEN
      DO ii = xs, i, incx
         ieq = (ii-1)*sord + 1
         jpsi(ieq:ieq+anord,n,oct) = xmat(:,ii)
      END DO
   END IF
!-----------------------------------------------------------------------
! If tpose = 0 --> make normal Jphi and Jpsi
ELSE
   ieq = (i-1)*sord + 1
   ! Diagonal Block
   plm1 = 0.0
   pl   = 1.0
   DO l = 0, anord
!      DO k = 0, anord
!         jmat(ieq+l,ieq+k) = jmat(ieq+l,ieq+k) + wt*pl*gaa(k)*tsxs(k)
!      END DO
      jmat(ieq+l,ieq:ieq+anord) = jmat(ieq+l,ieq:ieq+anord) + wt*pl*gaa*tsxs
      CALL legpoly(l,incx,mu,plm1,pl)
   END DO

   ! Off-diagonal contribution from X
   IF (i /= xs) THEN
      DO ii = xs, i, incx
         jeq = (ii-1)*sord + 1
         plm1 = 0.0
         pl   = 1.0
         DO l = 0, anord
!            DO k = 0, anord
!               jmat(ieq+l,jeq+k) = jmat(ieq+l,jeq+k) + wt*pl*gax*xold(k,ii)
!            END DO
            jmat(ieq+l,jeq:jeq+anord) = jmat(ieq+l,jeq:jeq+anord) + wt*pl*gax*xold(:,ii)
            CALL legpoly(l,incx,mu,plm1,pl)
         END DO
      END DO
   END IF

   ! Now check to see if the computational cell is at end
   ! Place an xmat set of values if it is
   IF (i == xe) THEN
      DO ii = xs, i, incx
         jeq = (ii-1)*sord + 1
         jpsi(n,jeq:jeq+anord,oct) = xmat(:,ii)
      END DO
   END IF
END IF

RETURN
END SUBROUTINE jima
