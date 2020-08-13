FUNCTION bnc(N,P,rank)

!-------------------------------------------------------------
!
!  Compute the number of cells for the subdomain, one
!   dimension at a time.
!  Take in the number of cells total for a dimension, N, the 
!   number of processes P in that dimension, and the rank.
!  Return the integer number of cells for that sub-domain, bnc.
!
!  When not evenly divided, keeps all but one even, then tries
!   tries to keep last sub-domain close
!
!-------------------------------------------------------------

IMPLICIT NONE
INTEGER, INTENT(IN) :: N, P, rank
INTEGER :: bnc, tmpi1, tmpi2
REAL :: tmpr

! Set the values to be checked against
tmpr = 0.5*REAL(P)

tmpi1 = MOD(N,P)
tmpi2 = N/P

! If it's evenly divisible, just use the division
IF (tmpi1 == 0) THEN
   bnc = tmpi2

! If not, then compare the MOD result to half the number of processes and set
ELSE IF (REAL(tmpi1) < tmpr) THEN
   IF (rank /= P-1) THEN
      bnc = tmpi2
   ELSE
      bnc = N - (P-1)*tmpi2
   END IF

ELSE IF (REAL(tmpi1) >= tmpr) THEN
   IF (rank /= P-1) THEN
      bnc = tmpi2 + 1
   ELSE
      bnc = N - (P-1)*(tmpi2+1)
   END IF

END IF

END FUNCTION bnc
