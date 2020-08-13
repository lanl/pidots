FUNCTION bnc(N,P,rank)

!-------------------------------------------------------------
!
!  Compute the number of cells for the subdomain, one
!   dimension at a time.
!  Take in the number of cells total for a dimension, N, the 
!   number of processes P in that dimension, and the rank.
!  Return the integer number of cells for that sub-domain, bnc.
!
!  When not evenly divided, keeps all but one even, then puts
!   remainder in with final sub-domain
!
!-------------------------------------------------------------

USE invar
IMPLICIT NONE
INTEGER, INTENT(IN) :: N, P, rank
INTEGER :: bnc, tmpi1, tmpi2

tmpi1 = MOD(N,P)
tmpi2 = N/P

! If it's evenly divisible, just use the division
IF (tmpi1 == 0) THEN
   bnc = tmpi2

! If not, then set all but the last process equal to division
! Last process set to remainder
ELSE
   IF (rank == 0) THEN
      warn = warn + 1
      WRITE (8,'(/,A,/)') "WARNING: Number of cells not evenly divided by processes."
   END IF
   IF (rank < tmpi1) THEN
      bnc = tmpi2 + 1
   ELSE
      bnc = tmpi2
   END IF
END IF

END FUNCTION bnc
