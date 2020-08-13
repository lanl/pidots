SUBROUTINE angle(iexit)

!-------------------------------------------------------------
!
! Prepares the angular quadrature with directions from input
!
!-------------------------------------------------------------

USE invar
IMPLICIT NONE
INTEGER, INTENT(INOUT) :: iexit 
REAL*8 :: wtsum
INTEGER :: n

! Check that the type is now 0 -- 1 has been ruled out already
IF (qdtyp /= 0) THEN
   WRITE(8,'(/,3x,A)') "ERROR: Illegal value for the quadrature type."
   iexit = iexit + 1
END IF

! If qdtyp is 0, use the Legendre Quadrature Set for Slab Geometry (From Lewis & Miller)
! Set the angles (cosines)
IF (qdord == 2) THEN
   ang(1) = 0.5773502691
ELSE IF (qdord == 4) THEN
   ang(1) = 0.3399810435
   ang(2) = 0.8611363115
ELSE IF (qdord == 6) THEN
   ang(1) = 0.2386191860
   ang(2) = 0.6612093864
   ang(3) = 0.9324695142
ELSE IF (qdord == 8) THEN
   ang(1) = 0.1834346424
   ang(2) = 0.5255324099
   ang(3) = 0.7966664774
   ang(4) = 0.9602898564
ELSE IF (qdord == 10) THEN
   ang(1) = 0.1488743389
   ang(2) = 0.4333953941
   ang(3) = 0.6794095682
   ang(4) = 0.86506336663
   ang(5) = 0.9739065285
ELSE IF (qdord == 12) THEN
   ang(1) = 0.1252334085
   ang(2) = 0.3678314989
   ang(3) = 0.5873179542
   ang(4) = 0.7699026741
   ang(5) = 0.9041172563
   ang(6) = 0.9815606342
ELSE
   WRITE(8,'(/,3x,A)') "ERROR: Illegal value for LQn qdord (N). Must be 2, 4, 6, 8, 10, or 12."
   iexit = iexit + 1  
END IF

! Set the weights
IF (qdord == 2) THEN
   w(1) = 1.0000000000
ELSE IF (qdord == 4) THEN
   w(1) = 0.6521451549
   w(2) = 0.3478548451
ELSE IF (qdord == 6) THEN
   w(1) = 0.4679139346
   w(2) = 0.3607615730
   w(3) = 0.1713244924
ELSE IF (qdord == 8) THEN
   w(1) = 0.3626837834
   w(2) = 0.3137066459
   w(3) = 0.2223810344
   w(4) = 0.1012285363
ELSE IF (qdord == 10) THEN
   w(1) = 0.2955242247
   w(2) = 0.2692667193
   w(3) = 0.2190863625
   w(4) = 0.1494513492
   w(5) = 0.0666713443
ELSE IF (qdord == 12) THEN
   w(1) = 0.2491470458
   w(2) = 0.2334925365
   w(3) = 0.2031674267
   w(4) = 0.1600783286
   w(5) = 0.1069393260
   w(6) = 0.0471753364
END IF

! Renormalize all the weights
wtsum = SUM(w)
DO n = 1, apo
   w(n) = w(n) * 0.5/wtsum
END DO

RETURN
END SUBROUTINE angle
