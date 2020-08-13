SUBROUTINE stages(iexit)

!-------------------------------------------------------------
!
! Using number of nodes input, determine the number of
!  V-Cycle stages and the coarsening factor (1,2,4) for each
!  stage in each direction.
!
!-------------------------------------------------------------

USE invar
IMPLICIT NONE
INTEGER, INTENT(INOUT) :: iexit
INTEGER :: i, t, tt

! First compute the total number of stages, ns
t = MAXVAL([npx,npy,npz])
i = 1
tt = 1
DO
   IF (i < t) THEN
      i = 2*i
      tt = tt + 1
      CYCLE
   ELSE IF (i == t) THEN
      EXIT
   ELSE
      WRITE(8,'(/,3X,A)') "ERROR: Pi/Pj/Pk must be a power of 2."
      iexit = iexit + 1
      EXIT
   END IF
END DO
ns = 1 + tt/2

! Continue by determining the coarsening factors in each direction at each stage
! (if no error from before)
IF (iexit == 0) THEN
   ! Allocate the arrays that store the coarsening factors
   ALLOCATE(xcf(ns-1),ycf(ns-1),zcf(ns-1))

   ! x-direction
   t = npx
   DO i = 1, ns
      ! Coarsen by 4 if evenly divisible by 4
      IF (MOD(t,4) == 0) THEN
         xcf(i) = 4
         t = t/4
      ! Coarsen by 2 if evenly divisible by 2
      ELSE IF (MOD(t,2) == 0) THEN
         xcf(i) = 2
         t = t/2
      ELSE
         ! No need to coarse if already down to P=1
         IF (t == 1) THEN
            xcf(i) = 1
         ! Exit loop if P is not a power of 2
         ELSE
            WRITE(8,'(/,3X,A)') "ERROR: Pi must be a power of 2."
            iexit = iexit + 1
            EXIT
         END IF
      END IF
   END DO

   ! y-direction
   t = npy
   DO i = 1, ns
      ! Coarsen by 4 if evenly divisible by 4
      IF (MOD(t,4) == 0) THEN
         ycf(i) = 4
         t = t/4
      ! Coarsen by 2 if evenly divisible by 2
      ELSE IF (MOD(t,2) == 0) THEN
         ycf(i) = 2
         t = t/2
      ELSE
         ! No need to coarse if already down to P=1
         IF (t == 1) THEN
            ycf(i) = 1
         ! Exit loop if P is not a power of 2
         ELSE
            WRITE(8,'(/,3X,A)') "ERROR: Pj must be a power of 2."
            iexit = iexit + 1
            EXIT
         END IF
      END IF
   END DO

   ! z-direction
   t = npz
   DO i = 1, ns
      ! Coarsen by 4 if evenly divisible by 4
      IF (MOD(t,4) == 0) THEN
         zcf(i) = 4
         t = t/4
      ! Coarsen by 2 if evenly divisible by 2
      ELSE IF (MOD(t,2) == 0) THEN
         zcf(i) = 2
         t = t/2
      ELSE
         ! No need to coarse if already down to P=1
         IF (t == 1) THEN
            zcf(i) = 1
         ! Exit loop if P is not a power of 2
         ELSE
            WRITE(8,'(/,3X,A)') "ERROR: Pk must be a power of 2."
            iexit = iexit + 1
            EXIT
         END IF
      END IF
   END DO
END IF

RETURN
END SUBROUTINE stages
