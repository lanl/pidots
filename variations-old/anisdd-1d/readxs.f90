SUBROUTINE readxs(xsfile,iexit)

!-------------------------------------------------------------
!
! Reads the cross sections from a file
!  Cross sections read by group and by material number
!  Read in the total cross section followed by the
!   full scattering matrix.
!  Limit the scattering to down-scatter only for now
!  
!   sigt(m,g) => total cross section of material m, group g
!   sigs(m,l,g,g') => l-th order scattering cross section from g' to g
!
!-------------------------------------------------------------

USE invar
IMPLICIT NONE
INTEGER, INTENT(INOUT) :: iexit
INTEGER :: m, g, gp, l
CHARACTER(8), INTENT(IN) :: xsfile
REAL*8 :: temp

! Set up the size of the arrays for the cross sections
ALLOCATE(sigt(nm,ng), sigs(nm,0:anord,ng,ng), ssum(nm,ng))

! Open the cross-section file for reading
OPEN (UNIT = 11, FILE=xsfile)

READ (11, *)
! Loop over material overall
DO m = 1, nm
   READ (11,*)
   ! Next loop is the FROM group
   DO g = 1, ng
      READ(11,*)
      READ(11,*) sigt(m,g)
      ! Separate multigroup scattering by anisotropic order
      DO l = 0, anord
         ! Read by scattering to other groups from group g
         READ(11,*) (sigs(m,l,gp,g), gp = 1, ng)
      END DO
   END DO

END DO

! Perform checks on the cross section data
IF (MINVAL(sigt) < 0.) THEN
   WRITE(8,'(/,3X,A)') "ERROR: total cross section must be zero or greater"
   iexit = iexit + 1
ELSE IF (MINVAL(sigs) < 0.) THEN
   WRITE(8,'(/,3X,A)') "ERROR: scattering cross sections must be zero or greater"
   iexit = iexit + 1
END IF

DO m = 1, nm
   DO g = 1, ng
      DO gp = 1, ng
         ssum(m,g) = ssum(m,g) + sigs(m,0,gp,g)
      END DO
      IF (ssum(m,g) > sigt(m,g)) THEN
         WRITE(8,'(/,3X,A)') "ERROR: Scattering XS must be less than total XS"
         iexit = iexit + 1
      END IF
   END DO
END DO

! Check the anistropic orders
IF (anord > 0) THEN
   DO m = 1, nm
      DO g = 1, ng
         DO gp = 1, ng
            temp = 0.0
            DO l = 1, anord
               temp = temp + sigs(m,l,gp,g)
            END DO
            IF (temp > sigs(m,0,gp,g)) THEN
               WRITE(8,'(/,3X,A)') "ERROR: Anisotropic order scattering XSs must sum to less than 0-order."
               iexit = iexit + 1
            END IF
         END DO
      END DO
   END DO
END IF

RETURN
END SUBROUTINE readxs
