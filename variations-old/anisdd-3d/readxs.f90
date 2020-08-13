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
!   sigs(m,g,g') => scattering cross section from g' to g
!
!-------------------------------------------------------------

USE invar
IMPLICIT NONE
INTEGER, INTENT(INOUT) :: iexit
INTEGER :: m, g, gp, l
CHARACTER(8), INTENT(IN) :: xsfile
REAL*8 :: temp

! Set up the size of the arrays for the cross sections
ALLOCATE(sigt(nm,ng), sigs(0:anord,nm,ng,ng), ssum(nm,ng))

! Open the cross-section file for reading
OPEN (UNIT = 11, FILE=xsfile)

READ (11, *)
! Loop over material overall
DO m = 1, nm
   READ (11,*)
   ! Next loop is the INTO group
   DO g = 1, ng
      READ(11,*)
      READ(11,*) sigt(m,g)
      ! Separate multigroup scattering by anisotropic order
      DO l = 0, anord
         ! Read by scattering to other groups from group g
         READ(11,*) (sigs(l,m,gp,g), gp = 1, ng)
      END DO
   END DO
END DO

! Perform checks on the cross section data
IF (MINVAL(sigt) < 0.) THEN
   WRITE(8,'(/,3X,A)') "ERROR: total cross section must be zero or greater"
   iexit = iexit + 1
ELSE IF (MINVAL(sigs) < 0.) THEN
   WRITE(8,'(/,3X,A)') "ERROR: scattering cross section must be zero or greater"
   iexit = iexit + 1
END IF

ssum = 0.0
DO m = 1, nm
   DO gp = 1, ng
      DO g = 1, ng
         ssum(m,gp) = ssum(m,gp) + sigs(0,m,g,gp)
      END DO
      IF (ssum(m,gp) > sigt(m,gp)) THEN
         WRITE(8,'(/,3X,A)') "ERROR: Scattering XS must be less than total XS"
         iexit = iexit + 1
      END IF
   END DO
END DO

! Check the anisotropic orders
IF (anord > 0) THEN
   DO m = 1, nm
      DO g = 1, ng
         DO gp = 1, ng
            temp = 0.0
            DO l = 1, anord
               temp = temp + sigs(l,m,gp,g)
            END DO
            IF (temp > sigs(0,m,gp,g)) THEN
               WRITE(8,'(/,3X,A)') "ERROR: Anisotropic order scattering XSs must sum to less than 0-order."
               iexit = iexit + 1
            END IF
         END DO
      END DO
   END DO
END IF

RETURN
END SUBROUTINE readxs
