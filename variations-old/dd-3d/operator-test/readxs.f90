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
INTEGER :: m
CHARACTER(8), INTENT(IN) :: xsfile

! Set up the size of the arrays for the cross sections
ALLOCATE(sigt(nm), sigs(nm))

! Open the cross-section file for reading
OPEN (UNIT = 11, FILE=xsfile)

READ (11, *)
! Loop over material overall
DO m = 1, nm
   READ(11,*)
   READ(11,*) sigt(m)
   READ(11,*) sigs(m)
END DO

! Perform checks on the cross section data
IF (MINVAL(sigt) < 0.) THEN
   WRITE(8,'(/,3X,A)') "ERROR: total cross section must be zero or greater"
   iexit = iexit + 1
ELSE IF (MINVAL(sigs) < 0.) THEN
   WRITE(8,'(/,3X,A)') "ERROR: scattering cross section must be zero or greater"
   iexit = iexit + 1
END IF

DO m = 1, nm
   IF (sigs(m) > sigt(m)) THEN
      WRITE(8,'(/,3X,A)') "ERROR: Scattering XS must be less than total XS"
      iexit = iexit + 1
   END IF
END DO

RETURN
END SUBROUTINE readxs
