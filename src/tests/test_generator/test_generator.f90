PROGRAM PIDOTS_TEST_GEN
  IMPLICIT NONE

  !Takes a list of desired tests and creates all test-specific input files files
  ! List must contain one desired test per line as
  ! <# processors per dimension>
  ! Code will build corresponding test for 2,4,8,16 subdomains per axis

  CHARACTER(200):: input_file=""
  INTEGER:: in_unit = 20, arg_count, ioerr
  INTEGER:: p_count

  arg_count = COMMAND_ARGUMENT_COUNT()

  IF (arg_count .EQ. 0) THEN
    WRITE(*,*) "Must provide input file name"
    STOP
  END IF

  CALL GET_COMMAND_ARGUMENT(1, input_file)
  WRITE(*,*) " Using test list: ", TRIM(input_file)

  OPEN(UNIT = in_unit, FILE = input_file, ACTION = "READ", STATUS = "OLD", IOSTAT = ioerr)
  IF (ioerr .NE. 0) THEN
    STOP 'Error opening input file'
  END IF

  DO WHILE (ioerr .EQ. 0)

    READ(in_unit,'(I4)',IOSTAT= ioerr) p_count
    IF (ioerr .NE. 0) CYCLE

    CALL BUILD_SRC(p_count)
    CALL BUILD_TEST(p_count, 2)
    CALL BUILD_TEST(p_count, 4)
    CALL BUILD_TEST(p_count, 8)
    CALL BUILD_TEST(p_count, 16)

  END DO

END PROGRAM

SUBROUTINE BUILD_SRC(p)
IMPLICIT NONE
  INTEGER, INTENT(IN):: p
  INTEGER:: ioerr
  CHARACTER(50):: src_filename

  WRITE(src_filename,'(A,I0)') "src",p

  OPEN(UNIT = 21, FILE = src_filename, ACTION = "WRITE", STATUS = "REPLACE", IOSTAT = ioerr)
  IF (ioerr .NE. 0) THEN
    STOP 'Error creating src file'
  END IF

  WRITE(21, '(A,I0)') "! Source data input file for p = ", p
  WRITE(21, '(I0)') 0
  WRITE(21, '(I0)') 1
  WRITE(21, '(3(I0,x))') 1, 1, 1
  WRITE(21, '(3(I0,x))') 16*p, 16*p, 16*p
  WRITE(21, '(F3.1)') 1.0

  CLOSE(21)

END SUBROUTINE


SUBROUTINE BUILD_TEST(p, s)
IMPLICIT NONE
  INTEGER, INTENT(IN):: p, s
  INTEGER:: ioerr
  CHARACTER(50):: test_filename, src_filename

  WRITE(src_filename,'(A,I0)') "src",p
  WRITE(test_filename,'(I0,A,I0,A)') p, "-", s, ".i"

  OPEN(UNIT = 21, FILE = test_filename, ACTION = "WRITE", STATUS = "REPLACE", IOSTAT = ioerr)
  IF (ioerr .NE. 0) THEN
    STOP 'Error creating src file'
  END IF

  WRITE(21, '(A,I0,A,I0)') "!Input file for case: ", p, "-", s
  WRITE(21, '(A)') "8 1                       Order & Qudrature Type: 0/1/2 => TWOTRAN/EQN/Read-in"
  WRITE(21, '(3(I0,x),A)') 16*p, 16*p, 16*p, "                       Number of Nodes in x- y- & z-Directions"
  WRITE(21, '(3(I0,x),A)') p, p, p, "                       Number of Processes in X Y and Z"
  WRITE(21, '(3(I0,x),A)') 16/s, 16/s, 16/s, "c                       Number of cells per sub-domain in X Y and Z"
  WRITE(21, '(A)') "1                         Number of Materials"
  WRITE(21,'(I0,A)') 16*p,"*1.0                       x-size of cells"
  WRITE(21,'(I0,A)') 16*p,"*1.0                       y-size of cells"
  WRITE(21,'(I0,A)') 16*p,"*1.0                       z-size of cells"
  WRITE(21, '(A)')"0 0                        +Z, -Z BCs"
  WRITE(21, '(A)')"0 0                        +Y, -Y BCs"
  WRITE(21, '(A)')"0 0                        +X, -X BCs"
  WRITE(21, '(A)')"mat                       Name of the material map file"
  WRITE(21, '(A)')"                          Name of the Quadrature input file"
  WRITE(21, '(A)')"xs                        Name of Cross Sections File"
  WRITE(21,'(A,A)') TRIM(src_filename),"                       Name of External Source File"
  WRITE(21,'(A)')"                          Name of BC Angular flux file"
  WRITE(21,'(A)')"1.e-6 1500 1.e-10         Conv.Crit., Max.Its., Tolerance"
  WRITE(21,'(A)')"1                         ITP"
  WRITE(21,'(A)')"1                         SFP"

END SUBROUTINE
