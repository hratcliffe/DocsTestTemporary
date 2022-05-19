!> Handle user interaction (command line etc)
!> \author Fred Blogs \date 27/1/2020
MODULE user_interaction

  USE kinds
  USE command_line

  IMPLICIT NONE

  !> Unit number to use for temporary files
  INTEGER, PARAMETER :: file_unit = 105

 CONTAINS

!> Print help for command line arguments
!> Help comes from a text file (see Main Page)
SUBROUTINE print_help

  INTEGER :: iost
  CHARACTER(LEN=100) :: line ! Hopefully this is long enough!

  iost = 0

  OPEN(UNIT=file_unit, FILE='input_readme.txt', IOSTAT=iost)

  DO
    IF(iost /= 0) EXIT
    READ(file_unit, '(A)', IOSTAT=iost) line
    PRINT*, TRIM(ADJUSTL(line))
  END DO

  CLOSE(file_unit)

end SUBROUTINE

!> Process and extract the command line arguments
!> \param N To be filled with N (grid size)
!> \param max_iters  To be filled with max_iters (timesteps)
!> \param init To be filled with init (initial problem spec)
!> \param beta To be filled with beta(inverse temperature)
!> \param J To be filled with J
!> \param err Returns 0 if all mandatory arguments are present, 1 else
!> \ext Beta and J need not be essential
!>
!> \author Heather Ratcliffe
SUBROUTINE fill_command_args(N, max_iters, init, beta, J, err)

  INTEGER(KIND=INT32), INTENT(OUT) :: N, max_iters, err
  REAL(KIND=REAL64), INTENT(OUT) :: beta, J
  CHARACTER, INTENT(INOUT) :: init
  LOGICAL :: success, all_success, help_printed
  CHARACTER(LEN=30) :: help

  ! Get command line arguments
  CALL parse_args
  ! ALL parameters are required!!

  help_printed = .FALSE.

  success = get_arg("help", help)
  IF(TRIM(ADJUSTL(help)) /= "F") THEN
    CALL print_help
    help_printed = .TRUE.
  END IF

  ! Best and simplest way to use cmd line module
  success = get_arg("max_iters", max_iters)
  all_success = .TRUE. .AND. success

  ! Repeat for N
  success = get_arg("N", N)
  all_success = all_success .AND. success

  success = get_arg("beta", beta)
  all_success = all_success .AND. success

  success = get_arg("J", J)
  all_success = all_success .AND. success

  ! I use a single-character flag for the initial state,
  ! either 'R' for random, 'A' for alternating or 'F' for flat (=1)
  ! and I default to random if not given or not one of the options
  success = get_arg("init", init)

  all_success = all_success .AND. success


  IF( .NOT. all_success) THEN

    IF( .NOT. help_printed) CALL print_help
    err = 1

  END IF

END SUBROUTINE

END MODULE
