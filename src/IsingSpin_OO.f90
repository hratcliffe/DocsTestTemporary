! DOXYHINT - The main page of docs can be filled like this
!> \mainpage
!> Simple Ising Spin model using Object Oriented Fortran
!> \section Intro Introduction
!> A basic program for magnetic domain simulations using the
!> Ising Model.
!> Based on my model solutions for PX913 and extended into a
!> terribly overbuilt Object Oriented version. This relies
!> on the helper modules kinds.f90, command_line.f90, ascii_display.f90,
!> and also random_mod.f90 and sleep_mod.f90
!> HELLO I HAVE CHANGED MY MAINPAGE
!> \section Inputs Inputs
!> \verbinclude input_readme.txt
!> \author Heather Ratcliffe \date 27/1/2020


!This code is overbuilt for this problem, but shows some nice
! fortran stuff
! Note that this uses extra memory in places to get
! nice elegant operations
! I have deliberately kept the core code as in the simple
! Ising and just wrapped things away.

!DOXYHINT Documents the module which immediately follows
! Note the arrow points forward in Fortran style doxygen
!> The core Ising types and functions
!> \author Heather Ratcliffe \date 1972
MODULE ising_functions

  USE random_mod
  USE kinds

  IMPLICIT NONE

  ! No array-of-pointers in Fortran
  !> Wraps a pointer so we can use an array of them
  TYPE :: ptr_wrapper
    INTEGER(KIND=INT16), POINTER :: ptr
  END TYPE ptr_wrapper


  ! This is VERY overbuilt... But it shows what is possible
  !> A single spin element
  !DOXYHINT End-of line docs use the left-angle i.e. the arrow points to the thing you're documenting
  TYPE :: ising_element
    ! Will point into the grids 'state' array
    ! DOXYHINT Doxygen will try and create a link if you mention another type, function etc, like ising_grid here
    INTEGER(KIND=INT16), POINTER :: state !< Pointer to the spin state of the element in ising_grid
    TYPE(ptr_wrapper), DIMENSION(:), ALLOCATABLE :: neighbours !< Pointers to elements neighbours
    INTEGER :: el_flips=0 !< count of flips of this element
    CONTAINS
    PROCEDURE :: flip !< Flip the spin of this element
    PROCEDURE :: magnetization_el !< Calculate and return the magnetiation of this element
    ! Suppose elements could have different ways to calculate magnetization
  END TYPE ising_element

  ! Since I display the grid etc, I have a plain array
  ! of the current state, as well as the array of elements
  ! The elements can track how
  ! often they are flipped, etc
  !> Holds the grid of elements and provides functions on it
  TYPE :: ising_grid
    INTEGER :: sz=0 !< Size of grid (for convenience)
    INTEGER(KIND=INT16), DIMENSION(:, :), POINTER :: state !<State array
    !DOXYHINT Sometimes you need/want to qualify the type you're mentioning. Use :: to explain the "scope" 
    TYPE(ising_element), DIMENSION(:,:), POINTER :: elements !< Element array. Each ising_element has a pointer into the ising_grid::state array
    REAL(KIND=REAL64) :: J=1 !< Ising model parameters - J
    REAL(KIND=REAL64) :: beta=1 !< Ising model parameters - beta (inverse temperature)
    CONTAINS
    !DOXYHINT Fortran OO is not (yet?) supported, but we can get a good approximation using an ALIAS
    !DOXYHINT which copies the docs from the relevant function. We just specify which that is
    !DOXYHINT and document the function in the source below as normal
    PROCEDURE :: allocate_grid !< \procdoc{ising_functions::allocate_grid}
    PROCEDURE :: iterate !< \procdoc{ising_functions::iterate}
    PROCEDURE :: setup_grid !< \procdoc{ising_functions::setup_grid}
    PROCEDURE :: magnetization !< \procdoc{ising_functions::magnetization}
    PROCEDURE :: flips !< \procdoc{ising_functions::flips}
    FINAL :: destroy_grid !< \procdoc{ising_functions::destroy_grid}
  END TYPE ising_grid

  CONTAINS

  ! Function(s) for individual elements
  !> Flip the given element
  !> \param this Element to flip
  SUBROUTINE flip(this)
    CLASS(ising_element) :: this
    this%state = - this%state
    this%el_flips = this%el_flips + 1
  END SUBROUTINE


  ! Functions for the grid
  !> Allocate the grid and setup neighbours
  !> \param this The grid to allocate
  !> \param N Size of grid
  !> \caveat Does 4 neighbours in a + formation, but at edges/corners the
  !> missing elements are omitted
  !>\see Derivations PDF \cite Derivations
  SUBROUTINE allocate_grid(this, N)
    CLASS(ising_grid) :: this
    INTEGER :: N, i, j

    this%sz = N
    !> \todo Example of an in-body todo!
    ! This time, I have the sized domain and limit the neighbours
    ALLOCATE(this%state(1:N, 1:N))
    ALLOCATE(this%elements(1:N, 1:N))

    ! Elements aren't ready yet!!


    ! All of the complexity is now here in the setup, and the iterate routine
    ! Doesn't need to know how we handle boundaries, or how many neighbours we have
    DO i = 1, N
      DO j = 1, N
        ! Do 4 neighbours in a + config. Could do other things
        ! But at the boundaries, only include the cells which exist
        ! We could also do periodic then all cells have exactly 4
        ! Note we could include more distance neignbours, diagonals etc
        this%elements(i,j)%state => this%state(i,j)

        IF ( (i == 1 .OR. i ==N) .AND. (j ==1 .OR. j == N)) THEN
          ! Corners:
          ALLOCATE(this%elements(i,j)%neighbours(2))
          IF(i == 1) THEN
            this%elements(i,j)%neighbours(1)%ptr => this%state(2, j)
          ELSE IF (i == N) THEN
            this%elements(i,j)%neighbours(1)%ptr => this%state(N-1, j)
          END IF
          IF (j == 1) THEN
            this%elements(i,j)%neighbours(2)%ptr => this%state(i, 2)
          ELSE IF (j == N) THEN
            this%elements(i,j)%neighbours(2)%ptr => this%state(i, N-1)
          END IF

        ELSE IF(i == 1 .OR. i ==N) THEN
          ! Non-corner Edges - i
          ALLOCATE(this%elements(i,j)%neighbours(3))
          IF(i == 1) THEN
            this%elements(i,j)%neighbours(1)%ptr => this%state(2, j)
          ELSE IF (i == N) THEN
            this%elements(i,j)%neighbours(1)%ptr => this%state(N-1, j)
          END IF
          this%elements(i,j)%neighbours(2)%ptr => this%state(i, j-1)
          this%elements(i,j)%neighbours(3)%ptr => this%state(i, j+1)

        ELSE IF(j ==1 .OR. j == N) THEN
          ! Non-corner Edges - j
          ALLOCATE(this%elements(i,j)%neighbours(3))
          IF(j == 1) THEN
            this%elements(i,j)%neighbours(1)%ptr => this%state(i, 2)
          ELSE IF (j == N) THEN
            this%elements(i,j)%neighbours(1)%ptr => this%state(i,N-1)
          END IF
          this%elements(i,j)%neighbours(2)%ptr => this%state(i-1, j)
          this%elements(i,j)%neighbours(3)%ptr => this%state(i+1, j)

        ELSE
          ALLOCATE(this%elements(i,j)%neighbours(4))
          this%elements(i,j)%neighbours(1)%ptr => this%state(i-1, j)
          this%elements(i,j)%neighbours(2)%ptr => this%state(i+1, j)
          this%elements(i,j)%neighbours(3)%ptr => this%state(i  , j-1)
          this%elements(i,j)%neighbours(4)%ptr => this%state(i  , j+1)
        END IF
      END DO
    END DO

  END SUBROUTINE allocate_grid

  !> Finalizer - called when grid is destroyed
  SUBROUTINE destroy_grid(this)
    ! For finalizer, this arg MUST be TYPE (not CLASS)
    TYPE(ising_grid) :: this

    DEALLOCATE(this%state, this%elements)
    this%sz = 0

  END SUBROUTINE


  ! This function just wraps around the original setup functions
  ! passing the grid through
  !> Setup the grid. \see setup_alternating_spin_grid and setup_random_spin_grid
  SUBROUTINE setup_grid(this, N,  init, beta, J)
    CLASS(ising_grid) :: this
    CHARACTER(LEN=1) :: init
    INTEGER :: N
    REAL(KIND=REAL64) :: beta, J

    CALL allocate_grid(this, N)

    this%beta = beta
    this%J = J

    IF(init == "F" .OR. init == "f") THEN
      ! We can set everything up '+1'
      this%state = 1
    ELSE IF(init == "A" .OR. init == "a") THEN
      ! This will set alternating
      CALL setup_alternating_spin_grid(this%state)
    ELSE
      ! This will set random spins
      CALL setup_random_spin_grid(this%state)
    END IF

  END SUBROUTINE setup_grid

  ! \todo Clean up refs to assignments and other code
  ! This is not changed a lot from the iterate function in IsingSpin
  !> Iterate the grid
  !> \param this The grid object
  !> \param max_iters The number of iterations to do
  !> \author Heather \date 25/12/19
  SUBROUTINE iterate(this, max_iters)

    CLASS(ising_grid) :: this
    TYPE(ising_element), POINTER :: test_element
    INTEGER(KIND=INT32) :: iter, max_iters, n_neighbour
    REAL(KIND=REAL64) :: rand, prob, delta_E
    INTEGER :: rand_x, rand_y

    !Iterate random flips
    DO iter = 1, max_iters

      ! Pick random pos'n
      rand_x = 1+INT(random()*(this%sz))
      rand_y = 1+INT(random()*(this%sz))

     ! For clarity of the rest of the code, use a temporary pointer
      test_element => this%elements(rand_x, rand_y)

      ! Calc energy change
      delta_E = 0.0
      ! Now we just sum over the neighbours, without worrying where they are
      DO n_neighbour = 1, SIZE(test_element%neighbours)
        delta_E = delta_E + this%J * test_element%neighbours(n_neighbour)%ptr * &
            test_element%state
      END DO

      ! Check for acceptance and perform flip if necessary
      IF (delta_E .LT. 0) THEN
        CALL test_element%flip()
      ELSE
        rand = random()
        prob = EXP(- this%beta * delta_E)
        IF (rand .LT. prob) THEN
          CALL test_element%flip()
        END IF
      END IF

    END DO

  END SUBROUTINE

  !> Calculate total flips
  FUNCTION flips(this)

    CLASS(ising_grid) :: this
    INTEGER :: flips

    flips = SUM(this%elements%el_flips)

  END FUNCTION flips



  !> Calculate magnetization
  FUNCTION magnetization(this)

    CLASS(ising_grid) :: this
    REAL(KIND=REAL64) :: magnetization

    magnetization = SUM(this%elements%magnetization_el())


  END FUNCTION magnetization

  !> Calculate magnetization for a single element
  ! ELEMENTAL means the function can be applied element-wise to an array
  !> \param[in] this Element to calculate from
  !> \ext Extend this to Joe Bloggs methods from \cite Bloggs2018
  PURE ELEMENTAL FUNCTION magnetization_el(this)

    ! We have to specify intent for elemental function
    CLASS(ising_element), INTENT(IN) :: this
    REAL(KIND=REAL64) :: magnetization_el

    ! Suppose we could have a more complex calculation here
    magnetization_el = this%state

  END FUNCTION magnetization_el



!-----------------------These are the original functions from IsingSpin ---------------
  !DOXYHINT A Doxygen group collects multiple functions, variables etc together. You don't need
  ! to use groups explicitly when you have modules.But for developer docs, it can be helpful to
  ! group more closely or broadly than your module structure
  !DOXYHINT the @{ wraps everything until @} NOTE the position of @ as part of the group we just defined
  !> \defgroup core Core Functions
  !> \brief The basic functions from the Ising model
  !>
  !> These are the core functions operating on basic data type. DO NOT call these
  !> except from the Ising spin member functions!
  !> @{
  ! Blank Doygen lines close environments
  !> Setup the core spin grid with random states
  !> \param[inout] grid The (Integer) grid to populate
  SUBROUTINE setup_random_spin_grid(grid)

    INTEGER(KIND=INT16), DIMENSION(:,:), INTENT(INOUT) :: grid
    INTEGER, DIMENSION(2) :: sz
    INTEGER(KIND=INT32) :: i, j
    REAL(KIND=REAL64) :: rand

    sz = SHAPE(grid)
    grid = 0

    ! Setup random spins
    DO j = 1, sz(2)
      DO i = 1, sz(1)
        rand = random()*2.0
        ! note our choice of Int16 gets a bit annoying because
        ! literals default to INT32 so we have to add explicit conversions
        grid(i, j) = INT(2*INT(rand) - 1, KIND=INT16)
      END DO
    END DO

  END SUBROUTINE setup_random_spin_grid


  !> Setup the core spin grid with alternating (checkerboard) states
  !> \param[inout] grid The (Integer) grid to populate
  SUBROUTINE setup_alternating_spin_grid(grid)

    INTEGER(KIND=INT16), DIMENSION(:,:), INTENT(INOUT) :: grid
    INTEGER, DIMENSION(2) :: sz
    INTEGER(KIND=INT32) :: i, j, offset

    sz = SHAPE(grid)

    ! Setup alternating spins:
    ! Set all to 1, then reverse every second
    grid = 1
    DO j = 1, sz(2)
      ! This makes sure rows alternate first and second cell flipped
      offset = 0
      IF(MOD(j, 2) .EQ. 1) offset = 1
      DO i = offset+1, sz(1), 2
        grid(i, j) = -1
      END DO
    END DO

  END SUBROUTINE setup_alternating_spin_grid

!> @}

END MODULE ising_functions


!> Main Program
!> \todo Output to file!
! While Doxygen can extract comments from inside FUNCTIONs/SUBROUTINEs
! it doesn't extract them in Main in my version, so I have to put todos
! here. But in functions I can put them where they fit.
PROGRAM main

  !Provided pretty printing module
  USE ascii_display
  USE sleep_mod
  USE user_interaction

  USE ising_functions

  IMPLICIT NONE

  INTEGER(KIND=INT32):: N
  INTEGER(KIND=INT32):: max_iters
  INTEGER(KIND=INT32) :: total_flips, err
  REAL(KIND=REAL64) :: beta, J, init_mag, final_mag
  TYPE(ising_grid) :: grid
  LOGICAL :: success
  CHARACTER :: init

  err = 0
  CALL fill_command_args(N, max_iters, init, beta, J, err)
  ! Get command line arguments

  IF(err > 0) THEN
    PRINT*, "Fatal Error, command line arguments incorrect"
    STOP
  END IF

  CALL grid%setup_grid(N, init, beta, J)

  init_mag = grid%magnetization()

  !I chose to show the initial state, and wait for Enter to run
  CALL display_3val_array(grid%state, .TRUE.)
  CALL wait_for_enter_key

  !Run the iterations and then display again
  CALL grid%iterate(max_iters)

  CALL display_3val_array(grid%state, .TRUE.)

  total_flips = grid%flips()
  PRINT*,'Flipped ', total_flips, ' states (', REAL(total_flips)/REAL(max_iters)*100, '%)'

  final_mag = grid%magnetization()
  PRINT*, "Initial magnetization", init_mag
  PRINT*, "Final magnetization  ", final_mag, "(", (final_mag/init_mag - 1.0)*100.0, "%)"

END PROGRAM
