module cdata
  type qvort !our main structure
    real :: x(3) !position
    real :: u(3), u1(3), u2(3) !stored velocities (adam bash)
    real :: ghosti(3), ghostb(3)
    integer :: infront, behind !to form a line/loop
  end type
  type(qvort), allocatable :: f(:) !main vector
  integer :: pcount !number of particles in the simulation
  real :: t=0. !hold the current time globally
  integer :: itime !current timestep
  integer :: recon_count=0 !total number of reconnections
  real :: total_length !total length of filaments
  real :: avg_sep !average separation of the particles
  !some constants - precompute for speed
  real, parameter :: pi=3.14159265358979324
  real :: one_half = (1./2.)
  real :: three_twos=(3./2.)
  real :: twenty_three_twelve=(23./12.)
  real :: four_thirds=(4./3.)
  real :: five_twelths=(5./12.)
  integer :: nstart=1 !integer loop starts from (altered by reading in stored data)
  !parameters from run.in, given protected status so treated like parameters
  !by routines in the rest of the code...
  integer, protected :: nsteps, shots, recon_shots=1
  integer, protected :: init_pcount
  real, protected :: dt, delta
  real, protected :: box_size=0.
  logical :: periodic_bc=.false.
  character(40), protected :: velocity, initf
  integer, protected :: line_count=0
  !do we want to dump 'f' at a specific time, i.e. before a reconnection etc.
  real, protected :: special_dump=0. !special dump time
  integer :: int_special_dump=0. !special dump time integer
  contains
  !*************************************************************************************************  
  subroutine read_run_file()
  ! this subroutine reads the file run.in obtaining key runtime parameters
    implicit none
    ! input related variables
    character(len=100) :: buffer, label
    integer :: pos
    integer, parameter :: fh = 15
    integer :: ios = 0
    integer :: line = 0

    open(fh, file='run.in')

    ! ios is negative if an end of record condition is encountered or if
    ! an endfile condition was detected.

    do while (ios == 0)
       read(fh, '(A)', iostat=ios) buffer
       if (ios == 0) then
          line = line + 1

          ! Find the first instance of whitespace.  Split label and data.
          pos = scan(buffer, '     ')
          label = buffer(1:pos)
          buffer = buffer(pos+1:)
          select case (label)
          case ('nsteps')
             read(buffer, *, iostat=ios) nsteps !number of steps to take
          case ('shots')
             read(buffer, *, iostat=ios) shots !how often to print to file
          case ('recon_shots')
             read(buffer, *, iostat=ios) recon_shots !how often to perform reconnection algorithm
          case ('init_pcount')
             read(buffer, *, iostat=ios) init_pcount !initial particle count
          case ('dt')
             read(buffer, *, iostat=ios) dt !timestep value
          case ('delta')
             read(buffer, *, iostat=ios) delta !spatial resolution
          case ('box_size')
             read(buffer, *, iostat=ios) box_size !size of periodic box
          case ('velocity')
             read(buffer, *, iostat=ios) velocity !BS/LIA
          case ('initf')
             read(buffer, *, iostat=ios) initf !initial setup of filaments
          case ('line_count')
             read(buffer, *, iostat=ios) line_count !used in certain intial conditions
          case ('special_dump')
             read(buffer, *, iostat=ios) special_dump !special dump
          case default
             !print *, 'Skipping invalid label at line', line
          end select
       end if
    end do
  end subroutine
  !*************************************************************************************************  
  subroutine fatal_error()
    implicit none
    print*, 'ABORTING RUN'
    stop
  end subroutine
  !*************************************************************************************************  
  SUBROUTINE init_random_seed()
     INTEGER :: i, n, clock
     INTEGER, DIMENSION(:), ALLOCATABLE :: seed
     
     CALL RANDOM_SEED(size = n)
     ALLOCATE(seed(n))
     
     CALL SYSTEM_CLOCK(COUNT=clock)
     
     seed = clock + 37 * (/ (i - 1, i = 1, n) /)
     CALL RANDOM_SEED(PUT = seed)
     
     DEALLOCATE(seed)
  END SUBROUTINE 
  !**************************************************************************************************  
end module
