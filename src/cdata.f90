module cdata
  !**********VORTEX FILAMENT******************************************************
  type qvort !our main structure
    real :: x(3) !position
    real :: u(3), u1(3), u2(3) !stored velocities (adam bash)
    real :: ghosti(3), ghostb(3)
    integer :: infront, behind !to form a line/loop
    integer :: closest !nearest particle, used in reconnection
    real :: closestd !distance to nearest particle
  end type
  type(qvort), allocatable :: f(:) !main vector
  integer :: pcount !number of particles in the simulation
  !**********MESH STRUCTURE********************************************************
  type grid
    real :: x(3) !position
    real :: u_sup(3) !velocity (superfluid)
    real :: u_norm(3) !velocity (normal fluid)
  end type
  type(grid), allocatable :: mesh(:,:,:)
  real :: mesh_delta !mesh resolution
  !**********PARTICLE STRUCTURE****************************************************
  type quasi !quasi particle structure
    real :: x(3) !position
    real :: u(3), u1(3), u2(3) !stored velocities (adam bash)
  end type
  type(quasi), allocatable :: g(:) !vector of particles
  !**************TIME PARAMS*******************************************************
  real :: t=0. !hold the current time globally
  integer :: itime !current timestep
  integer :: nstart=1 !integer loop starts from (altered by reading in stored data)

  !***********DIAGNOSTIC INFO******************************************************
  integer :: recon_count=0 !total number of reconnections
  real :: total_length !total length of filaments
  real :: avg_sep !average separation of the particles
  real :: maxu,maxdu !velocity information
  real :: energy !vortex energy

  !***********CONSTANTS************************************************************
  !some constants - precompute for speed
  real, parameter :: pi=3.14159265358979324
  real, parameter :: one_half = (1./2.)
  real, parameter :: three_twos=(3./2.)
  real, parameter :: twenty_three_twelve=(23./12.)
  real, parameter :: four_thirds=(4./3.)
  real, parameter :: five_twelths=(5./12.)

  !***********RUN.IN***************************************************************
  !parameters from run.in, given protected status so treated like parameters
  !by routines in the rest of the code...
  !--------main parameters-please set thes in run.in-----------------------------
  integer, protected :: nsteps, shots, recon_shots=1
  integer, protected :: init_pcount
  real, protected :: dt, delta
  real, protected :: box_size=0.
  real, protected :: quant_circ
  integer, protected :: mesh_size=0
  logical :: periodic_bc=.false.
  character(len=30), protected :: velocity, initf
  integer, protected :: line_count=0

  !--------the following parameters add special features-------------------------
  !---------these should all have default values which 'switch' them off---------
  !------------normal fluid component--------------------------------------------
  character(len=30), protected :: normal_velocity='zero'
  real, protected :: alpha(2)=0. !mutual friction coefficients
  !-----------------forcing------------------------------------------------------
  character(len=20), protected :: force='off'
  real, protected :: force_amp=0.
  real, protected :: force_freq=0.  
  !-----------------special data dumps-------------------------------------------
  !do we want to dump 'f' at a specific time, i.e. before a reconnection etc.
  real, protected :: special_dump=0. !special dump time
  integer :: int_special_dump=0. !special dump time integer
  !---------------------particles------------------------------------------------
  integer, protected :: quasi_pcount=0 !number of particles (quasi or fluid)
  character(len=20), protected :: particle_type='fluid' !fluid/interial/quasi particles
  character(len=20), protected :: initg='random' !initial particle configuration
  !---------------------tree-code------------------------------------------------
  real, protected :: tree_theta=0.
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
          case ('quant_circ')
             read(buffer, *, iostat=ios) quant_circ !quantum of circulation
          case ('box_size')
             read(buffer, *, iostat=ios) box_size !size of periodic box
          case ('mesh_size')
             read(buffer, *, iostat=ios) mesh_size !size of mesh
          case ('velocity')
             read(buffer, *, iostat=ios) velocity !BS/LIA
          case ('normal_velocity')
             read(buffer, *, iostat=ios) normal_velocity !zero/xflow/ABC/KS
          case ('alpha')
             read(buffer, *, iostat=ios) alpha !mutual friction
          case ('initf')
             read(buffer, *, iostat=ios) initf !initial setup of filaments
          case ('initg')
             read(buffer, *, iostat=ios) initg !initial setup of particles
          case ('line_count')
             read(buffer, *, iostat=ios) line_count !used in certain intial conditions
          case ('force')
             read(buffer, *, iostat=ios) force !force the vortices
          case ('force_amp')
             read(buffer, *, iostat=ios) force_amp !forcing amplitude
          case ('force_freq')
             read(buffer, *, iostat=ios) force_freq !forcing frequency
          case ('special_dump')
             read(buffer, *, iostat=ios) special_dump !special dump
          case ('quasi_pcount')
             read(buffer, *, iostat=ios) quasi_pcount !number of particles
          case ('particle_type')
             read(buffer, *, iostat=ios) particle_type !particle type (quasi/fluid)
          case ('tree_theta')
             read(buffer, *, iostat=ios) tree_theta !tree code, opening angle
          case default
             !print *, 'Skipping invalid label at line', line
          end select
       end if
    end do
  end subroutine
  !*************************************************************************************************  
  subroutine fatal_error(location,message)
    implicit none      
    character(len=*) :: location
    character(len=*) :: message
    write (*,*) 'FATAL ERROR'
    write (*,*) trim(location) , ": " , trim(message)
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
