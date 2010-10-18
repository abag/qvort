module cdata
  !THIS MODULE HOLDS ALL THE VARIABLES GLOBALLY, ALL MODULES WILL NEED THIS
  !THIS ROUTINE CONTAINS SOME IMPORTANT ROUTINES TO READ IN THE RUNTIME PARAMETERS
  !IN RUN.IN AS WELL AS SETTING THE RANDOM NUMBER GENERATOR
  !**********VORTEX FILAMENT******************************************************
  type qvort !our main structure
    real :: x(3) !position
    real :: u(3), u1(3), u2(3) !stored velocities (adam bash)
    real :: ghosti(3), ghostb(3)
    integer :: infront, behind !to form a line/loop
    integer :: closest !nearest particle, used in reconnection
    real :: closestd !distance to nearest particle
    real :: delta
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
  integer :: remove_count=0 !total number of particle removals
  real :: total_length !total length of filaments
  real :: avg_sep !average separation of the particles
  real :: maxu,maxdu !velocity information
  real :: energy !vortex energy
  real :: kappa_bar !mean curvature
  real :: kappa_min, kappa_max !min/max curvature
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
  real, protected ::  delta
  real :: dt
  real, protected :: box_size=0.
  real, protected :: quant_circ=9.97E-4
  logical :: periodic_bc=.false.
  logical :: mirror_bc=.false.
  character(len=30), protected :: velocity, initf, boundary
  logical, protected :: binary_print=.true.
  !--specific initial conditions selected in run.in - give these default values--
  integer, protected :: line_count=1
  integer, protected :: wave_count=1
  real, protected :: wave_slope=-1.5
  real, protected :: wave_amp=10.
  character(len=30), protected :: wave_type='planar' 
  !--------the following parameters add special features-------------------------
  !---------these should all have default values which 'switch' them off---------
  !----------------------mesh information----------------------------------------
  integer, protected :: mesh_size=0
  integer, protected :: mesh_shots=100
  !------------normal fluid component--------------------------------------------
  character(len=30), protected :: normal_velocity='zero'
  real, protected :: alpha(2)=0. !mutual friction coefficients
  real, protected :: normal_fluid_cutoff=1E8 !impossibly high time 
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
  logical, protected :: tree_print=.false.
  !--------------------additional diagnostics------------------------------------
  logical, protected :: curv_hist=.false. !dumps binned curvature information
  logical, protected :: mirror_print=.false. !prints the mirror filaments to file
  logical, protected :: vel_print=.false. !prints the full velocity information to file
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
             !number of steps, enter a natural number
             read(buffer, *, iostat=ios) nsteps !number of steps to take
          case ('shots')
             !number often to print to file, enter a natural number
             read(buffer, *, iostat=ios) shots !how often to print to file
          case ('recon_shots')
             !how often to try reconnection algorithm, enter a natural number
             read(buffer, *, iostat=ios) recon_shots !how often to perform reconnection algorithm
          case ('pcount')
             !initial number of particles, can be overwritten by intial condition, enter natural number
             read(buffer, *, iostat=ios) init_pcount !initial particle count
          case ('dt')
             !size of timestep, cheked based on delta
             read(buffer, *, iostat=ios) dt !timestep value
          case ('binary_print')
             !print to binary (T) or formatted data (F)
             read(buffer, *, iostat=ios) binary_print !print binary var data
          case ('delta')
             !resolution, real number
             read(buffer, *, iostat=ios) delta !spatial resolution
          case ('quant_circ')
             !quatum of circulation, not necessary to set, accepts real #
             read(buffer, *, iostat=ios) quant_circ !quantum of circulation
          case ('box_size')
             !size of box must be>0, real number
             read(buffer, *, iostat=ios) box_size !size of periodic box
          case ('mesh_size')
             !mesh for outputting veloctiy fields, by default is 0, enter
             !natural number
             read(buffer, *, iostat=ios) mesh_size !size of mesh
          case ('mesh_shots')
             read(buffer, *, iostat=ios) mesh_shots !how often to print mesh to file
          case ('velocity')
             !velocity field options are LIA, BS, Tree
             read(buffer, *, iostat=ios) velocity !BS/LIA/Tree
          case ('boundary')
             read(buffer, *, iostat=ios) boundary !open/periodic/mirror
          case ('normal_velocity')
             read(buffer, *, iostat=ios) normal_velocity !zero/xflow/ABC/KS
          case ('normal_fluid_cutoff')
             read(buffer, *, iostat=ios) normal_fluid_cutoff !turn off nf
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
          case ('tree_print')
             read(buffer, *, iostat=ios) tree_print !print the tree mesh
          case ('wave_count')
             read(buffer, *, iostat=ios) wave_count !for wave_spec initial conditions
          case ('wave_slope')
             read(buffer, *, iostat=ios) wave_slope !for wave_spec initial conditions
          case ('wave_amp')
             read(buffer, *, iostat=ios) wave_amp !for wave_spec initial conditions
          case ('wave_type')
             read(buffer, *, iostat=ios) wave_type !for wave_spec initial conditions
          case ('curv_hist')
             read(buffer, *, iostat=ios) curv_hist !do we want binned curvature info?
          case ('mirror_print')
             read(buffer, *, iostat=ios) mirror_print !print the mirror filaments
          case ('vel_print')
             read(buffer, *, iostat=ios) vel_print !print the velocity information
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
    write (*,*) '-------------------------FATAL ERROR-------------------------'
    write (*,*) trim(location) , ": " , trim(message)
    write (*,*) '-------------------------------------------------------------'
    stop
  end subroutine
  !*************************************************************************************************  
  subroutine warning_message(location,message)
    implicit none      
    character(len=*) :: location
    character(len=*) :: message
    write (*,*) '-------------------------WARNING----------------------------'
    write (*,*) trim(location) , ": " , trim(message)
    write (*,*) '------------------------------------------------------------'
  end subroutine
  !*************************************************************************************************  
  subroutine init_random_seed()
     !CREATE A NEW RANDOM SEED, UNLESS RESTARTING CODE
     integer :: i, n=8, clock
     integer, allocatable :: seed(:)
     logical :: seed_exists=.false.
     allocate(seed(n))
     
     call system_clock(count=clock)
     
     seed = clock + 37 * (/ (i - 1, i = 1, n) /)
     !read the seed in from file if possible
     inquire(file='./data/seed.dat',exist=seed_exists)
     if (seed_exists) then
       write(*,*)'reading in random seed from ./data'
       open(unit=37,file='./data/seed.dat',form='unformatted')
         read(37) seed
       close(37)
     else
       write(*,*)'generating new random seed saving to ./data'
       open(unit=37,file='./data/seed.dat',form='unformatted')
         write(37) seed
       close(37)
     end if
     call random_seed(put = seed)
     deallocate(seed)
  end subroutine 
  !**************************************************************************************************  
end module
