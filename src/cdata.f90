!>hold main variable globally, all other modules will need this module included
!>to access the filament/particle arrays and other important variables
!>also contains important routines for initialising random number generator and
!>reading in runfile
module cdata
  !**********VORTEX FILAMENT******************************************************
  !>our main structure which holds vortex points
  !!@param x position of the vortex poin  
  !!@param u velocity of the vortex point
  !!@param u1 @param u2 stored velocities for Adams-Bashforth
  !!@param u_sup nice to have the just the superfluid veloctity even with normal fluid/forcing
  !!@param ghosti @param ghostb ghost particles for periodic b.c's
  !!@param infront @param behind flag to make points an orientated filament
  !!@param closest closest particle, used in reconnections
  !!@param sph particle attached to, 0 if not
  !!@param sph2 neighbouring sph particles for hermite interpolation
  !!@param closestd separation between closest particle and particle 
  !!@param B magnetic field strength
  !!@param delta used for adaptive meshing along the filaments, used as a prefactor
  !!@param l1 l2 line length for magnetic field evolution 
  !!@param v1 v2 volume element for magnetic field evolution in a compressible field
  type qvort 
    real :: x(3)
    real :: u(3), u1(3), u2(3) 
    real :: u_sup(3), u_mf(3)
    real :: ghosti(3), ghostb(3)
    integer :: infront, behind 
    integer :: closest
    integer :: sph  
    logical :: pinnedi=.false.
    logical :: pinnedb=.false.
    real :: closestd
    real :: delta
    real :: B 
    real :: l1, l2
    real :: v1, v2
  end type
  !>main filament vector
  type(qvort), allocatable :: f(:) 
  !>number of vortex points in the simulation - size of f is pcount
  integer :: pcount 
  !**********MESH STRUCTURE********************************************************
  !>mesh structure - normal and superfluid velocities are calculated on the mesh
  !>which runs from -box_size /2 to box_size/2 to with mesh_size^3
  !>points, note by setting mesh_size>0 mesh is switched on
  !>mesh velocities calculated every mesh_shots timesteps
  !>@param x position of mesh point
  !>@param u_sup superfluid velocity at meshpoint
  !>@param u_norm normal velocity at meshpoint
  type grid
    real :: x(3) 
    real :: u_sup(3) 
    real :: u_norm(3) 
  end type
  !>3D allocatable mesh
  type(grid), allocatable :: mesh(:,:,:)
  !>the mesh resolution box_size/mesh_size^3
  real :: mesh_delta
  !**********QUASI PARTICLE STRUCTURE****************************************************
  !>particle structure - all quasi particle information held within
  !>@param x the position of the particle
  !>@param xold array of 6 previous positions used in backwards difference scheme, for stiff problems
  !>@param p momentum of particle - only used when a quasi particle
  !>@param pold array of 6 previous momenta used in backwards difference scheme
  !>param rdot_old old quasi particle velocity
  !>param energy energy of a quasi particle
  !>param rdot current change of position stored for diagnostic printing
  !>param pdot current change of momentum stored for diagnostic printing
  type quasi 
    real :: x(3) 
    real :: xold(6,3) 
    real :: p(3) 
    real :: pold(6,3) 
    real :: energy
    real :: rdot(3)
    real :: pdot(3)
  end type
  !>vector of quasi particles
  type(quasi), allocatable :: g(:)
  !**********PARTICLE STRUCTURE****************************************************
  !>particle structure - all particle information held within
  !>@param x the position of the particle
  !>@param u the velocity of the particle
  !>@param u1 @param u2 old velocities for Adams-Bashforth timestepping
  type particella
    real :: x(3) 
    real :: u(3), u1(3), u2(3) 
  end type
  !>vector of particles
  type(particella), allocatable :: p(:)  
  !**********SPH STRUCTURE****************************************************
  !>SPH particle structure
  !>@param x the position of the particle
  !>@param i i=s(i)%i used in tree neighbour list
  !>@param m the mass of the particle
  !>@param rho the density at the particle position
  !>@param drhodh \f$ \partial \rho_i / \partial h_i\f$
  !>@param P the pressure at i
  !>@param h the smoothing length associated with the particle
  !>@param f the correction to the smoothing length
  !>@param u the velocity of the particle
  !>@param divu the divergence of the velocity field
  !>@param a the acceleration of the particle
  !>@param u1 @param u2 old velocities for Adams-Bashforth timestepping
  !>@param a1 @param a2 old velocities for Adams-Bashforth timestepping
  type smooth_particle
    real :: x(3)
    integer :: i
    real :: m
    real :: rho, drhodh
    real :: P
    real :: h
    real :: f
    real :: divu
    real :: u(3), u1(3), u2(3) 
    real :: a(3), a1(3), a2(3)
    integer :: ncount !number of neighbours
  end type
  !>vector of SPH particles
  type(smooth_particle), allocatable :: s(:)
  !>nearest neighbour using a linked list
  !>linked list type1
  type neigh_ll
    integer :: i
    type( neigh_ll ), pointer :: next
  end type neigh_ll
  !>linked list type2, which needs above
  type use_neigh_ll
    integer :: ncount
    type( neigh_ll ), pointer :: list, current, previous
  end type use_neigh_ll
  type(use_neigh_ll), allocatable :: s_NN(:)
  !**************TIME PARAMS*******************************************************
  !>time held globally
  real :: t=0. 
  !>current timestep
  integer :: itime 
  !>integer loop starts from (altered by reading in stored data - i.e. restarting)
  integer :: nstart=1 
  !***********DIAGNOSTIC INFO******************************************************
  !>total number of reconnections
  integer :: recon_count=0 
  !>total number of particle removals due to contraction of filament
  integer :: remove_count=0 
  !>total length of filaments
  real :: total_length
  !>average separation of the vortex points 
  real :: avg_sep
  !>maximum velocity
  real :: maxu
  !>maximum velocity change - acceleration x dt 
  real :: maxdu
  !>vortex energy
  real :: energy 
  !>mean curvature
  real :: kappa_bar 
  real :: kappa_min, kappa_max !min/max curvature
  !>rms of magnetic field
  real :: Brms 
  !>self reconnection count
  integer :: self_rcount=0 
  !>vortex vortex reconnection count 
  integer :: vv_rcount=0
  !>linking number - http://en.wikipedia.org/wiki/Linking_number
  real :: linking_number
  !>writhing number - http://mathworld.wolfram.com/Writhe.html
  real :: writhing_number
  !***********CONSTANTS************************************************************
  !some constants - precompute for speed
  real, parameter :: pi=3.14159265358979324
  real, parameter :: rootpi=1.77245385090551
  real, parameter :: one_half = (1./2.)
  real, parameter :: three_twos=(3./2.)
  real, parameter :: twenty_three_twelve=(23./12.)
  real, parameter :: four_thirds=(4./3.)
  real, parameter :: five_twelths=(5./12.)
  !*********ANYTHING ELSE**************************************************
  logical :: nf_compressible=.false. !is the normal velocity field divergence free?
  !***********RUN.IN***************************************************************
  !parameters from run.in, given protected status so treated like parameters
  !by routines in the rest of the code...
  !--------main parameters-please set thes in run.in-----------------------------
  integer, protected :: nsteps, shots, recon_shots=1
  integer, protected :: init_pcount
  real, protected ::  delta
  !>timestep
  real :: dt
  real, protected :: box_size=0.
  real, protected :: quant_circ=9.97E-4
  real , protected :: corea=8.244023E-9
  !>are boundaries periodic?  
  logical :: periodic_bc=.false.
  logical :: periodic_bc_notx=.false.
  real :: xdim_scaling_factor=1.
  !>are boundaries solid?
  logical :: mirror_bc=.false.
  !key arguements that must be set
  character(len=30), protected :: velocity, initf, boundary
  !-----------------arguements used by initial.mod-------------------------
  integer, protected :: line_count=1
  real, protected :: lattice_ratio=1
  real, protected :: rotation_factor=1 !also used in injection routines
  integer, protected :: wave_count=1
  real, protected :: wave_slope=-1.5
  real, protected :: wave_amp=10.
  character(len=30), protected :: wave_type='planar' 
  !--------the following parameters add special features-------------------------
  !---------these should all have default values which 'switch' them off---------
  logical, protected :: binary_print=.true.
  logical, protected :: dt_adapt=.false.
  logical, protected :: delta_adapt=.false.
  logical, protected :: delta_adapt_print=.false.
  !--------------------simulate phonon emission at high k------------------------
  logical, protected :: phonon_emission=.false. !do we want it one?
  real, protected :: phonon_percent=0.95 !what percentage of 2/delta?
  !----------------------mesh information----------------------------------------
  integer, protected :: mesh_size=0
  integer, protected :: mesh_shots=100
  !------------normal fluid component--------------------------------------------
  character(len=30), protected :: normal_velocity='zero'
  real, protected :: alpha(2)=0. !mutual friction coefficients
  real, protected :: normal_fluid_cutoff=1E8 !impossibly high time
  !------------KS model--------------------------------------------
  integer,protected :: KS_rey_int=8
  real,protected :: KS_slope=-5./3.
  integer, protected :: KS_modes=50
  !-----------------forcing------------------------------------------------------
  character(len=20), protected :: force='off'
  real, protected :: force_amp=0.
  real, protected :: force_freq=0.  
  !-----------------special data dumps-------------------------------------------
  !do we want to dump 'f' at a specific timeu1, i.e. before a reconnection etc.
  real, protected :: special_dump=0. !special dump time
  integer :: int_special_dump=0. !special dump time integer
  !---------------------quasi particles------------------------------------------------
  integer :: quasi_pcount=0 !number of quasi particles
  character(len=20), protected :: initg='random' !initial quasi particle configuration
  !---------------------particles------------------------------------------------
  integer :: part_count !number of particles 
  real :: part_stokes=0. !the stokes number of the particles
  character(len=20), protected :: particle_type='fluid' !fluid/interial/quasi particles
  character(len=20), protected :: initp='random' !initial particle configuration
  logical, protected :: particles_only=.false. !only evolve particles in the code
  !----------------------------SPH-----------------------------------------------
  integer :: SPH_count !number of SPH particles in the code, this may reduce due to mergers
  real,protected :: SPH_mass=0. !initial mass of the particle set to 0 which stops code runnning
  character(len=20),protected :: SPH_init='random' !initial setup of SPH particles
  real,protected :: SPH_init_r=0.1 !used for some initial SPH conditions
  real,protected :: SPH_gamma=5./3. !adiabatic index
  real,protected :: SPH_G=0. !gravitational constant
  real,protected :: SPH_theta=0. !opening angle for SPH tree
  integer, protected :: SPH_mesh_size=0 !mesh size for SPH
  !---------------------tree-code------------------------------------------------
  real, protected :: tree_theta=0.
  logical, protected :: tree_print=.false.
  !--------------------additional diagnostics------------------------------------
  logical, protected :: curv_hist=.false. !dumps binned curvature information
  logical, protected :: topo_inf=.false. !calculate topological information
  logical, protected :: energy_inf=.false. !calculate energy of vortex 
  integer, protected :: one_dim=0 !size of 1d velocity information printed to file
  integer, protected :: two_dim=0 !size of 2d velocity information printed to file
  logical, protected :: vapor_print=.false. !dumps raw mesh data for vapor 
  logical, protected :: mirror_print=.false. !prints the mirror filaments to file
  logical, protected :: vel_print=.false. !prints the full velocity information to file
  logical, protected :: vel_print_extra=.false. !prints extra velocity information to file
  logical, protected :: full_B_print=.false. !prints the full magnetic field to file
  logical, protected :: recon_info=.false. !more in depth reconnection information
  logical, protected :: boxed_vorticity=.false. !smoothed vorticity in a box
  integer, protected :: boxed_vorticity_size=32 !how big is the mesh for boxed vorticity
  logical, protected :: simple_plots=.false. !call scripts from command line to plot on the fly
  !-------------------------------smoothing-------------------------------------------
  !gaussian smoothing of vorticity/B field
  real, protected :: smoothing_length=1. !length we smooth over
  integer, protected :: sm_size=0 !size of smoothing mesh - 0 by default which deactivates smoothing
  !----------------------------magnetic field-------------------------------------
  !----------------ENABLE THE FILAMENTS TO ACT AS MAGNETIC FLUX TUBES-------------
  logical, protected :: magnetic=.false. !no by defult
  real, protected :: B_init=1. !initial field strength 
  real, protected :: B_nu=0. !switched off by default
  real, protected :: B_tension=0. !magnetic tension coeff.
  logical, protected :: B_3D_nu=.false. !1/3D diffusion 
  !------------------------------filament injection-------------------------------
  integer, protected :: inject_skip=10000000!how often we insert the vortice
  integer, protected :: inject_size=0 !number of points used
  real, protected :: inject_stop=1E8 !when to stop injection - arbitrarily high
  character(len=20),protected :: inject_type='off' !how we inject loops
  !----------------------------code testing---------------------------------------
  logical, protected :: switch_off_recon=.false.!turns of reconnection algorithm
  logical, protected :: seg_fault=.false.!use print statements to try and isolate segmentation faults
  logical, protected :: NAN_test=.true.!test for NANs in arrays
  contains
  !*************************************************************************************************  
  !>read the file run.in obtaining all parameters at runtime, avoiding the need to recompile the code
  subroutine read_run_file()
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
          case ('dt_adapt')
             !adpative timestep, primarily for magnetic runs with tension
             read(buffer, *, iostat=ios) dt_adapt !adpative dt
          case ('delta')
             !resolution, real number
             read(buffer, *, iostat=ios) delta !spatial resolution
          case ('quant_circ')
             !quatum of circulation, not necessary to set, accepts real #
             read(buffer, *, iostat=ios) quant_circ !quantum of circulation
          case ('corea')
             read(buffer, *, iostat=ios) corea !size of vortex core
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
          case ('xdim_scaling_factor')
             read(buffer, *, iostat=ios) xdim_scaling_factor !scale x axis
          case ('normal_velocity')
             read(buffer, *, iostat=ios) normal_velocity !zero/xflow/ABC/KS
          case ('normal_fluid_cutoff')
             read(buffer, *, iostat=ios) normal_fluid_cutoff !turn off nf
          case ('alpha')
             read(buffer, *, iostat=ios) alpha !mutual friction
          case ('initf')
             read(buffer, *, iostat=ios) initf !initial setup of filaments
          case ('initg')
             read(buffer, *, iostat=ios) initg !initial setup of quasi particles
          case ('initp')
             read(buffer, *, iostat=ios) initp !initial setup of particles
          case ('line_count')
             read(buffer, *, iostat=ios) line_count !used in certain intial conditions
          case ('lattice_ratio')
             read(buffer, *, iostat=ios) lattice_ratio !used in lattice initial conditions
          case ('rotation_factor')
             read(buffer, *, iostat=ios) rotation_factor !used in random loops initial cond.
          case ('force')
             read(buffer, *, iostat=ios) force !force the vortices
          case ('force_amp')
             read(buffer, *, iostat=ios) force_amp !forcing amplitude
          case ('force_freq')
             read(buffer, *, iostat=ios) force_freq !forcing frequency
          case ('phonon_emission')
             read(buffer, *, iostat=ios) phonon_emission !phonon emission on or off
          case ('phonon_percent')
             read(buffer, *, iostat=ios) phonon_percent !percentage of max curv we cutoff at
          case ('special_dump')
             read(buffer, *, iostat=ios) special_dump !special dump
          case ('quasi_pcount')
             read(buffer, *, iostat=ios) quasi_pcount !number of quasi particles
          case ('part_count')
             read(buffer, *, iostat=ios) part_count !number of particles
          case ('particle_type')
             read(buffer, *, iostat=ios) particle_type !particle type (fluid/intertial)
          case ('particles_only')
             read(buffer, *, iostat=ios) particles_only !only evolve particles
          case ('part_stokes')
             read(buffer, *, iostat=ios) part_stokes !stokes number of inertial particles
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
          case ('boxed_vorticity')
             read(buffer, *, iostat=ios) boxed_vorticity !do we want vorticity on a mesh?
          case ('boxed_vorticity_size')
             read(buffer, *, iostat=ios) boxed_vorticity_size !what is the size of this mesh?
          case ('mirror_print')
             read(buffer, *, iostat=ios) mirror_print !print the mirror filaments
          case ('vel_print')
             read(buffer, *, iostat=ios) vel_print !print the velocity information
          case ('vel_print_extra')
             read(buffer, *, iostat=ios) vel_print_extra !print extra velocity information
          case ('vapor_print')
             read(buffer, *, iostat=ios) vapor_print !print the velocity field for vapor
          case ('KS_slope')
             read(buffer, *, iostat=ios) KS_slope !KS velocity field spectrum
          case ('KS_rey_int')
             read(buffer, *, iostat=ios) KS_rey_int !KS Reynolds number proxy
          case ('KS_modes')
             read(buffer, *, iostat=ios) KS_modes !the number of KS modes
          case ('one_dim')
             read(buffer, *, iostat=ios) one_dim !size of 1D print
          case ('two_dim')
             read(buffer, *, iostat=ios) two_dim !size of 2D print
          case ('recon_info')
             read(buffer, *, iostat=ios) recon_info !extra reconnection information
          case ('topo_inf')
             read(buffer, *, iostat=ios) topo_inf !topological information
          case ('energy_inf')
             read(buffer, *, iostat=ios) energy_inf !vortex energy - only open boundaries
          case ('switch_off_recon')
             read(buffer, *, iostat=ios) switch_off_recon !for test cases only!
          case ('seg_fault')
             read(buffer, *, iostat=ios) seg_fault !use print statements to find segmentation faults
          case ('simple_plots')
             read(buffer, *, iostat=ios) simple_plots !perform plots on the fly
          case ('NAN_test')
             read(buffer, *, iostat=ios) NAN_test !test for NANs
          case ('smoothing_length')
             read(buffer, *, iostat=ios) smoothing_length !length we are smoothing over (delta)
          case ('sm_size')
             read(buffer, *, iostat=ios) sm_size !size of smoothing mesh
          case ('magnetic')
             read(buffer, *, iostat=ios) magnetic !act as a magnetic field
          case ('B_init')
             read(buffer, *, iostat=ios) B_init !initial field strength
          case ('B_nu')
             read(buffer, *, iostat=ios) B_nu !magnetic diffusivity
          case ('B_3D_nu')
             read(buffer, *, iostat=ios) B_3D_nu !is diffusion 3D?
          case ('B_tension')
             read(buffer, *, iostat=ios) B_tension !magnetic tension coeff.
          case ('full_B_print')
             read(buffer, *, iostat=ios) full_B_print !print full B info
           case ('SPH_count')
             read(buffer, *, iostat=ios) SPH_count !how many SPH particles in the code
          case ('SPH_mass')
             read(buffer, *, iostat=ios) SPH_mass!initial mass of SPH particles
          case ('SPH_init')
             read(buffer, *, iostat=ios) SPH_init!initial SPH setup
          case ('SPH_gamma')
             read(buffer, *, iostat=ios) SPH_gamma!adiabatic index in SPH sims.
          case ('SPH_G')
             read(buffer, *, iostat=ios) SPH_G !gravitational constant
          case ('SPH_init_r')
             read(buffer, *, iostat=ios) SPH_init_r !initial sphere radius
          case ('SPH_theta')
             read(buffer, *, iostat=ios) SPH_theta !tree opening angle
          case ('SPH_mesh_size')
             read(buffer, *, iostat=ios) SPH_mesh_size !SPH mesh size
          case ('delta_adapt')
             read(buffer, *, iostat=ios) delta_adapt !is the discretisation adaptive
          case ('delta_adapt_print')
             read(buffer, *, iostat=ios) delta_adapt_print !print apative discretisation
          case ('inject_size')
             read(buffer, *, iostat=ios) inject_size !size of injected filaments
          case ('inject_skip')
             read(buffer, *, iostat=ios) inject_skip !how often we inject
          case ('inject_type')
             read(buffer, *, iostat=ios) inject_type !how we inject new vortices
          case ('inject_stop')
             read(buffer, *, iostat=ios) inject_stop !when (if ever) we stop injecting 
          case default
             !print *, 'Skipping invalid label at line', line
          end select
       end if
    end do
  end subroutine
  !*************************************************************************************************  
  !>reload the file run.in and show differences
  subroutine reload_run_file()
    implicit none
    ! input related variables
    character(len=100) :: buffer, label
    integer :: pos
    integer, parameter :: fh = 15
    integer :: ios = 0
    integer :: line = 0
    !****************************************
    !*******RELOAD DUMMY VARIABLES***********
    integer :: reload_nsteps
    integer :: reload_shots,reload_recon_shots, reload_mesh_shots
    real ::  reload_normal_fluid_cutoff, reload_inject_stop
    logical :: reload_curv_hist, reload_vel_print, reload_vapor_print
    logical :: reload_topo_inf, reload_energy_inf, reload_recon_info
    integer :: reload_one_dim, reload_two_dim
    open(fh, file='run.in')
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
             read(buffer, *, iostat=ios) reload_nsteps
             if (reload_nsteps/=nsteps) then
               nsteps=reload_nsteps
               write(*,'(a,i6.6)') 'RELOAD: changed nsteps to ', nsteps
             end if
          case ('shots')
             read(buffer, *, iostat=ios) reload_shots
             if (reload_shots/=shots) then
               shots=reload_shots
               write(*,'(a,i6.6)') 'RELOAD: changed shots to ', shots
             end if 
          case ('recon_shots')
             read(buffer, *, iostat=ios) reload_recon_shots 
             if (reload_recon_shots/=recon_shots) then
               recon_shots=reload_recon_shots
               write(*,'(a,i6.6)') 'RELOAD: changed recon_shots to ', recon_shots
             end if 
          case ('mesh_shots')
             read(buffer, *, iostat=ios) reload_mesh_shots 
             if (reload_mesh_shots/=mesh_shots) then
               mesh_shots=reload_mesh_shots
               write(*,'(a,i6.6)') 'RELOAD: changed mesh_shots to ', mesh_shots
             end if 
          case ('normal_fluid_cutoff')
             read(buffer, *, iostat=ios) reload_normal_fluid_cutoff 
             if ((reload_normal_fluid_cutoff>normal_fluid_cutoff).or. & 
                (reload_normal_fluid_cutoff<normal_fluid_cutoff)) then
               normal_fluid_cutoff=reload_normal_fluid_cutoff
               write(*,'(a,f10.4)') 'RELOAD: changed normal_fluid_cutoff to ', normal_fluid_cutoff
             end if 
          case ('curv_hist')
             read(buffer, *, iostat=ios) reload_curv_hist 
             if (reload_curv_hist.neqv.curv_hist) then
               curv_hist=reload_curv_hist
               write(*,*) 'RELOAD: changed curv_hist to ', curv_hist
             end if 
          case ('vel_print')
             read(buffer, *, iostat=ios) reload_vel_print 
             if (reload_vel_print.neqv.vel_print) then
               vel_print=reload_vel_print
               write(*,*) 'RELOAD: changed vel_print to ', vel_print
             end if 
          case ('vapor_print')
             read(buffer, *, iostat=ios) reload_vapor_print 
             if (reload_vapor_print.neqv.vapor_print) then
               vapor_print=reload_vapor_print
               write(*,*) 'RELOAD: changed vapor_print to ', vapor_print
             end if 
          case ('one_dim')
             read(buffer, *, iostat=ios) reload_one_dim
             if (reload_one_dim/=one_dim) then
               one_dim=reload_one_dim
               write(*,*) 'RELOAD: changed one_dim to ', one_dim
             end if 
          case ('two_dim')
             read(buffer, *, iostat=ios) reload_two_dim 
             if (reload_two_dim/=two_dim) then
               two_dim=reload_two_dim
               write(*,*) 'RELOAD: changed two_dim to ', two_dim
             end if 
          case ('recon_info')
             read(buffer, *, iostat=ios) reload_recon_info
             if (reload_recon_info.neqv.recon_info) then
               recon_info=reload_recon_info
               write(*,*) 'RELOAD: changed recon_info to ', recon_info
             end if 
          case ('topo_inf')
             read(buffer, *, iostat=ios) reload_topo_inf 
             if (reload_topo_inf.neqv.topo_inf) then
               topo_inf=reload_topo_inf
               write(*,*) 'RELOAD: changed topo_inf to ', topo_inf
             end if 
          case ('energy_inf')
             read(buffer, *, iostat=ios) reload_energy_inf 
             if (reload_energy_inf.neqv.energy_inf) then
               energy_inf=reload_energy_inf
               write(*,*) 'RELOAD: changed energy_inf to ', energy_inf
             end if 
          case ('inject_stop')
             read(buffer, *, iostat=ios) reload_inject_stop 
             if ((reload_inject_stop>inject_stop).or. & 
                (reload_inject_stop<inject_stop)) then
               inject_stop=reload_inject_stop
               write(*,'(a,f10.4)') 'RELOAD: changed inject_stop to ', inject_stop
             end if 
          case default
          end select
       end if
    end do
  end subroutine
  !*************************************************************************************************  
  !>print error message to screen and stop the run
  !!provide a location of error and the message to print
  subroutine fatal_error(location,message)
    implicit none      
    character(len=*) :: location
    character(len=*) :: message
    write (*,*) '-------------------------FATAL ERROR-------------------------'
    write (*,*) trim(location) , ": " , trim(message)
    write (*,*) "FYI: t= ", t, "iteration number= ", itime 
    write (*,*) '-------------------------------------------------------------'
    stop
  end subroutine
  !*************************************************************************************************  
  !>print a warning message to screen - will not stop the code
  !!provide a location of error and the message to print
  subroutine warning_message(location,message)
    implicit none      
    character(len=*) :: location
    character(len=*) :: message
    write (*,*) '-------------------------WARNING----------------------------'
    write (*,*) trim(location) , ": " , trim(message)
    write (*,*) '------------------------------------------------------------'
  end subroutine
  !*************************************************************************************************  
  !>generate a new random seed, or read one in from ./data if restating code
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
!>\page BUG Code testing and bug fixing
!>The code uses a system of messaging subroutines to kill the code
!>these are warning message, which will print a warning message to screen
!>and create an empty file called WARNING to let you know something has gone 
!>wrong, but will not end the run. fatal_error on the other hand should be used 
!>in places where the code is killed. These routines are located in cdata.mod.
!>
!>A routine called NAN_finder (located in general.mod) is called at the end of
!>each timestep and trys to find NANs (not a number) in any of the main arrays.
!>If one is called then fatal_error is called telling the user where the NAN is
!>located.
!>
!>Finally a few testing flags can be set in run.in at present the following are
!>available (set T to turn on):
!>
!>- \p switchoff_recon - turn off reconnection algorithm.
!>- \p seg_fault - will use print statements place throughout run.f90 to try and
!>help the user isolate a segmentation fault in the code.
!>\page RELOAD RELOAD and STOP files
!!Whilst the code is runnning you may wish to force the code to stop runnning. One
!!could kill the exectuable using ctrl+c, however if the code is running in the background
!!(i.e. you have set the code running with ./run.sh -q) then this is not possible.
!!To stop the code one can simply create a file STOP (the easiest way to do this is 
!!with the command touch STOP), the code will then stop. \n
!!
!!A more important issue is changing runtime parameters whilst the code is running. You
!!could simply kill the code, make your changes and then restart the code (./run.sh -r), however
!!this is not neat solution.\n
!!Instead make your changes to the run.in file and then create a reload file (touch RELOAD). This
!!will re-read the run.in file allowing you to make changes on the fly.
!!Note not all parameters can be changed below is a list of parameters that can be re-read
!!- \p nsteps - number of timesteps the code runs for
!!- \p shots - how often the code prints to file (filament only not meshes)
!!- \p mesh_shots - how often we print mesh's (3,2 and 1D to file)
!!- \p normal_fluid_cutoff - when (if) we turn off the normal fluid drive
!!- \p curv_hist - print binned curvature information?
!!- \p vel_print - print velocity information to run angela's velocity stat scripts
!!- \p vapor_print - print meshes as binary files speciffically for vapor?
!!- \p one_dim - 1D mesh size (runs in x-direction, y=0,z=0)
!!- \p two_dim - 2D mesh size (in xy-plane, z=0)
!!- \p recon_info - print specific reconnection information
!!- \p topo_inf - calculate topological information of the filament (linking and writhe) very slow!
!!- \p energy_inf - print energy information to file 
!!- \p inject_stop - if (when) the code stops injecting loops (if inject has been set)
