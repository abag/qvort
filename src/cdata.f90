!>hold main variable globally, all other modules will need this module included
!>to access the filament/particle arrays and other important variables
!>also contains important routines for initialising random number generator and
!>reading in runfile
module cdata
  use omp_lib
  !**********VORTEX FILAMENT******************************************************
  !>our main structure which holds vortex points
  !!@param x position of the vortex point  
  !!@param u velocity of the vortex point
  !!@param u1 @param u2 stored velocities for Adams-Bashforth
  !!@param u_sup nice to have the just the superfluid veloctity even with normal fluid/forcing
  !!@param ghosti @param ghostb ghost particles for periodic b.c's
  !!@param ghostii @param ghostbb ghost particles for periodic b.c's  
  !!@param infront @param behind flag to make points an orientated filament
  !!@param closest closest particle, used in reconnections
  !!@param closestd separation between closest particle and particle 
  !!@param closestd_loop as above but not on same loop 
  !!@param delta used for adaptive meshing along the filaments, used as a prefactor
  type qvort 
    real :: x(3)
    real :: u(3), u1(3), u2(3), u3(3)
    real :: u_sup(3), u_mf(3)
    real :: ghosti(3), ghostb(3)
    real :: ghostii(3), ghostbb(3)
    integer :: infront, behind 
    integer :: closest
    real :: closestd
    real :: closestd_loop
    logical :: pinnedi=.false.
    logical :: pinnedb=.false.
    real :: delta
    real :: t_recon(2)
  end type
  !>main filament vector
  type(qvort), allocatable :: f(:) 
  !>number of vortex points in the simulation - size of f is pcount
  integer :: pcount
  !----------------------------------------------------------------------
  !centre of vorticity data
  type centre_of_vorticity
    real :: x(3), x_old(3), u(3)
  end type
  type(centre_of_vorticity) :: cov
  !macro ring radii data - 2 types
  type macro_ring_radii_1
    real :: x_max(4), x_min(4), x_spread(4)
    real :: r, r_old, r_u, a, a_old, a_u, ra
  end type
  type(macro_ring_radii_1) :: mrr1 
  type macro_ring_radii_2
    real :: r(5), a(5)
  end type
  type(macro_ring_radii_2), allocatable :: mrr2(:)
  !other variables for macro ring radii data - type 2
  real :: avg_r(5), avg_r_old(5), avg_r_u(5)
  real :: avg_a(5), avg_a_old(5), avg_a_u(5)
  real :: avg_ra(5)
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
  type(grid), allocatable :: mesh2D(:,:)
  !>1D allocatable mesh
  type(grid), allocatable :: lat_mesh_1D(:,:)
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
  !**************TIME PARAMS*******************************************************
  !>time held globally
  real :: t=0. 
  !>current timestep
  integer :: itime 
  !>integer loop starts from (altered by reading in stored data - i.e. restarting)
  integer :: nstart=1 
  logical, protected :: restart_rewind=.false.
  !***********DIAGNOSTIC INFO******************************************************
  !>total number of reconnections
  integer :: recon_count=0 
  !>total number of particle removals due to contraction of filament
  integer :: remove_count=0 
  !>total number of removals due to phonon emission
  integer :: phonon_count=0
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
  real, protected :: quant_circ=9.97E-4 !He-4 by default
  real , protected :: corea=8.244023E-9 !He-4 by default
  !>are boundaries periodic?  
  logical :: periodic_bc=.false.
  logical :: periodic_bc_notx=.false.
  logical :: periodic_bc_notxy=.false.
  real :: xdim_scaling_factor=1.
  logical :: sticky_z_boundary=.false.
  !>which reconnection algorithm to use
  character(len=30), protected :: recon_type='original'  
  !>are boundaries solid?
  logical :: mirror_bc=.false.
  !key arguements that must be set
  character(len=30), protected :: velocity, initf, boundary
  !order of derivatives
  character(len=30), protected :: deriv_order='second'
  logical ,protected :: fixed_LIA_beta=.false.
  !-----------arguements used by initial.mod/initial_cond.mod-------------
  integer, protected :: line_count=1
  real, protected :: line_sigma=0.
  real, protected :: lattice_ratio=1
  !how much we translate random_loops initial condition by 
  real, protected :: loop_translate(3)=1.!for separate xyz components
  character(len=20), protected :: initial_distribution='uniform'
  real, protected :: rotation_factor=1 !also used in injection routines
  !--------for wave_loop/wave_line initf--------------
  integer, protected :: wave_count=1 !number of waves
  real, protected :: wave_slope=-1.5 !spectral slope
  integer, protected :: wave_start=1 !starting wavenumber
  integer, protected :: wave_skip=1 !the skip used
  real, protected :: wave_amp=10. !amplitude of 1st wave
  character(len=30), protected :: wave_type='planar' !planar or helical
  !--------for macro_ring initf--------------
  real, protected :: macro_ring_R=0. !major radius
  real, protected :: macro_ring_a=0. !minor radius
  real, protected :: nf_mra_factor=1. !factor to increase size of 'a' in NF
  !--------for hyperboloid initf--------------
  real, protected :: hyperboloid_r=2. !radius of bundle (terms of delta)
  real, protected :: hyperboloid_e=1. !effects curvature of bundle
  !---------for central_bundle------------------------
  character(len=30), protected :: bundle_type='polarised' !polarised or random
  !---------for criss-cross------------------------
  integer, protected :: criss_cross_bundle=1 !typical size of bundles
  real, protected :: criss_cross_width=1. !width in terms of \delta
  !---------for torus_knot initf---------------------- 
  integer, protected :: torus_p=1, torus_q=1 !for torus_knot initf
  real, protected :: torus_epsilon=1.
  !--------the following parameters add special features-------------------------
  !---------these should all have default values which 'switch' them off---------
  logical, protected :: binary_print=.true.
  logical, protected :: dt_adapt=.false.
  logical, protected :: delta_adapt=.false.
  logical, protected :: delta_adapt_print=.false.
  !--------------------simulate phonon emission at high k------------------------
  logical, protected :: phonon_emission=.false. !do we want it one?
  real, protected :: phonon_percent=0.95 !what percentage of 2/delta?
  logical, protected :: hyperviscosity=.false. !instead of smoothing use hyperviscosity
  integer, protected :: hyp_power=4 !degree of hyperviscosity
  real :: hyp_power_dissipate !power dissipated by friction
  real, protected :: hyp_curv=0. !allow some curvature undamped
  real, protected :: hyp_nu=1. !scaling parameter
  !----------------------mesh information----------------------------------------
  integer, protected :: mesh_size=0
  integer, protected :: mesh_shots=100
  logical, protected :: hollow_mesh_core=.false.
  !------------normal fluid component--------------------------------------------
  character(len=30), protected :: normal_velocity='zero'
  real, protected :: alpha(2)=0. !mutual friction coefficients
  real, protected :: normal_fluid_cutoff=1E8 !impossibly high time
  integer, protected :: normal_fluid_freq=1 !used for certain normal fluid flows
  real, protected :: norm_vel_xflow=0.5
  !------------KS model--------------------------------------------
  integer,protected :: KS_rey_int=8
  real,protected :: KS_slope=-5./3.
  integer, protected :: KS_modes=50
  logical, protected :: KS_maximise_rey=.false.
  real, protected :: KS_bubble=0.
  !-----------------forcing------------------------------------------------------
  character(len=30), protected :: force='off'
  real, protected :: force_amp=0.
  real, protected :: force_freq=0.
  real, protected :: force_cutoff=1E8 !impossibly high time 
  logical, protected :: forcing_mesh_fileprint=.false. 
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
  logical, protected ::particle_super_velocity=.false. !do particles feel superlfluid vel
  !---------------------tree-code------------------------------------------------
  real, protected :: tree_theta=0.
  logical, protected :: tree_print=.false.
  logical, protected :: tree_extra_correction=.true.
  !--------------------additional diagnostics------------------------------------
  logical, protected :: curv_hist=.false. !dumps binned curvature information
  logical, protected :: torsion_hist=.false. !dumps binned torsion information
  logical, protected :: topo_inf=.false. !calculate topological information
  logical, protected :: energy_inf=.false. !calculate energy of vortex 
  logical, protected :: sep_inf=.false. !calculate information and histogram of point separation
  logical, protected :: line_sep_inf=.false. !calculate information and histogram of closest line separation
  integer, protected :: one_dim=0 !size of 1d velocity information printed to file
  character(len=1), protected :: one_dim_direction='x' !component we calculate 1D spectra in
  integer, protected :: one_dim_lattice=0 !create a lattice in z direction to calculate vel info on
  integer, protected :: one_dim_lattice_count=1 !number of 1d spectra to draw - must be a square number
  integer, protected :: two_dim=0 !size of 2d velocity information printed to file
  logical, protected :: vapor_print=.false. !dumps raw mesh data for vapor 
  logical, protected :: mirror_print=.false. !prints the mirror filaments to file
  logical, protected :: vel_print=.false. !prints the full velocity information to file
  logical, protected :: vel_print_extra=.false. !prints extra velocity information to file
  logical, protected :: recon_info=.false. !more in depth reconnection information
  logical, protected :: boxed_vorticity=.false. !smoothed vorticity in a box
  integer, protected :: boxed_vorticity_size=32 !how big is the mesh for boxed vorticity
  logical, protected :: simple_plots=.false. !call scripts from command line to plot on the fly
  logical, protected :: anisotropy_params=.false. !get anistropy parameters
  logical, protected :: particle_plane_inf=.false. !2D point infomation in z=0 plane
  logical, protected :: closest_distance=.false. !what is the minimum separation between points?
  logical, protected :: full_loop_counter=.false. !check loop count and size
  logical, protected :: recon_time_info=.false. !print time difference between recon 
  !-------------------------------smoothing-------------------------------------------
  !gaussian smoothing of vorticity/B field
  real, protected :: smoothing_length=1. !length we smooth over
  integer, protected :: sm_size=0 !size of smoothing mesh - 0 by default which deactivates smoothing
  logical, protected :: smoothing_interspace=.false. !smooth using intervortex spacing?
  !------------------------------filament injection-------------------------------
  integer, protected :: inject_skip=10000000!how often we insert the vortice
  integer, protected :: inject_size=0 !number of points used
  integer, protected :: inject_freq=10000000!how often we switch injection direction
  real, protected :: inject_stop=1E8 !when to stop injection - arbitrarily high
  character(len=20),protected :: inject_type='off' !how we inject loops
  !----------------------------code testing---------------------------------------
  logical, protected :: switch_off_recon=.false.!turns of reconnection algorithm
  logical, protected :: seg_fault=.false.!use print statements to try and isolate segmentation faults
  logical, protected :: NAN_test=.true.!test for NANs in arrays
  logical, protected :: overide_timestep_check=.false.!do not perform initial dt test
  !----------------------------warnings---------------------------------------
  integer, private :: max_warn_count=5 !maximum warning count code can survive
  integer, private :: warn_count=0 !number of warnings
  !------------------------------batch mode---------------------------------------
  !this adds the ability to email a user to say when a code has finished or if there is a fatal error
  logical, protected :: batch_mode=.false. !set to true to enable messaging
  character(len=80),protected :: batch_name='qvort run' !what is the name of the run
  character(len=60),protected :: batch_email='a.w.baggaley@gmail.com' !who you gonna call?
  !------------------------------openmp--------------------------------------------
  logical,private :: serial_run=.false.
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
          case ('recon_type')
             read(buffer, *, iostat=ios) recon_type !reconnection algorithm used
          case ('binary_print')
             !print to binary (T) or formatted data (F)
             read(buffer, *, iostat=ios) binary_print !print binary var data
          case ('delta')
             !resolution, real number
             read(buffer, *, iostat=ios) delta !spatial resolution
          case ('quant_circ')
             !quatum of circulation, not necessary to set, accepts real #
             read(buffer, *, iostat=ios) quant_circ !quantum of circulation
          case ('corea')
             read(buffer, *, iostat=ios) corea !size of vortex core
          case ('deriv_order')
             read(buffer, *, iostat=ios) deriv_order !order of spatial derivatives   
          case ('box_size')
             !size of box must be>0, real number
             read(buffer, *, iostat=ios) box_size !size of periodic box
          case ('mesh_size')
             !mesh for outputting veloctiy fields, by default is 0, enter
             !natural number
             read(buffer, *, iostat=ios) mesh_size !size of mesh
          case('hollow_mesh_core')
             read(buffer, *, iostat=ios) hollow_mesh_core
          case ('mesh_shots')
             read(buffer, *, iostat=ios) mesh_shots !how often to print mesh to file
          case ('velocity')
             !velocity field options are LIA, BS, Tree
             read(buffer, *, iostat=ios) velocity !BS/LIA/Tree
          case ('fixed_LIA_beta')
             read(buffer, *, iostat=ios) fixed_LIA_beta !LIA beta is fixed?
          case ('boundary')
             read(buffer, *, iostat=ios) boundary !open/periodic/mirror
          case ('xdim_scaling_factor')
             read(buffer, *, iostat=ios) xdim_scaling_factor !scale x axis
          case ('normal_velocity')
             read(buffer, *, iostat=ios) normal_velocity !zero/xflow/ABC/KS
          case ('normal_fluid_cutoff')
             read(buffer, *, iostat=ios) normal_fluid_cutoff !turn off nf
          case ('force_cutoff')
             read(buffer, *, iostat=ios) force_cutoff !turn off forcing             
          case ('normal_fluid_freq')
             read(buffer, *, iostat=ios) normal_fluid_freq !frequency we "drive" nf
          case ('norm_vel_xflow')
             read(buffer, *, iostat=ios) norm_vel_xflow !counterflow velocity             
          case ('alpha')
             read(buffer, *, iostat=ios) alpha !mutual friction
          case ('initf')
             read(buffer, *, iostat=ios) initf !initial setup of filaments
          case ('initg')
             read(buffer, *, iostat=ios) initg !initial setup of quasi particles
          case ('initp')
             read(buffer, *, iostat=ios) initp !initial setup of particles
          case ('restart_rewind')
             read(buffer, *, iostat=ios) restart_rewind !restart with ntime=1
          case ('line_count')
             read(buffer, *, iostat=ios) line_count !used in certain intial conditions
          case ('initial_distribution')
             read(buffer, *, iostat=ios) initial_distribution !used in certain intial conditions             
          case ('line_sigma')
             read(buffer, *, iostat=ios) line_sigma !used in certain intial conditions
          case ('lattice_ratio')
             read(buffer, *, iostat=ios) lattice_ratio !used in lattice initial conditions
          case ('loop_translate')
             read(buffer, *, iostat=ios) loop_translate !used in random_loops conditions
          case ('bundle_type')
             read(buffer, *, iostat=ios) bundle_type !used in central_bundle initial conditions
          case ('rotation_factor')
             read(buffer, *, iostat=ios) rotation_factor !how much we rotate loops by
          case ('force')
             read(buffer, *, iostat=ios) force !force the vortices
          case ('force_amp')
             read(buffer, *, iostat=ios) force_amp !forcing amplitude
          case ('forcing_mesh_fileprint')
             read(buffer, *, iostat=ios) forcing_mesh_fileprint !print forcing array to file
          case ('closest_distance')
             read(buffer, *, iostat=ios) closest_distance !track minimum separation between points?
          case ('force_freq')
             read(buffer, *, iostat=ios) force_freq !forcing frequency
          case ('sticky_z_boundary')
             read(buffer, *, iostat=ios) sticky_z_boundary !particles stuck to top/bottom boundaries
          case ('phonon_emission')
             read(buffer, *, iostat=ios) phonon_emission !phonon emission on or off
          case ('phonon_percent')
             read(buffer, *, iostat=ios) phonon_percent !percentage of max curv we cutoff at
          case ('hyperviscosity')
             read(buffer, *, iostat=ios) hyperviscosity !use hyperviscosity to dissipate KWC
          case ('hyp_power')
             read(buffer, *, iostat=ios) hyp_power !use hyperviscosity to dissipate KWC
          case ('hyp_curv')
             read(buffer, *, iostat=ios) hyp_curv !use hyperviscosity to dissipate KWC
          case ('hyp_nu')
             read(buffer, *, iostat=ios) hyp_nu !use hyperviscosity to dissipate KWC
           case ('special_dump')
             read(buffer, *, iostat=ios) special_dump !special dump
          case ('quasi_pcount')
             read(buffer, *, iostat=ios) quasi_pcount !number of quasi particles
          case ('part_count')
             read(buffer, *, iostat=ios) part_count !number of particles
          case ('particle_type')
             read(buffer, *, iostat=ios) particle_type !particle type (fluid/intertial)
          case('particle_super_velocity')
             read(buffer, *, iostat=ios) particle_super_velocity !do particles feel superlfluid vel
          case ('particles_only')
             read(buffer, *, iostat=ios) particles_only !only evolve particles
          case ('part_stokes')
             read(buffer, *, iostat=ios) part_stokes !stokes number of inertial particles
          case ('tree_theta')
             read(buffer, *, iostat=ios) tree_theta !tree code, opening angle
          case ('tree_print')
             read(buffer, *, iostat=ios) tree_print !print the tree mesh
          case ('tree_extra_correction')
             read(buffer, *, iostat=ios) tree_extra_correction
          case ('anisotropy_params')
             read(buffer, *, iostat=ios) anisotropy_params
          case ('torus_q')
             read(buffer, *, iostat=ios) torus_q !q integer for torus knot initial condition
          case ('torus_p')
             read(buffer, *, iostat=ios) torus_p !p integer for torus knot initial condition
          case('criss_cross_bundle')
             read(buffer, *, iostat=ios) criss_cross_bundle !for criss-cross initial condition
          case('criss_cross_width')
             read(buffer, *, iostat=ios) criss_cross_width !for criss-cross initial condition
          case ('torus_epsilon')
             read(buffer, *, iostat=ios) torus_epsilon !epsilon real for torus knot initial condition
          case ('wave_count')
             read(buffer, *, iostat=ios) wave_count !for wave_spec initial conditions
          case ('wave_slope')
             read(buffer, *, iostat=ios) wave_slope !for wave_spec initial conditions
          case ('wave_start')
             read(buffer, *, iostat=ios) wave_start !for wave_spec initial conditions
          case ('wave_skip')
             read(buffer, *, iostat=ios) wave_skip !for wave_spec initial conditions
          case ('wave_amp')
             read(buffer, *, iostat=ios) wave_amp !for wave_spec initial conditions
          case ('wave_type')
             read(buffer, *, iostat=ios) wave_type !for wave_spec initial conditions
          case ('macro_ring_R')
             read(buffer, *, iostat=ios) macro_ring_R !for macro_ring initial conditions
          case ('macro_ring_a')
             read(buffer, *, iostat=ios) macro_ring_a !for macro_ring initial conditions
          case ('nf_mra_factor')
             read(buffer, *, iostat=ios) nf_mra_factor !for factor for minor radius in NF                 
          case ('hyperboloid_e')
             read(buffer, *, iostat=ios) hyperboloid_e !for hyperboloid initial conditions
          case ('hyperboloid_r')
             read(buffer, *, iostat=ios) hyperboloid_r !for hyperboloid initial conditions
          case ('curv_hist')
             read(buffer, *, iostat=ios) curv_hist !do we want binned curvature info?
          case ('torsion_hist')
             read(buffer, *, iostat=ios) torsion_hist !do we want binned torsion info?
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
          case ('KS_bubble')
             read(buffer, *, iostat=ios) KS_bubble !force wavenumbers apart
          case ('KS_rey_int')
             read(buffer, *, iostat=ios) KS_rey_int !KS Reynolds number proxy
          case ('KS_maximise_rey')
             read(buffer, *, iostat=ios) KS_maximise_rey !Force 1st wavenumber to be as small as possible
          case ('KS_modes')
             read(buffer, *, iostat=ios) KS_modes !the number of KS modes
          case ('recon_time_info')
             read(buffer, *, iostat=ios) recon_time_info !time diffn between recon
          case ('one_dim')
             read(buffer, *, iostat=ios) one_dim !size of 1D print
          case ('one_dim_direction')
             read(buffer, *, iostat=ios) one_dim_direction !which axis we perform 1D print    
          case ('one_dim_lattice')
             read(buffer, *, iostat=ios) one_dim_lattice !size of 1D lattice print
          case ('one_dim_lattice_count')
             read(buffer, *, iostat=ios) one_dim_lattice_count !number of lines in lattice       
          case ('two_dim')
             read(buffer, *, iostat=ios) two_dim !size of 2D print
          case ('recon_info')
             read(buffer, *, iostat=ios) recon_info !extra reconnection information
          case ('sep_inf')
             read(buffer, *, iostat=ios) sep_inf !point separation
          case ('line_sep_inf')
             read(buffer, *, iostat=ios) line_sep_inf !line separation
          case ('topo_inf')
             read(buffer, *, iostat=ios) topo_inf !topological information
          case ('energy_inf')
             read(buffer, *, iostat=ios) energy_inf !vortex energy - only open boundaries
          case ('particle_plane_inf')
             read(buffer, *, iostat=ios) particle_plane_inf !point vortex in z=0
          case ('switch_off_recon')
             read(buffer, *, iostat=ios) switch_off_recon !for test cases only!
          case ('seg_fault')
             read(buffer, *, iostat=ios) seg_fault !use print statements to find segmentation faults
          case ('simple_plots')
             read(buffer, *, iostat=ios) simple_plots !perform plots on the fly
          case ('NAN_test')
             read(buffer, *, iostat=ios) NAN_test !test for NANs
          case ('max_warn_count')
             read(buffer, *, iostat=ios) max_warn_count !maximum number of warning counts we survive
          case ('overide_timestep_check')
             read(buffer, *, iostat=ios) overide_timestep_check !no timestep check 
          case ('smoothing_length')
             read(buffer, *, iostat=ios) smoothing_length !length we are smoothing over (delta)
          case ('smoothing_interspace')
             read(buffer, *, iostat=ios) smoothing_interspace !use intervortex spacing as smoothing length?
          case ('sm_size')
             read(buffer, *, iostat=ios) sm_size !size of smoothing mesh
          case ('delta_adapt')
             read(buffer, *, iostat=ios) delta_adapt !is the discretisation adaptive
          case ('delta_adapt_print')
             read(buffer, *, iostat=ios) delta_adapt_print !print apative discretisation
          case ('inject_size')
             read(buffer, *, iostat=ios) inject_size !size of injected filaments
          case ('inject_skip')
             read(buffer, *, iostat=ios) inject_skip !how often we inject
          case ('inject_freq')
             read(buffer, *, iostat=ios) inject_freq !how often we switch injection
          case ('inject_type')
             read(buffer, *, iostat=ios) inject_type !how we inject new vortices
          case ('inject_stop')
             read(buffer, *, iostat=ios) inject_stop !when (if ever) we stop injecting 
          case ('full_loop_counter')
             read(buffer, *, iostat=ios) full_loop_counter !check loop count and size
          case ('batch_mode')
             read(buffer, *, iostat=ios) batch_mode !do we message user if there are messages?        
          case ('batch_name')
             read(buffer, *, iostat=ios) batch_name !what is the name of the run
          case ('batch_email')
             read(buffer, *, iostat=ios) batch_email !where do we message?
          case ('serial_run')
             read(buffer, *, iostat=ios) serial_run !where do we message?
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
    real ::  reload_normal_fluid_cutoff,reload_force_cutoff, reload_inject_stop
    logical :: reload_curv_hist,reload_torsion_hist, reload_vel_print, reload_vapor_print
    logical :: reload_topo_inf,reload_energy_inf,reload_recon_info,reload_sep_inf
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
          case ('force_cutoff')
             read(buffer, *, iostat=ios) reload_force_cutoff 
             if ((reload_force_cutoff>force_cutoff).or. & 
                (reload_force_cutoff<force_cutoff)) then
               force_cutoff=reload_force_cutoff
               write(*,'(a,f10.4)') 'RELOAD: changed force_cutoff to ', force_cutoff
             end if 
          case ('curv_hist')
             read(buffer, *, iostat=ios) reload_curv_hist 
             if (reload_curv_hist.neqv.curv_hist) then
               curv_hist=reload_curv_hist
               write(*,*) 'RELOAD: changed curv_hist to ', curv_hist
             end if
          case ('torsion_hist')
             read(buffer, *, iostat=ios) reload_torsion_hist 
             if (reload_torsion_hist.neqv.torsion_hist) then
               torsion_hist=reload_torsion_hist
               write(*,*) 'RELOAD: changed torsion_hist to ', torsion_hist
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
          case ('sep_inf')
             read(buffer, *, iostat=ios) reload_sep_inf
             if (reload_sep_inf.neqv.sep_inf) then
               sep_inf=reload_sep_inf
               write(*,*) 'RELOAD: changed sep_inf to ', sep_inf
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
    character(len=300) :: message_file
    write (*,*) '-------------------------FATAL ERROR-------------------------'
    write (*,*) trim(location) , ": " , trim(message)
    write (*,*) "FYI: t= ", t, "iteration number= ", itime 
    write (*,*) '-------------------------------------------------------------'
    if (batch_mode) then
      write(unit=message_file,fmt="(a,a,a,a,a,a,a)")'echo "','error message from ',trim(batch_name),&
      ': \n' ,trim(message), '" | mail -s qvort_error ',trim(batch_email)
      call system(message_file)
    end if
    stop
  end subroutine
  !*************************************************************************************************  
  !>if in batch mode send an email to say code has finished
  subroutine completion_message
    implicit none      
    character(len=300) :: message_file
    if (batch_mode) then
      write(unit=message_file,fmt="(a,a,a,a,a)")'echo "','completion notification from ',&
      trim(batch_name), '" | mail -s qvort_finished ',trim(batch_email)
      call system(message_file)
    end if
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
    !increment warning count
    warn_count=warn_count+1
    if (warn_count>max_warn_count) then
      call fatal_error('warning count','maximum number of warning counts exceeded')
    end if
  end subroutine
  !*************************************************************************************************  
  !>generate a new random seed, or read one in from ./data if restating code
  subroutine init_random_seed()
     !CREATE A NEW RANDOM SEED, UNLESS RESTARTING CODE
     integer :: i, n=16, clock
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
  !*************************************************
  subroutine init_openmp
    implicit none
    integer :: nthreads, thread_limit
    write(*,'(a)') ' ------------------------OPENMP-----------------------' 
    if (serial_run) then
      call omp_set_num_threads(1) 
      write(*,'(a)') 'serial run, running on one process'
    else
      nthreads=omp_get_max_threads()
      write(*,'(a,i2.2,a)') ' parallel run, running on ', nthreads, ' processes'
    end if
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
!!- \p sep_inf - print binned point separation information?
!!- \p vel_print - print velocity information to run angela's velocity stat scripts
!!- \p vapor_print - print meshes as binary files speciffically for vapor?
!!- \p one_dim - 1D mesh size (runs in x-direction, y=0,z=0)
!!- \p two_dim - 2D mesh size (in xy-plane, z=0)
!!- \p recon_info - print specific reconnection information
!!- \p topo_inf - calculate topological information of the filament (linking and writhe) very slow!
!!- \p energy_inf - print energy information to file 
!!- \p inject_stop - if (when) the code stops injecting loops (if inject has been set)
