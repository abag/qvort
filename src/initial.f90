!> Module which contains all the routines/call to routines which setup and
!!restart the code. Initial conditions for filament set in run.in for more 
!!information on the options see \ref INIT
module initial
  use cdata
  use normal_fluid
  use forcing
  use periodic
  use smoothing
    use sph
  contains
  !*************************************************************************
  !>Prints and sets up intial conditions - will give warnings/errors if there
  !>are any conflicting options set in run.in
  subroutine init_setup()
    use quasip
    use particles
    use mag !not strictly needed as mag used in timestep.mod used in quasip.mod 
    use sph
    implicit none
    logical :: restart
    write(*,'(a)') ' ---------------------VORTEX PARAMETERS------------------' 
    write(*,'(a,f9.7)') ' quantum of circulation is:', quant_circ 
    write(*,'(a,e9.3)') ' core size is:', corea 
    !check particle separation has been set
    if (delta<epsilon(0.)) call fatal_error('init.mod','delta must be set in run.in')
    !check that the particle count (pcount) has been set
    if (particles_only) then
      if (init_pcount>0) call fatal_error('init.mod','must have pcount=0 for particles_only')
    else
      if (init_pcount<5) then
        call fatal_error('init.mod','you must set enough intial particles')
      end if 
    end if
    !we must check the timestep is sufficient to resolve the motion
    !based on the smallest separation possible in the code
    write(*,'(a)') ' ---------------------TIME-STEP--------------------' 
    call timestep_check !initial.mod
    if (dt_adapt) then
      write(*,'(a)') ' using an adpative timestepping routine - in testing!'
    end if
    if (phonon_emission) then
      write(*,'(a)') ' ------------------PHONON EMISSION--------------------' 
      write(*,'(a,f6.2,a,f8.1)') ' simulating phonon emission, cutoff is ', 100*phonon_percent, '% of max:', 2/delta
    end if
    !how is data being outputted (binary or formatted)
    if (binary_print) then
    write(*,'(a)') ' ---------------------DATA FORMAT--------------------' 
      write(*,*) 'binary data output, formatted data can be selected in run.in'
    else
      write(*,*) 'formatted data output selected in run.in'
    end if
    write(*,'(a,i3.2,a)') ' outputting filament information every ', shots, ' time-steps'
    !periodic bounday conditions?
    write(*,'(a)') ' ---------------------BOUNDARY CONDITIONS--------------------' 
    if (box_size>0.) then
      !what is the boundary
      select case(boundary)
        case('periodic')
          periodic_bc=.true.
          write(*,'(a,f6.3)') ' running with periodic boundaries, box size:', box_size
        case('mirror')
          mirror_bc=.true.
          call warning_message('init.mod','mirror b.c.s are still in testing and will probably fail at some point in the run!')
          select case(velocity)
            case('LIA','Tree')
              call fatal_error('init_setup','mirror bcs are not set up to work with the LIA/Tree veloctity')
          end select
          write(*,'(a,f6.3)') ' running with mirrored boundaries, box size:', box_size
          if (mirror_print) write(*,*) 'printing mirror filaments to file'
        case('open')
          write(*,*) 'running with open boundaries'
        case default
          call fatal_error('init_setup:', 'incorrect boundary parameter')
      end select
    else
      call fatal_error('init_setup:', 'box size is less than zero')
    end if
    write(*,'(a)') ' --------------------NORMAL FLUID--------------------' 
    call setup_normal_fluid !normal_fluid.mod
    !sort out if there is a special dump time
    int_special_dump=int(special_dump/dt) !convert to integer
    write(*,'(a)') ' --------------------INITIAL CONDITIONS--------------------' 
    !check if we can restart the code
    inquire(file="./data/var.dat", exist=restart)
    if (restart) then
      call data_restore !init.mod
    else  
      pcount=init_pcount
      allocate(f(pcount)) !main vector allocated
      !choose the correct setup routine based on the value of initf in run.in
      select case(initf)
        case('single_loop')
          call setup_single_loop !init.mod
        case('single_line')
          call setup_single_line !init.mod
        case('random_loops')
          call setup_random_loops !init.mod
        case('crow')
          call setup_crow !init.mod
        case('leap-frog')
          call setup_leap_frog !init.mod
        case('linked_filaments')
          call setup_linked_filaments !init.mod
        case('colliding_loops')
          call setup_colliding_loops !init.mod
        case('kivotedes')
          call setup_kivotedes !init.mod
        case('cardoid')
          call setup_cardoid !init.mod
        case('wave_loop')
          call setup_wave_loop !init.mod
        case('linked_wave_loop')
          call setup_linked_wave_loop !init.mod
        case('wave_line')
          call setup_wave_line !init.mod
        case('tangle')
          call setup_tangle !init.mod
        case('criss-cross')
          call setup_criss_cross !init.mod
        case default
          call fatal_error('cdata.mod:init_setup', &
                         'invalid choice for initf parameter') !cdata.mod
      end select
    end if
    !test if we have a non-zero mesh size
    write(*,'(a)') ' ------------------------MESH-----------------------' 
    if (mesh_size>0) then
      if (box_size>0) then
        call setup_mesh !init.mod
      else 
        call fatal_error('cdata.mod:init_setup', &
        'running with a non-zero mesh size requires a non-zero box size') !cdata.mod
      end if
    else
      write(*,*) 'velocity fields not being stored on a mesh - (no spectra etc.)'
    end if
    if (mesh_shots<shots) call warning_message('init.mod',&
                          'mesh shots < shots which will create output anomolies')
    !do we employ forcing on the boundary?
    write(*,'(a)') ' -----------------------FORCING-----------------------' 
    call setup_forcing !forcing,mod
    !are there quasi particles in the code?
    write(*,'(a)') ' ------------------------QUASI PARTICLES-----------------------' 
    if (quasi_pcount>0) then
      if (restart) then
        write(*,*) 'quasi particles have been restored from dump file, see above'
      else 
        call setup_quasip !quasip.mod
      end if
    else
      write(*,*) 'no quasi particles in the code'
    end if
    !are there particles in the code?
    write(*,'(a)') ' ------------------------PARTICLES-----------------------' 
    if (part_count>0) then
      if (restart) then
        write(*,*) 'particles have been restored from dump file, see above'
      else 
        call setup_particles !particles.mod
      end if
    else
      if (particles_only) then
        if (SPH_count==0) then 
          call fatal_error('init.mod','must have part_count>0 &
                            or SPH_count>0 for particles_only')
        end if
      end if
      write(*,*) 'no particles in the code'
    end if
    write(*,'(a)') ' ---------------------VELOCITY CALCULATION----------------------' 
    !is the tree code being used?
    if (tree_theta>0) then
      if (box_size>0.) then
        write(*,*) 'using tree algorithms for reconnection routine - scales like O(NlogN)'
      else
        call fatal_error('cdata.mod:init_setup','tree algorithms require a positive box size')
      end if
    else
      write(*,*) 'using brute force reconnection routine - scales like O(N^2)'
    end if
    !print information about the velocity field to screen
    select case(velocity)
      case('Off')
        select case(normal_velocity)
          case('zero')
            select case(force)
              case ('off')
                call warning_message('initial.mod',&
                'no superfluid velocity, no normal velocity, no forcing')
            end select
        end select
        write(*,*) 'No superfluid velocity: simulation of passive lines'
      case('LIA')
        write(*,*) 'using local induction approximation - scales like O(N)'
      case('BS')
        write(*,*) 'using full Biot-Savart integral - scales like O(N^2)'
      case('Tree')
        if (tree_theta<epsilon(0.)) then 
          call fatal_error('init.mod:init_setup', & 
          'runnning with tree velocity but tree_theta is zero') !cdata.mod
        end if
        write(*,*) 'using tree approximation to Biot-Savart integral - scales like O(NlogN)'
        if (tree_theta>0.6) then
          call warning_message('init_setup','tree_theta>0.6, the vortices may be unstable')
        end if
     case default
       print*, 'correct value for velocity in run.in has not been set'
       print*, 'options are: LIA, BS, Tree'
       call fatal_error('init.mod:init_setup', & 
        'correct value for "velocity" in run.in has not been set') !cdata.mod
    end select
    !any special diagnostic information?
    write(*,'(a)') ' ---------------------FURTHER DIAGNOSTICS----------------------' 
    if (curv_hist) write(*,*) 'printing histograms of curvature to file'
    if (vel_print) write(*,'(a,i4.2,a)') ' printing full velocity information every: ', mesh_shots, ' timesteps'
    !final boundary conditions sanity check
    if (periodic_bc.and.mirror_bc) call fatal_error('init.mod','both periodic and mirror bcs are set')
    if (one_dim>0) write(*,'(a,i5.3)') ' printing 1D velocity info to file, mesh size: ', one_dim
    if (two_dim>0) write(*,'(a,i5.3)') ' printing 2D velocity info to file, mesh size: ', two_dim
    if (recon_info) write(*,*) 'printing extra reconnection information to file'
    if (switch_off_recon) call warning_message('init.mod','reconnections switched off: I HOPE YOU KNOW WHAT YOUR DOING!')
    !----------------------gaussian smoothing of field------------------------------
    if (sm_size>0) call setup_smoothing_mesh !smoothing.mod
    !----------------------------magnetic field-------------------------------------
    if (magnetic) call setup_mag !mag.mod
    !------------------------------SPH----------------------------------
    if (SPH_count>0) then
      call setup_SPH !sph.mod
    else
      !must check for printing to screen
      if (particles_only) then
        if (part_count==0) then 
          call fatal_error('init.mod','must have part_count>0 &
                            or SPH_count>0 for particles_only')
        end if
      end if
    end if
  end subroutine
  !**********************************************************************
  !>restart the code code periodically writes all the main variables to a file
  !>./data/dump.dat, if this is detected at startup then the code will restart
  !>from that data file
  subroutine data_restore
    !restart the code
    use stiff_solver
    implicit none
    integer :: dummy_itime 
    open(unit=63,file="./data/var.dat",FORM='unformatted')
      read(63) pcount
      read(63) recon_count
      read(63) dummy_itime
      read(63) t
      allocate(f(pcount))
      read(63) f
      write(*,*) 'restored vortex filament'
      read(63) quasi_pcount
      if (quasi_pcount>0) then
        allocate(g(quasi_pcount))
        read(63) g
        write(*,'(a,i4.1,a)') ' restored ', quasi_pcount,' quasi particles'
        !finally intialise the backwards difference coefficients array
        call set_BDF_coeff !stiff_solver.mod
      end if
      read(63) part_count
      if (part_count>0) then
        allocate(p(part_count))
        read(63) p
        write(*,'(a,i4.1,a)') ' restored ', part_count,' particles'
      end if
    close(63)
    nstart=dummy_itime+1
    write(*,*) 'data read in from dump file at t=', t
    if (quasi_pcount>0) then
      write(*,*)'note that the quasi particles will not be restored'
    end if
  end subroutine
  !*************************************************************************
  !>set up a single loop in the x-y plane, it's size is dictated by the initial 
  !>number of particles set and the size of \f$\delta\f$
  subroutine setup_single_loop
    implicit none
    real :: velocity
    real :: radius
    integer :: i 
    radius=(0.75*pcount*delta)/(2*pi) !75% of potential size 
    velocity=(quant_circ/(4*pi*radius))*log(8E8*radius)
    write(*,*) 'initf: single loop, radius of loop:', radius
    write(*,*) 'velocity should be:', velocity
    
    !loop over particles setting spatial and 'loop' position
    do i=1, pcount
      f(i)%x(1)=radius*sin(pi*real(2*i-1)/pcount)
      f(i)%x(2)=radius*cos(pi*real(2*i-1)/pcount)
      f(i)%x(3)=0.
      if (i==1) then
        f(i)%behind=pcount ; f(i)%infront=i+1
      else if (i==pcount) then 
        f(i)%behind=i-1 ; f(i)%infront=1
      else
        f(i)%behind=i-1 ; f(i)%infront=i+1
      end if
      !zero the stored velocities
      f(i)%u1=0. ; f(i)%u2=0.
    end do   
  end subroutine
  !****************************************************************************
  !>set up a single line from -z to z, number of particles is automatically adjusted
  !>to box size and \f$\delta\f$
  subroutine setup_single_line
    implicit none
    integer :: pcount_required
    integer :: i
    if (periodic_bc) then
      !work out the number of particles required for single line
      !given the box size specified in run.i
      pcount_required=nint(box_size/(0.75*delta)) !75%
      write(*,*) 'changing size of pcount to fit with box_length and delta'
      write(*,*) 'pcount is now', pcount_required
      deallocate(f) ; pcount=pcount_required ; allocate(f(pcount))
    else
      call fatal_error('init.mod:setup_single_line', &
      'periodic boundary conditions required')
    end if
    do i=1, pcount
      f(i)%x(1)=0.
      f(i)%x(2)=0.
      f(i)%x(3)=-box_size/2.+box_size*real(2*i-1)/(2.*pcount)
      if (i==1) then
        f(i)%behind=pcount ; f(i)%infront=i+1
      else if (i==pcount) then 
        f(i)%behind=i-1 ; f(i)%infront=1
      else
        f(i)%behind=i-1 ; f(i)%infront=i+1
      end if
      !zero the stored velocities
      f(i)%u1=0. ; f(i)%u2=0.
    end do
  end subroutine
  !*************************************************************************
  !>draw cardoid (http://en.wikipedia.org/wiki/Cardioid) in the x-y plane
  subroutine setup_cardoid
    implicit none
    real :: velocity
    real :: radius
    integer :: i 
    radius=(0.75*pcount*delta)/(6*pi) !75% of potential size 
    write(*,*) 'initf: cardoid, cusp to probe kelvin waves'
    !loop over particles setting spatial and 'loop' position
    do i=1, pcount
      f(i)%x(1)=radius*(2*sin(pi*real(2*i-1)/pcount)-sin(2*pi*real(2*i-1)/pcount))
      f(i)%x(2)=radius*(2*cos(pi*real(2*i-1)/pcount)-cos(2*pi*real(2*i-1)/pcount))
      f(i)%x(3)=0.0
      if (i==1) then
        f(i)%behind=pcount ; f(i)%infront=i+1
      else if (i==pcount) then 
        f(i)%behind=i-1 ; f(i)%infront=1
      else
        f(i)%behind=i-1 ; f(i)%infront=i+1
      end if
      !zero the stored velocities
      f(i)%u1=0. ; f(i)%u2=0.
    end do   
  end subroutine
  !*************************************************************************
  !>two loops in the x-y plane their sizes are dictated by the initial 
  !>number of particles set and the size of \f$\delta\f$
  subroutine setup_leap_frog
    implicit none
    real :: radius
    integer :: i 
    if (mod(pcount,2)/=0) then
      call fatal_error('init.mod:setup_leap_frog', &
      'pcount is not a multiple of 2-aborting')
    end if
    radius=(0.75*pcount*delta)/(4*pi) !75% of potential size
    write(*,*) 'initf: leap-frog, radius of loops:', radius 
    !loop over particles setting spatial and 'loop' position
    do i=1, pcount/2
      f(i)%x(1)=radius*sin(pi*real(2*i-1)/(pcount/2))
      f(i)%x(2)=radius*cos(pi*real(2*i-1)/(pcount/2)) ; f(i)%x(3)=0.
      if (i==1) then
        f(i)%behind=pcount/2 ; f(i)%infront=i+1
      else if (i==pcount/2) then 
        f(i)%behind=i-1 ; f(i)%infront=1
      else
        f(i)%behind=i-1 ; f(i)%infront=i+1
      end if
      !zero the stored velocities
      f(i)%u1=0. ; f(i)%u2=0.
    end do
    !second loop
    do i=pcount/2+1, pcount
      f(i)%x(1)=radius*sin(pi*real(2.*(i-pcount/2)-1)/(pcount/2))
      f(i)%x(2)=radius*cos(pi*real(2.*(i-pcount/2)-1)/(pcount/2)) 
      f(i)%x(3)=2.*delta
      if (i==(pcount/2+1)) then
        f(i)%behind=pcount ; f(i)%infront=i+1
      else if (i==pcount) then 
        f(i)%behind=i-1 ; f(i)%infront=1+pcount/2
      else
        f(i)%behind=i-1 ; f(i)%infront=i+1
      end if
      !zero the stored velocities
      f(i)%u1=0. ; f(i)%u2=0.
    end do    
  end subroutine
  !*************************************************************************
  !>anti parallel lines from -z to z which drives the crow instability
  !>particle count automatically adjusted
  subroutine setup_crow
    implicit none
    integer :: pcount_required
    integer :: i 
    if (periodic_bc) then
      !work out the number of particles required for single line
      !given the box size specified in run.i
      pcount_required=2*nint(box_size/(0.75*delta)) !75%
      print*, 'changing size of pcount to fit with box_length and delta'
      print*, 'pcount is now', pcount_required
      deallocate(f) ; pcount=pcount_required ; allocate(f(pcount))
    else
      call fatal_error('init.mod:setup_crow', &
      'periodic boundary conditions required')
    end if
    write(*,*) 'initf: crow, separation of lines is:', 2.*delta 
    !loop over particles setting spatial and 'loop' position
    do i=1, pcount/2
      f(i)%x(1)=0.
      f(i)%x(2)=delta
      f(i)%x(3)=-box_size/2.+box_size*real(2*i-1)/(pcount)
      if (i==1) then
        f(i)%behind=pcount/2 ; f(i)%infront=i+1
      else if (i==pcount/2) then 
        f(i)%behind=i-1 ; f(i)%infront=1
      else
        f(i)%behind=i-1 ; f(i)%infront=i+1
      end if
      !zero the stored velocities
      f(i)%u1=0. ; f(i)%u2=0.
    end do
    !second line
    do i=pcount/2+1, pcount
      f(i)%x(1)=0.
      f(i)%x(2)=-delta
      f(i)%x(3)=box_size/2.-box_size*real(2*(i-pcount/2)-1)/(pcount)
      if (i==(pcount/2+1)) then
        f(i)%behind=pcount ; f(i)%infront=i+1
      else if (i==pcount) then 
        f(i)%behind=i-1 ; f(i)%infront=1+pcount/2
      else
        f(i)%behind=i-1 ; f(i)%infront=i+1
      end if
      !zero the stored velocities
      f(i)%u1=0. ; f(i)%u2=0.
    end do    
  end subroutine
  !*************************************************************************
  !>two linked loops which drive a reconnection
  subroutine setup_linked_filaments
    !two linked loops 
    implicit none
    real :: radius
    integer :: i 
    if (mod(pcount,2)/=0) then
      call fatal_error('init.mod:setup_linked_filaments', &
      'pcount is not a multiple of 2-aborting')
    end if
    radius=(0.75*pcount*delta)/(4*pi) !75% of potential size
    write(*,'(a,f9.6)') ' initf: linked filaments, radius of loops:', radius 
    !loop over particles setting spatial and 'loop' position
    do i=1, pcount/2
      f(i)%x(1)=radius*sin(pi*real(2*i-1)/(pcount/2))
      f(i)%x(2)=radius*cos(pi*real(2*i-1)/(pcount/2))
      f(i)%x(3)=0.
      if (i==1) then
        f(i)%behind=pcount/2 ; f(i)%infront=i+1
      else if (i==pcount/2) then 
        f(i)%behind=i-1 ; f(i)%infront=1
      else
        f(i)%behind=i-1 ; f(i)%infront=i+1
      end if
      !zero the stored velocities
      f(i)%u1=0. ; f(i)%u2=0.
    end do
    !second loop
    do i=pcount/2+1, pcount
      f(i)%x(1)=radius*sin(pi*real(2.*(i-pcount/2)-1)/(pcount/2))+radius/0.51
      f(i)%x(2)=0.
      f(i)%x(3)=radius*cos(pi*real(2.*(i-pcount/2)-1)/(pcount/2)) 

      if (i==(pcount/2+1)) then
        f(i)%behind=pcount ; f(i)%infront=i+1
      else if (i==pcount) then 
        f(i)%behind=i-1 ; f(i)%infront=1+pcount/2
      else
        f(i)%behind=i-1 ; f(i)%infront=i+1
      end if
      !zero the stored velocities
      f(i)%u1=0. ; f(i)%u2=0.
    end do    
  end subroutine
  !*************************************************************************
  !>two unlinked loops which move towards each other driving a reconnection
  subroutine setup_colliding_loops
    !two linked loops 
    implicit none
    real :: radius
    integer :: i 
    if (mod(pcount,2)/=0) then
      call fatal_error('init.mod:setup_linked_filaments', &
      'pcount is not a multiple of 2-aborting')
    end if
    radius=(0.75*pcount*delta)/(4*pi) !75% of potential size
    write(*,'(a,f9.6)') ' initf: colliding loops, radius of loops:', radius 
    !loop over particles setting spatial and 'loop' position
    do i=1, pcount/2
      f(i)%x(1)=radius*sin(pi*real(2*i-1)/(pcount/2))
      f(i)%x(2)=radius*cos(pi*real(2*i-1)/(pcount/2))
      f(i)%x(3)=box_size/4.
      if (i==1) then
        f(i)%behind=pcount/2 ; f(i)%infront=i+1
      else if (i==pcount/2) then 
        f(i)%behind=i-1 ; f(i)%infront=1
      else
        f(i)%behind=i-1 ; f(i)%infront=i+1
      end if
      !zero the stored velocities
      f(i)%u1=0. ; f(i)%u2=0.
    end do
    !second loop
    do i=pcount/2+1, pcount
      f(i)%x(1)=radius*sin(pi*real(2*i-1)/(pcount/2))
      f(i)%x(2)=-radius*cos(pi*real(2*i-1)/(pcount/2))
      f(i)%x(3)=-box_size/4.

      if (i==(pcount/2+1)) then
        f(i)%behind=pcount ; f(i)%infront=i+1
      else if (i==pcount) then 
        f(i)%behind=i-1 ; f(i)%infront=1+pcount/2
      else
        f(i)%behind=i-1 ; f(i)%infront=i+1
      end if
      !zero the stored velocities
      f(i)%u1=0. ; f(i)%u2=0.
    end do    
  end subroutine
  !*************************************************************************
  !>four unlinked loops as in Kivotedes 2001, PRL
  subroutine setup_kivotedes
    implicit none
    real :: radius
    integer :: i 
    if (mod(pcount,4)/=0) then
      call fatal_error('init.mod:setup_linked_filaments', &
      'pcount is not a multiple of 4-aborting')
    end if
    radius=(0.75*pcount*delta)/(8*pi) !75% of potential size
    write(*,'(a,f9.6)') ' initf: 4 unlinked loops as in kivotedes, radius of loops:', radius 
    !loop over particles setting spatial and 'loop' position
    do i=1, pcount/4
      f(i)%x(1)=-radius*sin(pi*real(2*i-1)/(pcount/4))
      f(i)%x(2)=radius*cos(pi*real(2*i-1)/(pcount/4))
      f(i)%x(3)=0.
      if (i==1) then
        f(i)%behind=pcount/4 ; f(i)%infront=i+1
      else if (i==pcount/4) then 
        f(i)%behind=i-1 ; f(i)%infront=1
      else
        f(i)%behind=i-1 ; f(i)%infront=i+1
      end if
      !zero the stored velocities
      f(i)%u1=0. ; f(i)%u2=0.
    end do
    !second loop
    do i=pcount/4+1, 2*pcount/4
      f(i)%x(1)=-radius*sin(pi*real(2.*(i-pcount/4)-1)/(pcount/4))
      f(i)%x(2)=radius/0.9
      f(i)%x(3)=radius*cos(pi*real(2.*(i-pcount/4)-1)/(pcount/4))+radius/0.9

      if (i==(pcount/4+1)) then
        f(i)%behind=2*pcount/4 ; f(i)%infront=i+1
      else if (i==2*pcount/4) then 
        f(i)%behind=i-1 ; f(i)%infront=1+pcount/4
      else
        f(i)%behind=i-1 ; f(i)%infront=i+1
      end if
      !zero the stored velocities
      f(i)%u1=0. ; f(i)%u2=0.
    end do    
    !third loop
    do i=2*pcount/4+1, 3*pcount/4
      f(i)%x(1)=radius*sin(pi*real(2.*(i-2*pcount/4)-1)/(pcount/4))
      f(i)%x(2)=radius*cos(pi*real(2.*(i-2*pcount/4)-1)/(pcount/4))
      f(i)%x(3)=2*radius/0.9

      if (i==(2*pcount/4+1)) then
        f(i)%behind=3*pcount/4 ; f(i)%infront=i+1
      else if (i==3*pcount/4) then 
        f(i)%behind=i-1 ; f(i)%infront=1+2*pcount/4
      else
        f(i)%behind=i-1 ; f(i)%infront=i+1
      end if
      !zero the stored velocities
      f(i)%u1=0. ; f(i)%u2=0.
    end do
    !final loop
    do i=3*pcount/4+1, pcount
      f(i)%x(1)=radius*sin(pi*real(2.*(i-3*pcount/4)-1)/(pcount/4))
      f(i)%x(2)=-radius/0.9
      f(i)%x(3)=radius*cos(pi*real(2.*(i-3*pcount/4)-1)/(pcount/4))+radius/0.9

      if (i==(3*pcount/4+1)) then
        f(i)%behind=pcount ; f(i)%infront=i+1
      else if (i==pcount) then 
        f(i)%behind=i-1 ; f(i)%infront=1+3*pcount/4
      else
        f(i)%behind=i-1 ; f(i)%infront=i+1
      end if
      !zero the stored velocities
      f(i)%u1=0. ; f(i)%u2=0.
    end do
  end subroutine
  !*************************************************************************
  !>clustered loops, helical/planar waves added with
  !>specific spectrum all set in run.in
  subroutine setup_wave_loop
    implicit none
    real :: wave_number, prefactor
    real :: amp, random_shift
    real :: zpos
    real :: loop_radius
    integer :: pcount_required
    integer :: loop_size, loop_position
    integer :: i, j, k
    !test run.in parameters, if wrong program will exit
    if (line_count==0) then
      call fatal_error('init.mod:setup_wave_loops', &
      'you have not set a value for line_count in run.in')
    end if
    if (mod(pcount,line_count)/=0) then
      call fatal_error('init.mod:setup_wave_loops', &
      'pcount/line_count is not an integer')
    end if
    write(*,'(a,i3.1,a)') ' drawing', line_count, ' loops in the x-y plane'
    write(*,'(i4.1,a,a,a,f9.5)') wave_count, ' ',trim(wave_type),' wave pertubations, with spectral slope:', wave_slope
    loop_size=int(pcount/line_count)
    loop_radius=loop_size*(0.75*delta)/(2*pi) !75% of potential size
    !START THE LOOP
    do i=1, line_count
      call random_number(zpos)
      zpos=(zpos-.5)*box_size*0.1 !1/10th of the box
      do j=1, loop_size
        loop_position=j+(i-1)*loop_size
        f(loop_position)%x(1)=loop_radius*sin(pi*real(2*j-1)/loop_size)
        f(loop_position)%x(2)=loop_radius*cos(pi*real(2*j-1)/loop_size)
        f(loop_position)%x(3)=zpos
        if(j==1) then
          f(loop_position)%behind=i*loop_size
          f(loop_position)%infront=loop_position+1
        else if (j==loop_size) then
          f(loop_position)%behind=loop_position-1
          f(loop_position)%infront=(i-1)*loop_size+1
        else
          f(loop_position)%behind=loop_position-1
          f(loop_position)%infront=loop_position+1
        end if
        f(loop_position)%u1=0. ; f(loop_position)%u2=0.
      end do
      !we have now drawn the basic loop, now we add the wave pertubations
      prefactor=wave_amp/(2.**wave_slope) !our starting wavenumber is 2
      if (i==1) then !only write this once
        write(*,*)'wave information recorded in ./data/wave_info.log'
        open(unit=34,file='./data/wave_info.log')
        write(34,*) '%------------------WAVE INFORMATION-------------------'
      end if
      do k=1, wave_count !wave_count set in run.in
        call ghostp !we must call this routine at the start of adding every wave 
                    !the routines normalf, binormalf rely on correct ghostpoints
        !on a loop so wavenumbers must be an integer
        wave_number=2+k !starting wavenumber is 2
        amp=prefactor*(wave_number**wave_slope)
        call random_number(random_shift) !help things along with a 
        random_shift=random_shift*2*pi   !random shift \in (0,2\pi)
        if (k==1) then
          write(34,'(a,i3.1)') '%loop number ', i
        end if
        write(34,'(f9.5,f9.5,f9.5)') wave_number, amp, random_shift
        do j=1, loop_size !reloop over particles
          loop_position=j+(i-1)*loop_size
          select case(wave_type) !what wave type do we want?
            case('planar')       !set this in run.in
              f(loop_position)%x(:)=f(loop_position)%x(:)+&
              binormalf(loop_position)*amp*delta*sin(random_shift+wave_number*2.*pi*real(2.*j-1)/(2.*loop_size))
            case('helical')
              f(loop_position)%x(:)=f(loop_position)%x(:)+&
              normalf(loop_position)*amp*delta*sin(random_shift+wave_number*2.*pi*real(2.*j-1)/(2.*loop_size))+&
              binormalf(loop_position)*amp*delta*cos(random_shift+wave_number*2.*pi*real(2.*j-1)/(2.*loop_size))
            case default
              call fatal_error('initial.mod:wave_line','incorrect wave type parameter')
          end select
        end do
      end do !close k loop
    end do !closes the i loop
    close(34)
  end subroutine
  !*************************************************************************
  !>linked loops, helical/planar waves added with specific spectrum all set in run.in
  subroutine setup_linked_wave_loop
    implicit none
    real :: wave_number, prefactor
    real :: amp, random_shift
    real :: loop_radius
    integer :: pcount_required
    integer :: loop_size, loop_position
    integer :: i, j, k
    !test run.in parameters, if wrong program will exit
    if (mod(pcount,2)/=0) then
      call fatal_error('init.mod:setup_linked_wave_loops', &
      'pcount/2 is not an integer')
    end if
    write(*,'(a)') ' drawing 2 linked loops'
    write(*,'(i4.1,a,a,a,f9.5)') wave_count, ' ',trim(wave_type),' wave pertubations, with spectral slope:', wave_slope
    loop_size=int(pcount/2)
    loop_radius=loop_size*(0.75*delta)/(2*pi) !75% of potential size
    !START THE LOOP
    do i=1, 2
      do j=1, loop_size
        if (i==1) then
          loop_position=j+(i-1)*loop_size
          f(loop_position)%x(1)=loop_radius*sin(pi*real(2*j-1)/loop_size)
          f(loop_position)%x(2)=loop_radius*cos(pi*real(2*j-1)/loop_size)
          f(loop_position)%x(3)=0.
          if(j==1) then
            f(loop_position)%behind=i*loop_size
            f(loop_position)%infront=loop_position+1
          else if (j==loop_size) then
            f(loop_position)%behind=loop_position-1
            f(loop_position)%infront=(i-1)*loop_size+1
          else
            f(loop_position)%behind=loop_position-1
            f(loop_position)%infront=loop_position+1
          end if
          f(loop_position)%u1=0. ; f(loop_position)%u2=0.
        else
          loop_position=j+(i-1)*loop_size
          f(loop_position)%x(1)=loop_radius*sin(pi*real(2*j-1)/loop_size)+loop_radius/0.5
          f(loop_position)%x(2)=0.
          f(loop_position)%x(3)=loop_radius*cos(pi*real(2*j-1)/loop_size)
          if(j==1) then
            f(loop_position)%behind=i*loop_size
            f(loop_position)%infront=loop_position+1
          else if (j==loop_size) then
            f(loop_position)%behind=loop_position-1
            f(loop_position)%infront=(i-1)*loop_size+1
          else
            f(loop_position)%behind=loop_position-1
            f(loop_position)%infront=loop_position+1
          end if
          f(loop_position)%u1=0. ; f(loop_position)%u2=0.
        end if
      end do
      !we have now drawn the basic loop, now we add the wave pertubations
      prefactor=wave_amp/(2.**wave_slope) !our starting wavenumber is 2
      if (i==1) then !only write this once
        write(*,*)'wave information recorded in ./data/wave_info.log'
        open(unit=34,file='./data/wave_info.log')
        write(34,*) '%------------------WAVE INFORMATION-------------------'
      end if
      do k=1, wave_count !wave_count set in run.in
        call ghostp !we must call this routine at the start of adding every wave 
                    !the routines normalf, binormalf rely on correct ghostpoints
        !on a loop so wavenumbers must be an integer
        wave_number=2+.25*k !starting wavenumber is 2
        amp=prefactor*(wave_number**wave_slope)
        call random_number(random_shift) !help things along with a 
        random_shift=random_shift*2*pi   !random shift \in (0,2\pi)
        if (k==1) then
          write(34,'(a,i3.1)') '%loop number ', i
        end if
        write(34,'(f9.5,f9.5,f9.5)') wave_number, amp, random_shift
        do j=1, loop_size !reloop over particles
          loop_position=j+(i-1)*loop_size
          select case(wave_type) !what wave type do we want?
            case('planar')       !set this in run.in
              f(loop_position)%x(:)=f(loop_position)%x(:)+&
              binormalf(loop_position)*amp*delta*sin(random_shift+wave_number*2.*pi*real(2.*j-1)/(2.*loop_size))
            case('helical')
              f(loop_position)%x(:)=f(loop_position)%x(:)+&
              normalf(loop_position)*amp*delta*sin(random_shift+wave_number*2.*pi*real(2.*j-1)/(2.*loop_size))+&
              binormalf(loop_position)*amp*delta*cos(random_shift+wave_number*2.*pi*real(2.*j-1)/(2.*loop_size))
            case default
              call fatal_error('initial.mod:wave_line','incorrect wave type parameter')
          end select
        end do
      end do !close k loop
    end do !closes the i loop
    close(34)
  end subroutine
  !*************************************************************************
  !>lines from the lop of the box to the bottom, helical/planar waves added with
  !>specific spectrum all set in run.in
  subroutine setup_wave_line
    implicit none
    real :: wave_number, prefactor
    real :: amp, random_shift
    real :: xpos, ypos
    integer :: pcount_required
    integer :: line_size, line_position
    integer :: i, j, k
    !test run.in parameters, if wrong program will exit
    if (line_count==0) then
      call fatal_error('init.mod:setup_wave_line', &
      'you have not set a value for line_count in run.in')
    end if
    if (periodic_bc) then
      !work out the number of particles required for single line
      !given the box size specified in run.i
      pcount_required=line_count*nint(box_size/delta) !100% as waves are added
      write(*,*) 'changing size of pcount to fit with box_length and delta'
      write(*,'(a,i7.1)') ' pcount is now:', pcount_required
      deallocate(f) ; pcount=pcount_required ; allocate(f(pcount))
    else
      call fatal_error('init.mod:setup_wave_line', &
      'periodic boundary conditions required')
    end if
    write(*,'(a,i3.1,a)') ' drawing', line_count, ' lines from -z to +z'
    write(*,'(i4.1,a,a,a,f9.5)') wave_count, ' ',trim(wave_type),' wave pertubations, with spectral slope:', wave_slope
    line_size=int(pcount/line_count)
    !START THE LOOP
    do i=1, line_count
      call random_number(xpos) ; call random_number(ypos)
      xpos=(xpos-.5)*box_size*0.1 !1/10th of the box
      ypos=(ypos-.5)*box_size*0.1 
      do j=1, line_size
        line_position=j+(i-1)*line_size
        f(line_position)%x(1)=xpos
        f(line_position)%x(2)=ypos
        f(line_position)%x(3)=-box_size/2.+box_size*real(2*j-1)/(2.*line_size)
        if(j==1) then
          f(line_position)%behind=i*line_size
          f(line_position)%infront=line_position+1
        else if (j==line_size) then
          f(line_position)%behind=line_position-1
          f(line_position)%infront=(i-1)*line_size+1
        else
          f(line_position)%behind=line_position-1
          f(line_position)%infront=line_position+1
        end if
        f(line_position)%u1=0. ; f(line_position)%u2=0.
      end do
      !we have now drawn the basic line, now we add the wave pertubations
      prefactor=wave_amp/(2.**wave_slope) !our starting wavenumber is 2
      if (i==1) then !only write this once
        write(*,*)'wave information recorded in ./data/wave_info.log'
        open(unit=34,file='./data/wave_info.log')
        write(34,*) '%------------------WAVE INFORMATION-------------------'
      end if
      do k=1, wave_count !wave_count set in run.in
        !wave_number=2+.1*k !starting wavenumber is 2, step in .1
        wave_number=2+2*k !starting wavenumber is 2, step in .1
        amp=prefactor*(wave_number**wave_slope)
        call random_number(random_shift) !help things along with a 
        random_shift=random_shift*2*pi   !random shift \in (0,2\pi)
        if (i==1) then
          write(34,'(f9.5,f9.5,f9.5)') wave_number, amp, random_shift
        end if
        do j=1, line_size !reloop over particles
          line_position=j+(i-1)*line_size
          select case(wave_type) !what wave type do we want?
            case('planar')       !set this in run.in
              if (k==1) then !normal undefined for a straight line
                f(line_position)%x(1)=f(line_position)%x(1)+&
                amp*delta*sin(random_shift+wave_number*2.*pi*real(2.*j-1)/(2.*line_size))
              else
                f(line_position)%x(:)=f(line_position)%x(:)+&
                normalf(line_position)*amp*delta*sin(random_shift+wave_number*2.*pi*real(2.*j-1)/(2.*line_size))
              end if
            case('helical')
              if (k==1) then !normal/binormal undefined for a straight line
                f(line_position)%x(1)=f(line_position)%x(1)+&
                amp*delta*sin(random_shift+wave_number*2.*pi*real(2.*j-1)/(2.*line_size))
                f(line_position)%x(2)=f(line_position)%x(2)+&
                amp*delta*cos(random_shift+wave_number*2.*pi*real(2.*j-1)/(2.*line_size))
              else
                f(line_position)%x(:)=f(line_position)%x(:)+&
                normalf(line_position)*amp*delta*sin(random_shift+wave_number*2.*pi*real(2.*j-1)/(2.*line_size))+&
                binormalf(line_position)*amp*delta*cos(random_shift+wave_number*2.*pi*real(2.*j-1)/(2.*line_size))
              end if
            case default
              call fatal_error('initial.mod:wave_line','incorrect wave type parameter')
          end select
        end do
        call ghostp !we must call this routine at the end of adding every wave 
                    !the routines normalf, binormalf rely on correct ghostpoints
      end do !close k loop
      if (i==1) then
        close(34)
      end if 
    end do !closes the i loop
  end subroutine
  !*************************************************************************
  subroutine setup_random_loops
    implicit none
    real :: loop_radius
    real :: anglex,angley,anglez
    real,dimension(3)::translate
    real,dimension(3)::dummy_xp_1, dummy_xp_2, dummy_xp_3, dummy_xp_4
    integer :: loop_size
    integer:: loop_position
    integer :: i,j
    !test run.in parameters, if wrong program will exit
    if (line_count==0) then
      call fatal_error('init.mod:setup_random_loops', &
      'you have not set a value for line_count in run.in')
    end if
    if (mod(pcount,line_count)/=0) then
      call fatal_error('init.mod:setup_random_loops', &
      'pcount/line_count is not an integer')
    end if
    loop_size=int(pcount/line_count)
    loop_radius=loop_size*(0.75*delta)/(2*pi) !75% of potential size
    do i=1, line_count
      call random_number(anglex)
      call random_number(angley)
      call random_number(anglez)
      call random_number(translate)
      anglex=anglex*2*pi
      angley=angley*2*pi
      anglez=anglez*2*pi
      translate=(box_size*translate-box_size/2.)-loop_radius
        
      do j=1, loop_size

        loop_position=j+(i-1)*loop_size

        dummy_xp_1(1)=loop_radius*sin(pi*real(2*j-1)/loop_size)
        dummy_xp_1(2)=loop_radius*cos(pi*real(2*j-1)/loop_size)
        dummy_xp_1(3)=0.

        dummy_xp_2(1)=dummy_xp_1(1)
        dummy_xp_2(2)=dummy_xp_1(2)*cos(anglex)+dummy_xp_1(3)*sin(anglex)
        dummy_xp_2(3)=dummy_xp_1(2)*(-sin(anglex))+dummy_xp_1(3)*cos(anglex)

        dummy_xp_3(1)=dummy_xp_2(1)*cos(angley)+dummy_xp_2(3)*(-sin(angley))
        dummy_xp_3(2)=dummy_xp_2(2)
        dummy_xp_3(3)=dummy_xp_2(1)*sin(angley)+dummy_xp_2(3)*cos(angley)

        dummy_xp_4(1)=dummy_xp_3(1)*cos(anglez)+dummy_xp_3(2)*sin(anglez)
        dummy_xp_4(2)=dummy_xp_3(1)*(-sin(anglez))+dummy_xp_3(2)*cos(anglez)
        dummy_xp_4(3)=dummy_xp_3(3)
    
        f(loop_position)%x(1)=dummy_xp_4(1)+translate(1)
        f(loop_position)%x(2)=dummy_xp_4(2)+translate(2)
        f(loop_position)%x(3)=dummy_xp_4(3)+translate(3)

        if(j==1) then
          f(loop_position)%behind=i*loop_size
          f(loop_position)%infront=loop_position+1
        else if (j==loop_size) then
          f(loop_position)%behind=loop_position-1
          f(loop_position)%infront=(i-1)*loop_size+1
        else
          f(loop_position)%behind=loop_position-1
          f(loop_position)%infront=loop_position+1
        end if
        f(loop_position)%u1=0. ; f(loop_position)%u2=0.
      end do
    end do
  end subroutine
!*************************************************************************
  subroutine setup_tangle
    implicit none
    real :: rand1, rand2, rand3, rand4
    integer :: pcount_required
    integer :: line_size
    integer:: line_position
    integer :: i,j
    !test run.in parameters, if wrong program will exit
    if (line_count==0) then
      call fatal_error('init.mod:setup_tangle', &
      'you have not set a value for line_count in run.in')
    end if
    if (periodic_bc) then
      !work out the number of particles required for our lines
      !given the box size specified in run.i
      pcount_required=line_count*nint(box_size/(0.75*delta)) !75%
      write(*,*) 'changing size of pcount to fit with box_length and delta'
      write(*,*) 'pcount is now', pcount_required
      deallocate(f) ; pcount=pcount_required ; allocate(f(pcount))
    else
      call fatal_error('init.mod:setup_tangle',&
                       'periodic boundary conditions required for this initial condition') !cdata.mod
    end if
    line_size=int(pcount/line_count)
    do i=1, line_count
      call random_number(rand1)
      call random_number(rand2)
      call random_number(rand3)
      call random_number(rand4)
      rand1=(rand1-.5)*box_size
      rand3=(rand3-.5)*box_size
      if (rand4<0.5) then
        rand4=-1.
      else
        rand2=1.
      end if
      do j=1, line_size
        line_position=j+(i-1)*line_size
        if (rand2<0.5) then
          f(line_position)%x(1)=rand4*(box_size/10.)*sin(pi*real(2*j-1)/(2.*line_size))+rand3
          f(line_position)%x(2)=rand1
          f(line_position)%x(3)=-box_size/2.+box_size*real(2*j-1)/(2.*line_size)
        else
          f(line_position)%x(1)=rand1
          f(line_position)%x(2)=rand4*(box_size/10.)*sin(pi*real(2*j-1)/(2.*line_size))+rand3
          f(line_position)%x(3)=-box_size/2.+box_size*real(2*j-1)/(2.*line_size)
        end if
        if(j==1) then
          f(line_position)%behind=i*line_size
          f(line_position)%infront=line_position+1
        else if (j==line_size) then
          f(line_position)%behind=line_position-1
          f(line_position)%infront=(i-1)*line_size+1
        else
          f(line_position)%behind=line_position-1
          f(line_position)%infront=line_position+1
        end if
        f(line_position)%u1=0. ; f(line_position)%u2=0.
      end do
    end do
  end subroutine
!*************************************************************************
  subroutine setup_criss_cross
    implicit none
    real :: rand1, rand2, rand3, rand4
    integer :: pcount_required
    integer :: line_size
    integer:: line_position
    integer :: i,j
    !test run.in parameters, if wrong program will exit
    if (line_count==0) then
      call fatal_error('init.mod:setup_tangle', &
      'you have not set a value for line_count in run.in')
    end if
    if (periodic_bc) then
      !work out the number of particles required for our lines
      !given the box size specified in run.i
      pcount_required=line_count*nint(box_size/(0.75*delta)) !75%
      write(*,*) 'changing size of pcount to fit with box_length and delta'
      write(*,*) 'pcount is now', pcount_required
      deallocate(f) ; pcount=pcount_required ; allocate(f(pcount))
    else
      call fatal_error('init.mod:setup_criss_cross',&
                       'periodic boundary conditions required for this initial condition') !cdata.mod
    end if
    line_size=int(pcount/line_count)
    do i=1, line_count
      call random_number(rand1)
      call random_number(rand2)
      call random_number(rand3)
      call random_number(rand4)
      rand1=(rand1-.5)*box_size
      rand2=(rand2-.5)*box_size
      if (rand4<0.5) then
        rand4=-1
      else
        rand4=1
      end if
      if (rand3<0.33333333) then
        do j=1, line_size
          line_position=j+(i-1)*line_size
          f(line_position)%x(1)=rand1
          f(line_position)%x(2)=rand2
          f(line_position)%x(3)=rand4*(-box_size/2.+box_size*real(2*j-1)/(2.*line_size))
          if(j==1) then
            f(line_position)%behind=i*line_size
            f(line_position)%infront=line_position+1
          else if (j==line_size) then
            f(line_position)%behind=line_position-1
            f(line_position)%infront=(i-1)*line_size+1
          else
            f(line_position)%behind=line_position-1
            f(line_position)%infront=line_position+1
          end if
          f(line_position)%u1=0. ; f(line_position)%u2=0.
        end do
      else if (rand3<0.6666666) then
        do j=1, line_size
          line_position=j+(i-1)*line_size
          f(line_position)%x(1)=rand1
          f(line_position)%x(2)=rand4*(-box_size/2.+box_size*real(2*j-1)/(2.*line_size))
          f(line_position)%x(3)=rand2
          if(j==1) then
            f(line_position)%behind=i*line_size
            f(line_position)%infront=line_position+1
          else if (j==line_size) then
            f(line_position)%behind=line_position-1
            f(line_position)%infront=(i-1)*line_size+1
          else
            f(line_position)%behind=line_position-1
            f(line_position)%infront=line_position+1
          end if
          f(line_position)%u1=0. ; f(line_position)%u2=0.
        end do
      else
        do j=1, line_size
          line_position=j+(i-1)*line_size
          f(line_position)%x(1)=rand4*(-box_size/2.+box_size*real(2*j-1)/(2.*line_size))
          f(line_position)%x(2)=rand1
          f(line_position)%x(3)=rand2
          if(j==1) then
            f(line_position)%behind=i*line_size
            f(line_position)%infront=line_position+1
          else if (j==line_size) then
            f(line_position)%behind=line_position-1
            f(line_position)%infront=(i-1)*line_size+1
          else
            f(line_position)%behind=line_position-1
            f(line_position)%infront=line_position+1
          end if
          f(line_position)%u1=0. ; f(line_position)%u2=0.
        end do
      end if
    end do
  end subroutine
!****************************************************************
  subroutine setup_mesh
    implicit none
    integer :: i,j,k
    real :: x,y,z
    if (mod(mesh_size,2)/=0) then
      call fatal_error('init.mod:setup_mesh', &
      'mesh size must be a multiple of 2')
    end if
    if (mesh_size<16) call warning_message('setup_mesh','warning mesh size is small')
    mesh_delta=real(box_size)/mesh_size
    write(*,'(a,i3.2,a,f7.6)') ' creating an n^3 mesh, n=', mesh_size, ' resolution=', mesh_delta
    write(*,'(a,i5.2,a)') ' mesh information will be printed every ', mesh_shots, ' time-steps'
    allocate(mesh(mesh_size,mesh_size,mesh_size))
    do k=1, mesh_size
      do j=1, mesh_size
        do i=1, mesh_size
          x=mesh_delta*real(2*i-1)/2.-(box_size/2.)
          y=mesh_delta*real(2*j-1)/2.-(box_size/2.)
          z=mesh_delta*real(2*k-1)/2.-(box_size/2.)
          mesh(k,j,i)%x(1)=x ; mesh(k,j,i)%x(2)=y ; mesh(k,j,i)%x(3)=z
          !clear the velocity slots - for safety
          mesh(k,j,i)%u_sup=0. ; mesh(k,j,i)%u_norm=0.
        end do
      end do
    end do
  end subroutine
  !******************************************************************
  subroutine timestep_check
    implicit none
    real :: delta_min, dt_max
    select case(velocity)
      case('Off')
        return
    end select
    delta_min=delta/2.
    dt_max=((delta_min)**2)/(quant_circ*log(delta_min*1E8/pi))
    if (dt<dt_max) then
      write(*,'(a,e10.4)') ' dt is below maximum possible dt:', dt_max
    else
      write(*,'(a,e10.4)') ' warning set dt below ', dt_max
      call fatal_error('initial.mod:timestep_check','dt is too large')
    end if
  end subroutine
end module
!******************************************************************
!>\page INIT Initial condition for filament
!!Initial condition for filament is set in run.in throught the parameter initf\n
!!Options are:\n
!!- \p single_loop - Single loop in the x-y plane (z=0). The size of the loop is dictated by the number of intial particles set.
!!- \p single_line - Single line from -z to +z (x=y=0). The number of particles needed is set by the box size and delta, pcount is automatically adjusted to fit this.
!!- \p random_loops - multiple random loops (number set by the parameter \p line_count) distributed throughout the box. The size of the loops is dictated by the number of loops and the number of intial particles set. You must ensure \p pcount is a multiple of the number of loops
!!- \p crow - Two lines from -z to +z (x=y=0), anti-parallel to trigger the
!!crow instability. The number of particles needed is set by the box size and
!!delta, pcount is automatically adjusted to fit this.
!!- \p leap-frog - two loops in the x-y plane at z=-delta, z=delta this sets
!!the loops leapfrogging (not with LIA). 

