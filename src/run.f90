!>The main program, runs the qvort code
program run
  use cdata 
  use initial
  use output
  use periodic
  use timestep
  use quasip
  use particles
  use line
  use diagnostic
  use quasip
  use tree
  use inject
  use mirror
  use smoothing 
  use mag
  use sph_interface
  implicit none
  integer :: i
  logical :: can_stop=.false.,can_reload=.false. 
  call banner_print !output.mod
  call init_random_seed !cdata.mod
  !read in parameters
  call read_run_file !cdata.mod
  !initial conditions - checks for restart?
  call init_setup !initial.mod
  if (seg_fault) write(*,*) 'here1'
  !print dimensions info for matlab plotting
  call print_dims !output.mod
  !begin time loop
  write(*,*) 'setup complete: beginning time loop'
  do itime=nstart, nsteps
    !----------------vortex injection------------------
    if ((mod(itime,inject_skip)==0).and.(t<inject_stop)) then
      call vortex_inject !inject.mod
    end if
    call ghostp !periodic.mod
    if (magnetic) call set_B_length1 !mag.mod
    !---------------------build tree routines----------------------
    if (tree_theta>0) then
      call construct_tree !tree.mod
    end if
    if (seg_fault) write(*,*) 'here2'
    !---------------------create mirror array----------------------
    if (mirror_bc) call mirror_init !mirror.mod
    !---------------------velocity operations----------------------
    call pmotion !timestep.mod
    if (seg_fault) write(*,*) 'here3'
    if (mod(itime,mesh_shots)==0) then
      call mesh_velocity !timestep.mod
      if (sm_size>0) call get_smoothed_field !smoothing.mod
    end if
    if (seg_fault) write(*,*) 'here4'
    !---------------------line operations--------------------------
    if (magnetic) then
      !we must redo ghost particles (i think!)
      call ghostp !periodic.mod
      call set_B_strength !mag.mod
    end if 
    call pinsert !line.mod
    if (magnetic) then
      !magnetic diffusion
      if (B_nu>epsilon(0.)) call B_diffusion !mag.mod
    end if
    if (seg_fault) write(*,*) 'here5'
    if (mod(itime, recon_shots)==0) then
      if (tree_theta>0) then
        !\todo we may need to empty the tree and then redraw it at this point
        call pclose_tree !tree.mod
      else
        call pclose !line.mod
      end if
      if (switch_off_recon.eqv..false.) call precon !line.mod
      call premove !line.mod  \todo switchoff premove in run.in
      select case(velocity) ; case('SPH')
          !if we have SPH particles can we latch onto them?
          call SPH_f_latch !sph_interface.mod
      end select
    end if
    if (seg_fault) write(*,*) 'here6'
    !-------------------boundary conditions------------------------
    select case(boundary)
      case('periodic')
        call enforce_periodic !periodic.mod
      case('openx')
        call enforce_periodic_yz !periodic.mod
      case('open-remove')
        call enforce_open_removal !periodic.mod
      case ('mirror')
        call mirror_pinning !mirror.mod
    end select
    !---------------once all algorithms have been run--------------
    t=t+dt  !increment the time
    !---------------------diagnostic info--------------------------
    if (mod(itime, shots)==0) then
      call velocity_info !diagnostics.mod
      call curv_info !diagnostics.mod
      if (energy_inf) call energy_info !diagnostics.mod
      if (topo_inf) call get_topo_info !diagnostics.mod
      if (mod(itime, mesh_shots)==0) then
        if (boxed_vorticity) call get_boxed_vorticity !diganostics.mod
        if (sep_inf) call get_sep_inf !diganostics.mod
      end if 
    end if
    if (seg_fault) write(*,*) 'here7'
    !--------------now do all data output--------------------------
    if (mod(itime, shots)==0) then
      !store data to a binary file for reload
      call data_dump !output.mod
      if(particles_only.eqv..false.) call print_info !output.mod
      if (mod(itime, mesh_shots)==0) then
        if (magnetic.and.full_B_print) call print_full_B(itime/mesh_shots)
!mag.f90
        !print the mesh to a binary file
        call print_mesh(itime/mesh_shots) !output.mod
        if (delta_adapt_print) call print_delta_adapt(itime/mesh_shots)!output.mod
        !print the smoothed mesh to binary file
        if (sm_size>0) call print_smooth_mesh(itime/mesh_shots)!smoothing.mod
        !can also print full velocity field for statistics 
        if (vel_print) call print_velocity(itime/mesh_shots)!output.mod
        call one_dim_vel(itime/mesh_shots) !diagnostics.mod
        call two_dim_vel(itime/mesh_shots) !diagnostics.mod
      end if
      if (magnetic) call B_ts !mag.mod
    end if
    if (seg_fault) write(*,*) 'here8'
    !---do we have particles in the code - if so evolve them
    if (quasi_pcount>0) call quasip_evolution !quasip.mod
    if (part_count>0) call particles_evolution !particles.mod
    if (SPH_count>0) call SPH_evolution !sph.mod
    if (seg_fault) write(*,*) 'here9'
    !-------------------remove the tree----------------------------
    if (tree_theta>0) then
      call empty_tree(vtree) !empty the tree to avoid a memory leak
      if (seg_fault) write(*,*) 'here9a'
      deallocate(vtree%parray) ; deallocate(vtree)
      if (seg_fault) write(*,*) 'here9b'
      nullify(vtree) !just in case!
    end if
    if (seg_fault) write(*,*) 'here10'
    !---------------------close mirror array-----------------------
    if (mirror_bc) then
      !call mirror_checker !mirror.mod THIS DOESNT EXIST!!!!?
      call mirror_close !mirror.mod
    end if
    if (seg_fault) write(*,*) 'here11'
    !---------------------special data dump------------------------
    if ((itime/=0).and.(itime==int_special_dump)) call sdata_dump
    !---------------------can we exit early------------------------
    if (mod(itime,shots)==0) then
      inquire(file='./STOP', exist=can_stop)
      if (can_stop) then
        write(*,*) 'Encountered stop file, ending simulation'
        stop
      end if
    end if
    !---------------------do we need to reload run.in------------------------
    if (mod(itime,shots)==0) then
      inquire(file='./RELOAD', exist=can_reload)
      if (can_reload) then
        write(*,*) 'Encountered reload file, reloading run.in'
        call reload_run_file !cdata.mod
        can_reload=.false.
        call system('rm ./RELOAD') !remove the RELOAD file
      end if
    end if
    if (seg_fault) write(*,*) 'here12'
    !--------------------final sanity checks----------------------
    if (NAN_test) call NAN_finder !general.mod
  end do
  !deallocate(f,mesh) !tidy up
end program
!------------MAIN PAGE FOR DOXYGEN IS HERE----------------
!>\mainpage Qvort Main Page
!>\author Andrew Baggaley
!>\date 2010-2011
!>\section intro_sec Running the code
!>to compile the code type (note > denotes terminal input)
!>
!>>make
!>
!>to run the code type
!>
!>>./run.sh
!>
!>more options for starting the code can be found by typing
!>
!>>./run.sh -h
!>\section Plotting
!>The main visualisation routines for the code are written in matlab,
!>to visualise first open matlab, type
!>
!>>matlab -nosplash -nojvm
!>
!>note this confines matlab to a terminal as opposed to the full matlab window.
!>To visualise snapshots of the simulation in the matlab window type
!>
!>>vortex_plot(n)
!>
!>where n is the snapshot you wish to visualise (1,10,999,etc.). This is a little slow
!>for a large number of vortex points, type 
!>
!>>help vortex_plot 
!>
!>for other options, this also includes information on making a movie from the output.
!>To visualise time series information type the following in the matlab window
!>
!>>ts
!>
!>More information on matlab visualisation can be found at \ref MATLAB
!>\subsection Gnuplot
!>To do very simply visualisation of the filaments set the following arguement
!>in run.in, binary_print F, this will now output formatted data, if there are
!>then x number of files to process run the following command
!>
!>>./scripts/qvort_animate_gnuplot.sh x
!>
!>this will create x very rough snapshots which can be animated with the
!>following line (if ImageMagick is installed):
!>
!>>animate data/*.png
!>
!>In a similar manner the user can create a simple time series plot using gnuplot
!>by typing the following into a terminal (assuming you have added the qvort scripts
!>into your path, see below)
!>
!>>qvort_ts_gnuplot.sh 
!>
!>One thing that I find very useful is to quickly check the line length of the filament
!>whilst connected to a remote computer over ssh. To save time/bandwith over a poor connection
!>the user can simply run the following script to get a rough plot of the evolution of the line
!>length inside the terminal:
!>
!>>qvort_ll_term.sh
!>
!>this script can easily be copied and modified to inform the user of other variables of interest.
!>\subsection Octave
!>If you do not have matlab but are a user of octave there is a single octave script to plot the filament,
!>the qvort octave path is automatically added to the path, so simply run:
!>
!>octave_plot(filenumber)
!>\section Scripts
!>I recommend linking the scripts located in ./scripts into you bin, e.g. 
!>
!>> cp  ./scripts/qvort_* ~/Bin
!>
!>this will allow you to use the qvort_nr (new run) script which is the optimal way to 
!>launch a new run. Once the script is in your path simply type
!>
!>>qvort_nr dir
!>
!>where dir is the directory you wish to create the newrun in. 
!>\page MATLAB Matlab visualisation
!>Eventually all the matlab routines should be documented here.
!>vortex_plot.m - the main vortex plotting routine, run this with a filenumber as an input, e.g.\n
!>
!>>vortex_plot(1) \n
!>
!>the routine will also take the follwing options in the form vortex_plot(n,'option1','option2') e.g.\n
!>
!>options are:\n
!>
!>   - rough: plots particles only advised for a large number of particles \n
!>   - line: looks better than above a nice compromise\n
!>   - dark: add a night time theme!\n
!>   - rainbow: colour code the vortex according to velocity\n
!>   - magnetic: colour code a magnetic flux tube according to log2(B) -this automatically switches on rainbow\n
!>   - overhead: angle the plot overhead\n
!>   - show_points: will only work if line is set, shows points as well as lines, ignored if rainbow set\n
!>   - print: print to file rather than screen\n
!>   - eps: if print is set and eps is set then output eps files\n
!>   - movie: make a movie by outputting lots of pngs - for batch mode use vortex_anim.m\n
!>
!>ts.m - read in the time series file (data/ts.log) and plot various dignostic information 
!>if given the option print will print to .eps file rather than screen, at present the following information
!>is plotted:
!>
!> - figure 1 - vortex point information
!>   -# t vs particle count
!>   -# t vs reconnection count
!>   -# t vs mean particle separation (between neighbours)
!>   -# t vs total line length
!> - figure 2 - velocity information
!>   -# t vs maximum velocity
!>   -# t vs maximum change of velocity
!> - figure 3 - number of evaluations (per point)
!> - figure 4 - mean curvature
!> - figure 5 - number of particles removed (not in reconnections)
!>
!> if fluid particle are in the code additional information will be plotted



