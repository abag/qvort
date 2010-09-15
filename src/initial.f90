module initial
  !INITIAL/RESTART CONDITIONS
  use cdata
  use normal_fluid
  use forcing
  contains
  !*************************************************************************
  subroutine init_setup()
    !prints and sets up intial conditions
    use quasip
    implicit none
    logical :: restart
    write(*,'(a,f9.7)') ' quantum of circulation is:', quant_circ 
    !we must check the timestep is sufficient to resolve the motion
    !based on the smallest separation possible in the code
    call timestep_check !initial.mod
    !how is data being outputted (binary or formatted)
    if (binary_print) then
      write(*,*) 'binary data output, formatted data can be selected in run.in'
    else
      write(*,*) 'formatted data output selected in run.in'
    end if
    !periodic bounday conditions?
    if (box_size>0.) then
      periodic_bc=.true.
      write(*,'(a,f6.3)') ' running with periodic boundaries, box size:', box_size
    else 
      write(*,*) 'running with open boundaries'
    end if
    call setup_normal_fluid !normal_fluid.mod
    !sort out if there is a special dump time
    int_special_dump=int(special_dump/dt) !convert to integer

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
        case('kivotedes')
          call setup_kivotedes !init.mod
        case('wave_spec')
          call setup_wave_spec !init.mod
        case('wave_line')
          call setup_wave_line !init.mod
        case('line_motion')
          call setup_line_motion !init.mod
        case('tangle')
          call setup_tangle !init.mod
        case default
          call fatal_error('cdata.mod:init_setup', &
                         'invalid choice for initf parameter') !cdata.mod
       end select
     end if
     !test if we have a non-zero mesh size
     if (mesh_size>0) then
       if (periodic_bc) then
         call setup_mesh !init.mod
       else 
         call fatal_error('cdata.mod:init_setup', &
         'running with a non-zero mesh size requires periodic BCs') !cdata.mod
       end if
     else
       write(*,*) 'velocity fields not being stored on a mesh - (no spectra etc.)'
     end if
     !do we employ forcing on the boundary?
     call setup_forcing !forcing,mod
     !are there particles in the code?
     if (quasi_pcount>0) then
       call setup_quasip !quasip.mod
     else
       write(*,*) 'No particles in the code'
     end if
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
             call fatal_error('initial.mod',&
                  'no superfluid velocity and no normal velocity')
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
      case default
        print*, 'correct value for velocity in run.in has not been set'
        print*, 'options are: LIA, BS, Tree'
        call fatal_error('init.mod:init_setup', & 
        'correct value for "velocity" in run.in has not been set') !cdata.mod
    end select
    !any special diagnostic information?
    if (curv_hist) write(*,*) 'printing histograms of curvature to file'
  end subroutine
  !**********************************************************************
  subroutine data_restore
    !restart the code
    implicit none
    integer :: dummy_itime 
    open(unit=63,file="./data/var.dat",FORM='unformatted')
      read(63) pcount
      read(63) recon_count
      read(63) dummy_itime
      read(63) t
      allocate(f(pcount))
      read(63) f
    close(63)
    nstart=dummy_itime+1
    write(*,*) 'data read in from dump file at t=', t
    if (quasi_pcount>0) then
      write(*,*)'note that the quasi particles will not be restored'
    end if
  end subroutine
  !*************************************************************************
  subroutine setup_single_loop
    !a loop in the x-y plane
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
  subroutine setup_single_line
    !a line from the lop of the box to the bottom, no curvature
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
      call fatal_error('init.mod:setup_line_motion', &
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
  subroutine setup_leap_frog
    !two loops in the x-y plane
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
  subroutine setup_crow
    !anti parallel lines (instability)
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
      f(i)%x(1)=radius*sin(pi*real(2.*(i-pcount/2)-1)/(pcount/2))+radius/0.52
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
  subroutine setup_kivotedes
    !four unlinked loops as in Kivotedes 2001, PRL
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
  subroutine setup_wave_spec
    !a loop in the x-y plane with wave pertubations
    implicit none
    real :: wave_number, prefactor
    real :: amp, random_shift
    real :: radius
    integer :: i, j
    radius=(0.75*pcount*delta)/(2*pi) !75% of potential size 
    write(*,*) ' initf: single loop, radius of loop:', radius 
    write(*,'(i2.1,a,f9.5)') wave_count,' wave pertubations, with spectral slope:', wave_slope
    !loop over particles setting spatial and 'loop' position
    do i=1, pcount
      f(i)%x(1)=radius*sin(pi*real(2*i-1)/pcount)
      f(i)%x(2)=radius*cos(pi*real(2*i-1)/pcount)
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
    !now add pertubations in Z
    prefactor=wave_amp/(2**wave_slope)
    write(*,*) '------------------WAVE INFORMATION-------------------'
    do j=2, wave_count+1 
      wave_number=j
      amp=prefactor*(wave_number**wave_slope)
      call random_number(random_shift)
      random_shift=random_shift*2*pi
      write(*,'(a,f9.5,a,f9.5,a,f9.5)'), ' wavenumber',wave_number,', amplitude', amp, ', shift', random_shift
      do i=1, pcount
        f(i)%x(3)=f(i)%x(3)+amp*delta*sin(random_shift+wave_number*2.*pi*real(2*i-1)/pcount)!+amp*delta*cos(wave_number*2.*pi*real(2*i-1)/pcount)
      end do
    end do
    write(*,*) '-----------------------------------------------------'

  end subroutine
  !*************************************************************************
  subroutine setup_wave_line
    !a line from the lop of the box to the bottom, helical waves added with
    !specific spectrum
    implicit none
    real :: wave_number, prefactor
    real :: amp, random_shift
    integer :: pcount_required
    integer :: i, j
    character(*), parameter :: wave_type='planar' !helical or planar waves
    if (periodic_bc) then
      !work out the number of particles required for single line
      !given the box size specified in run.i
      pcount_required=nint(box_size/delta) !100% as waves are added
      write(*,*) 'changing size of pcount to fit with box_length and delta'
      write(*,*) 'pcount is now', pcount_required
      deallocate(f) ; pcount=pcount_required ; allocate(f(pcount))
    else
      call fatal_error('init.mod:setup_line_motion', &
      'periodic boundary conditions required')
    end if
    write(*,'(i4.1,a,a,f9.5)') wave_count, wave_type,' wave pertubations, with spectral slope:', wave_slope
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
    !now we add the spiral waves
    prefactor=wave_amp/(2**wave_slope)
    write(*,*) ' wave information recorded in ./data/wave_info.log'
    open(unit=34,file='./data/wave_info.log')
    write(*,*) '%------------------WAVE INFORMATION-------------------'
    do j=1, wave_count
      wave_number=2+0.05*j
      amp=prefactor*(wave_number**wave_slope)
      call random_number(random_shift)
      random_shift=random_shift*2*pi
      write(34,'(a,f9.5,a,f9.5,a,f9.5)'), '%wavenumber',wave_number,', amplitude', amp, ', shift', random_shift
      do i=1, pcount
        select case(wave_type)
          case('planar')
            f(i)%x(1)=f(i)%x(1)+amp*delta*sin(random_shift+wave_number*4.*pi*f(i)%x(3)/box_size)
          case('helical')
            f(i)%x(1)=f(i)%x(1)+amp*delta*sin(random_shift+wave_number*4.*pi*f(i)%x(3)/box_size)
            f(i)%x(2)=f(i)%x(2)+amp*delta*cos(random_shift+wave_number*4.*pi*f(i)%x(3)/box_size)
        end select
      end do
    end do
    close(34)
  end subroutine
  !*************************************************************************
  subroutine setup_line_motion
    !a line from the lop of the box to the bottom,
    !given none zero curvature so it moves
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
      call fatal_error('init.mod:setup_line_motion', &
      'periodic boundary conditions required')
    end if
    do i=1, pcount
      f(i)%x(1)=(box_size/10.)*sin(pi*real(2*i-1)/(2.*pcount))
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
      translate=(2*box_size*translate-box_size)
        
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
!****************************************************************
  subroutine setup_mesh
    implicit none
    integer :: i,j,k
    real :: x,y,z
    if (mod(mesh_size,2)/=0) then
      call fatal_error('init.mod:setup_mesh', &
      'mesh size must be a multiple of 2')
    end if
    if (mesh_size<16) write(*,*) 'warning mesh size is small'
    mesh_delta=real(box_size)/mesh_size
    write(*,'(a,i3.2,a,f7.6)') ' creating an n^3 mesh, n=', mesh_size, ' resolution=', mesh_delta
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
    delta_min=delta/2.
    dt_max=((delta_min)**2)/(quant_circ*log(delta_min*1E8/pi))
    if (dt<dt_max) then
      write(*,'(a,f14.13)') ' dt is below maximum possible dt:', dt_max
    else
      write(*,'(a,f14.13)') ' warning set dt below ', dt_max
      call fatal_error('initial.mod:timestep_check','dt is too large')
    end if
  end subroutine
end module
