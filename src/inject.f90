!>the routines contained within this module will inject loops into the code every nth timestep
!>(n is selected in run.in via inject_skip). This is useful to model specific
!>experiments.
!>The size of the loop is set using inject_size and the 
!>topology of the injected loop is set via loop_type, again in run.in.
module inject
  use cdata
  use general
  use periodic
  integer,private :: injection_direction=1
  contains
  !>check all the necessary conditions to inject vortices are set in run.in
  subroutine setup_vortex_injection()
    implicit none
    integer :: inject_size_required !lucy added this
    select case(inject_type)
      case('off')
        return !leave the routine
    end select
    !check inject size is not too small
    if (inject_size<5) call fatal_error('inject.mod','inject size is too small')
    !print to screen pertinent details
    write(*,'(a)') ' ------------------VORTEX INJECTION-------------------'
    !check requirements of specific conditions
    select case(inject_type)
      case('edge_pulse')
        if (mod(inject_size,6)/=0) then
          call fatal_error('inject.mod','inject size must be a multiple of 6')
        end if
      case('loop_stream')
        if (mod(inject_size,2)/=0) then
          call fatal_error('inject.mod','inject size must be a multiple of 2')
        end if
        write(*,*) 'I recommend injecting', nint(0.245*box_size*4*pi/(0.75*delta)), ' particles'
        write(*,*) 'switching injection direction every ', inject_freq, ' timesteps'
      !******************lucy new code******************************************
      case('random_bundles')
        !work out number of particles required and change inject_size to suit
        inject_size_required=bundle_line_count*nint(box_size/(0.75*delta))
        write(*,*) 'changing size of inject_size to fit with box_length/delta and bundle_line_count'
        write(*,*) 'inject_size is now', inject_size_required
        inject_size=inject_size_required
        select case(bundle_distribution)
          case('uniform')
            write(*,*) 'lines will be uniformly distributed within bundle'
          case('normal')
            write(*,*) 'lines will be normaly distributed within bundle'
          case default
            call fatal_error('inject.mod','bundle_distribution must be normal or uniform')
        end select
        if (bundle_width<0 .or. bundle_width>1) then
          call fatal_error('inject.mod','bundle width must be between 0 and 1')
        end if
        if (bundle_line_count<4) then
          call fatal_error('inject.mod','bundle line count must be at least 4')
        end if
      !*************************************************************************
    end select 
    write(*,'(a,i4.1,a,i3.1,a)') ' loops will be injected every ', inject_skip, ' timesteps with ', inject_size, ' points'
    write(*,'(a,a)') ' inject type is set to: ', trim(inject_type)
    select case(inject_type)
      case('rand-xyz-loop')
        write(*,'(a,f7.4,a)') ' rotation applied to loops: ', rotation_factor, '*2\pi'
      case('rand-yz-loop-rot')
        write(*,'(a,f7.4,a)') ' rotation applied to loops: ', rotation_factor, '*2\pi'
        write(*,'(a,f7.4,a)') ' loops translated by: ', lattice_ratio, 'box_size'
      case('vibrating_mesh')
        write(*,'(a,i4.1,a)') ' using normal distn., mean loop size is ', inject_size, ' points'
        write(*,'(a,f9.6)') ' standard deviation is ', line_sigma
        write(*,*) 'switching injection direction every ', inject_freq, ' timesteps'
        write(*,*) 'mesh size is ', lattice_ratio, ' of box' 
        if (randomise_injection) write(*,*) ' randomising loop injection'
    end select
    !check if we have set an injection stop time in run.in
    if (inject_stop<1E6) then
      write(*,'(a,f10.4)') ' inject will stop after t= ', inject_stop
    end if
  end subroutine
  !*******************************************************************
  !>the routine which injects the vortices into the code
  subroutine vortex_inject
    implicit none    
    type(qvort), allocatable, dimension(:) :: tmp !used to copy array to
    integer, allocatable, dimension(:) :: empty_points_array !location of 0s in f(:)%infront
    integer, allocatable, dimension(:) :: inject_loc !the location where we insert new points
    real :: radius !used for loops
    real :: rand1, rand2, rand3, rand4, theta, R !random numbers lucy changes here
    real :: anglex,angley,anglez !for rotating
    real,dimension(3)::dummy_xp_1, dummy_xp_2, dummy_xp_3, dummy_xp_4
    integer :: old_pcount, dummy_inject_size, bundle_line_size, bundle_line_position !lucy again
    integer :: i, j, counter !to loop
    !check that inject_size is not 0
    if (inject_size==0) return !leave routine if it is
    !check if inject is set to off
    select case(inject_type)
      case('off')
        return !leave routine
      case('vibrating_mesh')
        dummy_inject_size=nint(rnorm(real(inject_size),line_sigma**2))
      case default
        dummy_inject_size=inject_size
    end select
    if (dummy_inject_size<3) then
      call warning_message('inject','skipping injection, very small or negative loop size, reduce line_sigma')
      return
    end if
    allocate(inject_loc(dummy_inject_size)) !we always need this array
    !--------------see if we can insert points into existing gaps----------
    if (count(mask=f(:)%infront==0)>=dummy_inject_size) then
      !print*, 'I have found some space'
      !print*, 'space is ', count(mask=f(:)%infront==0), ' long'
      !we have room find the location of the 0's
      allocate(empty_points_array(pcount))
      where (f(:)%infront==0)
        empty_points_array=1
      elsewhere
        empty_points_array=0
      end where
      counter=1 !initialise the counter
      do i=1, pcount
        if (empty_points_array(i)==1) then
          inject_loc(counter)=i
          !print*, 'location ', counter , ' in the inject loc array: ', inject_loc(counter)
          counter=counter+1
          !once we have found dummy_inject_size slots leave the loop
          if (counter==dummy_inject_size+1) exit
        end if
      end do
      deallocate(empty_points_array) !deallocate this dummy array here 
    else
      !expand the f(:) array
      old_pcount=pcount !first store the old point count
      !now we must resize the array 
      !increment by inject_size using a dummy array tmp
      !print*, 'I am expanding the array'
      allocate(tmp(size(f)+dummy_inject_size)) ; tmp(:)%infront=0 !0 the infront array
      !copy accross information
      tmp(1:size(f)) = f
      !now deallocate tmp and transfer particles back to f
      !move_alloc is new intrinsic function (fortran 2003)
      call move_alloc(from=tmp,to=f)
      pcount=pcount+dummy_inject_size !increase the particle count
      do i=1, dummy_inject_size !now set the location of these points in inject_loc
        inject_loc(i)=old_pcount+i 
      end do
    end if
    !---------------now check what the inject_type is--------------
    select case(inject_type)
      case('edge_pulse')!six rings at each face of the computational box
        radius=(0.75*dummy_inject_size*delta)/(6*2*pi) !75% of potential size 
        !loop over particles (6 loops) setting spatial and 'loop' position
        !first loop - xy plane top of box
        do i=1, dummy_inject_size/6
          f(inject_loc(i))%x(1)=radius*sin(pi*real(2*i-1)/(dummy_inject_size/6))
          f(inject_loc(i))%x(2)=radius*cos(pi*real(2*i-1)/(dummy_inject_size/6))
          f(inject_loc(i))%x(3)=box_size/2.01
          if (i==1) then
            f(inject_loc(i))%behind=inject_loc(dummy_inject_size/6)
            f(inject_loc(i))%infront=inject_loc(i+1)
          else if (i==dummy_inject_size/6) then 
            f(inject_loc(i))%behind=inject_loc(i-1)
            f(inject_loc(i))%infront=inject_loc(1)
          else
            f(inject_loc(i))%behind=inject_loc(i-1)
            f(inject_loc(i))%infront=inject_loc(i+1)
          end if
          !zero the stored velocities
          f(inject_loc(i))%u1=0. ; f(inject_loc(i))%u2=0. ; f(inject_loc(i))%u3=0.
        end do
        !second loop- xy plane bottom of box
        do i=dummy_inject_size/6+1, dummy_inject_size/3
          f(inject_loc(i))%x(1)=-radius*sin(pi*real(2*i-1)/(dummy_inject_size/6))
          f(inject_loc(i))%x(2)=radius*cos(pi*real(2*i-1)/(dummy_inject_size/6))
          f(inject_loc(i))%x(3)=-box_size/2.01
          if (i==(dummy_inject_size/6+1)) then
            f(inject_loc(i))%behind=inject_loc(dummy_inject_size/3)
            f(inject_loc(i))%infront=inject_loc(i+1)
          else if (i==dummy_inject_size/3) then 
            f(inject_loc(i))%behind=inject_loc(i-1)
            f(inject_loc(i))%infront=inject_loc(dummy_inject_size/6+1)
          else
            f(inject_loc(i))%behind=inject_loc(i-1)
            f(inject_loc(i))%infront=inject_loc(i+1)
          end if
          !zero the stored velocities
          f(inject_loc(i))%u1=0. ; f(inject_loc(i))%u2=0. ; f(inject_loc(i))%u3=0.
        end do    
        !third loop- yz plane pos side of box 
        do i=dummy_inject_size/3+1, dummy_inject_size/2
          f(inject_loc(i))%x(1)=box_size/2.01
          f(inject_loc(i))%x(2)=radius*cos(pi*real(2*i-1)/(dummy_inject_size/6))
          f(inject_loc(i))%x(3)=-radius*sin(pi*real(2*i-1)/(dummy_inject_size/6))
          if (i==(dummy_inject_size/3+1)) then
            f(inject_loc(i))%behind=inject_loc(dummy_inject_size/2)
            f(inject_loc(i))%infront=inject_loc(i+1)
          else if (i==dummy_inject_size/2) then 
            f(inject_loc(i))%behind=inject_loc(i-1)
            f(inject_loc(i))%infront=inject_loc(dummy_inject_size/3+1)
          else
            f(inject_loc(i))%behind=inject_loc(i-1)
            f(inject_loc(i))%infront=inject_loc(i+1)
          end if
          !zero the stored velocities
          f(inject_loc(i))%u1=0. ; f(inject_loc(i))%u2=0. ; f(inject_loc(i))%u3=0.
        end do
        !fourth loop- yz plane neg side of box 
        do i=dummy_inject_size/2+1, 2*dummy_inject_size/3
          f(inject_loc(i))%x(1)=-box_size/2.01
          f(inject_loc(i))%x(2)=radius*cos(pi*real(2*i-1)/(dummy_inject_size/6))
          f(inject_loc(i))%x(3)=radius*sin(pi*real(2*i-1)/(dummy_inject_size/6))
          if (i==(dummy_inject_size/2+1)) then
            f(inject_loc(i))%behind=inject_loc(2*dummy_inject_size/3)
            f(inject_loc(i))%infront=inject_loc(i+1)
          else if (i==2*dummy_inject_size/3) then 
            f(inject_loc(i))%behind=inject_loc(i-1)
            f(inject_loc(i))%infront=inject_loc(inject_size/2+1)
          else
            f(inject_loc(i))%behind=inject_loc(i-1)
            f(inject_loc(i))%infront=inject_loc(i+1)
          end if
          !zero the stored velocities
          f(inject_loc(i))%u1=0. ; f(inject_loc(i))%u2=0. ; f(inject_loc(i))%u3=0.
        end do
        !fifth loop- xz plane pos side box
        do i=2*dummy_inject_size/3+1, 5*dummy_inject_size/6
          f(inject_loc(i))%x(1)=-radius*sin(pi*real(2*i-1)/(dummy_inject_size/6))
          f(inject_loc(i))%x(2)=box_size/2.01
          f(inject_loc(i))%x(3)=radius*cos(pi*real(2*i-1)/(dummy_inject_size/6))
          if (i==(2*dummy_inject_size/3+1)) then
            f(inject_loc(i))%behind=inject_loc(5*dummy_inject_size/6)
            f(inject_loc(i))%infront=inject_loc(i+1)
          else if (i==5*dummy_inject_size/6) then 
            f(inject_loc(i))%behind=inject_loc(i-1)
            f(inject_loc(i))%infront=inject_loc(2*inject_size/3+1)
          else
            f(inject_loc(i))%behind=inject_loc(i-1)
            f(inject_loc(i))%infront=inject_loc(i+1)
          end if
          !zero the stored velocities
          f(inject_loc(i))%u1=0. ; f(inject_loc(i))%u2=0. ; f(inject_loc(i))%u3=0.
        end do
        !final loop- xz plane neg side of box
        do i=5*dummy_inject_size/6+1, dummy_inject_size
          f(inject_loc(i))%x(1)=radius*sin(pi*real(2*i-1)/(dummy_inject_size/6))
          f(inject_loc(i))%x(2)=-box_size/2.01
          f(inject_loc(i))%x(3)=radius*cos(pi*real(2*i-1)/(dummy_inject_size/6))
          if (i==(5*dummy_inject_size/6+1)) then
            f(inject_loc(i))%behind=inject_loc(dummy_inject_size)
            f(inject_loc(i))%infront=inject_loc(i+1)
          else if (i==dummy_inject_size) then 
            f(inject_loc(i))%behind=inject_loc(i-1)
            f(inject_loc(i))%infront=inject_loc(5*dummy_inject_size/6+1)
          else
            f(inject_loc(i))%behind=inject_loc(i-1)
            f(inject_loc(i))%infront=inject_loc(i+1)
          end if
          !zero the stored velocities
          f(inject_loc(i))%u1=0. ; f(inject_loc(i))%u2=0. ; f(inject_loc(i))%u3=0.
        end do
      !*********************lucy new code*****************************************
      case('random_bundles')
        bundle_line_size=int(dummy_inject_size/bundle_line_count)
        call random_number(rand1) !determine direction
        call random_number(rand3) !line location in 1st plane
        call random_number(rand4) !line location in 2nd plane
        rand3=(rand3*box_size)-box_size/2.
        rand4=(rand4*box_size)-box_size/2.
        do i=1, bundle_line_count
          if (bundle_distribution=='normal') then
           !andrew help - what random number generator to use
           call random_number(R) !for now until I find normal generator
           R=R*bundle_width*box_size
          else if (bundle_distribution=='uniform') then
           call random_number(R) !line offset from centre of bundle
           R=R*bundle_width*box_size
          end if
          call random_number(theta) !line angle
          theta=theta*2*pi
          call random_number(rand2) !determine orientation i.e. + to - or - to +
          do j=1, bundle_line_size
            bundle_line_position=j+(i-1)*bundle_line_size
           !------------- lines in z direction-----------------
            if (rand1<=0.33333333333) then
             !-----from -z to z--------------
             if (rand2<=0.5) then
               f(inject_loc(bundle_line_position))%x(1)=rand3+R*cos(theta)
               f(inject_loc(bundle_line_position))%x(2)=rand4+R*sin(theta)
               f(inject_loc(bundle_line_position))%x(3)=-box_size/2.+box_size*real(2*j-1)/(2.*bundle_line_size)             
             !---- from z to -z-------------
             else
               f(inject_loc(bundle_line_position))%x(1)=rand3+R*cos(theta)
               f(inject_loc(bundle_line_position))%x(2)=rand4+R*sin(theta)
               f(inject_loc(bundle_line_position))%x(3)=box_size/2.-box_size*real(2*j-1)/(2.*bundle_line_size)
             end if
           !-------------lines in y direction-----------------
            else if (rand1<=0.6666666666) then
             !-----from -y to y-------------
             if (rand2<=0.5) then
               f(inject_loc(bundle_line_position))%x(1)=rand3+R*cos(theta)
               f(inject_loc(bundle_line_position))%x(2)=-box_size/2.+box_size*real(2*j-1)/(2.*bundle_line_size)
               f(inject_loc(bundle_line_position))%x(3)=rand4+R*sin(theta) 
             !-----from y to -y-------------
             else
               f(inject_loc(bundle_line_position))%x(1)=rand3+R*cos(theta)
               f(inject_loc(bundle_line_position))%x(2)=box_size/2.-box_size*real(2*j-1)/(2.*bundle_line_size)
               f(inject_loc(bundle_line_position))%x(3)=rand4+R*sin(theta)
             end if
           !------------lines in x direction------------------
           else
             !-----from -x to x----------------
             if (rand2<=0.5) then 
               f(inject_loc(bundle_line_position))%x(1)=-box_size/2.+box_size*real(2*j-1)/(2.*bundle_line_size)
               f(inject_loc(bundle_line_position))%x(2)=rand3+R*cos(theta)
               f(inject_loc(bundle_line_position))%x(3)=rand4+R*sin(theta)
             !-----from x to -x---------------
             else
               f(inject_loc(bundle_line_position))%x(1)=box_size/2.-box_size*real(2*j-1)/(2.*bundle_line_size)
               f(inject_loc(bundle_line_position))%x(2)=rand3+R*cos(theta)
               f(inject_loc(bundle_line_position))%x(3)=rand4+R*sin(theta)
             end if
           end if
           if(j==1) then
             f(inject_loc(bundle_line_position))%behind=inject_loc(i*bundle_line_size)
             f(inject_loc(bundle_line_position))%infront=inject_loc(bundle_line_position+1)
           else if (j==bundle_line_size) then
             f(inject_loc(bundle_line_position))%behind=inject_loc(bundle_line_position-1)
             f(inject_loc(bundle_line_position))%infront=inject_loc((i-1)*bundle_line_size+1)
           else
             f(inject_loc(bundle_line_position))%behind=inject_loc(bundle_line_position-1)
             f(inject_loc(bundle_line_position))%infront=inject_loc(bundle_line_position+1)
           end if
           f(inject_loc(bundle_line_position))%u1=0.
           f(inject_loc(bundle_line_position))%u2=0.
           f(inject_loc(bundle_line_position))%u3=0.
          end do
        end do
      !***************************************************************************
      case('loop_stream')!two rings top and bottom in corners
        radius=(0.75*dummy_inject_size*delta)/(4*pi) !75% of potential size 
        !loop over particles (2 loops) setting spatial and 'loop' position
        !first loop - xy plane top of box
        do i=1, dummy_inject_size/2
          if (injection_direction==1) then
            f(inject_loc(i))%x(1)=radius*sin(pi*real(2*i-1)/(dummy_inject_size/2))+box_size/5
            f(inject_loc(i))%x(2)=radius*cos(pi*real(2*i-1)/(dummy_inject_size/2))+box_size/5
            f(inject_loc(i))%x(3)=box_size/2.01
          else if (injection_direction==2) then
            f(inject_loc(i))%x(2)=radius*sin(pi*real(2*i-1)/(dummy_inject_size/2))+box_size/5
            f(inject_loc(i))%x(3)=radius*cos(pi*real(2*i-1)/(dummy_inject_size/2))+box_size/5
            f(inject_loc(i))%x(1)=box_size/2.01
          else if (injection_direction==3) then
            f(inject_loc(i))%x(3)=radius*sin(pi*real(2*i-1)/(dummy_inject_size/2))+box_size/5
            f(inject_loc(i))%x(1)=radius*cos(pi*real(2*i-1)/(dummy_inject_size/2))+box_size/5
            f(inject_loc(i))%x(2)=box_size/2.01
          end if
          if (i==1) then
            f(inject_loc(i))%behind=inject_loc(dummy_inject_size/2)
            f(inject_loc(i))%infront=inject_loc(i+1)
          else if (i==dummy_inject_size/2) then 
            f(inject_loc(i))%behind=inject_loc(i-1)
            f(inject_loc(i))%infront=inject_loc(1)
          else
            f(inject_loc(i))%behind=inject_loc(i-1)
            f(inject_loc(i))%infront=inject_loc(i+1)
          end if
          !zero the stored velocities
          f(inject_loc(i))%u1=0. ; f(inject_loc(i))%u2=0. ; f(inject_loc(i))%u3=0.
        end do
        !second loop- xy plane bottom of box
        do i=dummy_inject_size/2+1, dummy_inject_size
          if (injection_direction==1) then
            f(inject_loc(i))%x(1)=-radius*sin(pi*real(2*i-1)/(dummy_inject_size/2))-box_size/5
            f(inject_loc(i))%x(2)=radius*cos(pi*real(2*i-1)/(dummy_inject_size/2))-box_size/5
            f(inject_loc(i))%x(3)=-box_size/2.01
          else if (injection_direction==2) then
            f(inject_loc(i))%x(2)=-radius*sin(pi*real(2*i-1)/(dummy_inject_size/2))-box_size/5
            f(inject_loc(i))%x(3)=radius*cos(pi*real(2*i-1)/(dummy_inject_size/2))-box_size/5
            f(inject_loc(i))%x(1)=-box_size/2.01
          else if (injection_direction==3) then
            f(inject_loc(i))%x(3)=-radius*sin(pi*real(2*i-1)/(dummy_inject_size/2))-box_size/5
            f(inject_loc(i))%x(1)=radius*cos(pi*real(2*i-1)/(dummy_inject_size/2))-box_size/5
            f(inject_loc(i))%x(2)=-box_size/2.01
          end if
          if (i==(dummy_inject_size/2+1)) then
            f(inject_loc(i))%behind=inject_loc(dummy_inject_size)
            f(inject_loc(i))%infront=inject_loc(i+1)
          else if (i==dummy_inject_size) then 
            f(inject_loc(i))%behind=inject_loc(i-1)
            f(inject_loc(i))%infront=inject_loc(dummy_inject_size/2+1)
          else
            f(inject_loc(i))%behind=inject_loc(i-1)
            f(inject_loc(i))%infront=inject_loc(i+1)
          end if
          !zero the stored velocities
          f(inject_loc(i))%u1=0. ; f(inject_loc(i))%u2=0. ; f(inject_loc(i))%u3=0.
        end do    
      case('xy-loop')!loops in xy-plane injected at bottom of box
        radius=(0.75*dummy_inject_size*delta)/(2*pi) !75% of potential size
        !loop over particles setting spatial and 'loop' position
        do i=1, dummy_inject_size
          f(inject_loc(i))%x(1)=radius*sin(pi*real(2*i-1)/dummy_inject_size)
          f(inject_loc(i))%x(2)=radius*cos(pi*real(2*i-1)/dummy_inject_size)
          f(inject_loc(i))%x(3)=-box_size/2.1
          if (i==1) then
            f(inject_loc(i))%behind=inject_loc(dummy_inject_size)
            f(inject_loc(i))%infront=inject_loc(i+1)
          else if (i==dummy_inject_size) then 
            f(inject_loc(i))%behind=inject_loc(i-1)
            f(inject_loc(i))%infront=inject_loc(1)
          else
            f(inject_loc(i))%behind=inject_loc(i-1)
            f(inject_loc(i))%infront=inject_loc(i+1)
          end if
          !zero the stored velocities
          f(inject_loc(i))%u1=0. ; f(inject_loc(i))%u2=0.
        end do
      case('yz-loop') !loops in yz plane injected at side of box
        radius=(0.75*dummy_inject_size*delta)/(2*pi) !75% of potential size
        !loop over particles setting spatial and 'loop' position
        do i=1, dummy_inject_size
          f(inject_loc(i))%x(1)=-box_size/2.1
          f(inject_loc(i))%x(2)=radius*cos(pi*real(2*i-1)/dummy_inject_size)
          f(inject_loc(i))%x(3)=radius*sin(pi*real(2*i-1)/dummy_inject_size)
          if (i==1) then
            f(inject_loc(i))%behind=inject_loc(dummy_inject_size)
            f(inject_loc(i))%infront=inject_loc(i+1)
          else if (i==dummy_inject_size) then 
            f(inject_loc(i))%behind=inject_loc(i-1)
            f(inject_loc(i))%infront=inject_loc(1)
          else
            f(inject_loc(i))%behind=inject_loc(i-1)
            f(inject_loc(i))%infront=inject_loc(i+1)
          end if
          !zero the stored velocities
          f(inject_loc(i))%u1=0. ; f(inject_loc(i))%u2=0.
        end do
      case('vibrating_mesh') !loops in yz plane injected at side of box
        radius=(0.75*dummy_inject_size*delta)/(2*pi) !75% of potential size
        rand1=runif(-box_size/2,box_size/2)*lattice_ratio
        rand2=runif(-box_size/2,box_size/2)*lattice_ratio
        if (randomise_injection) then
          call random_number(rand3) !if randomising generate a random
          if (rand3>0.5) then       !number either -1 or 1
            rand3=1.
          else
            rand3=-1.
          end if
        else
          rand3=1. !else if off always set to be 1
        end if
        !loop over particles setting spatial and 'loop' position
        do i=1, dummy_inject_size
          if (injection_direction==1) then
            f(inject_loc(i))%x(1)=rand3*(-box_size/2.1)
            f(inject_loc(i))%x(2)=rand3*radius*cos(pi*real(2*i-1)/dummy_inject_size)+rand1
            f(inject_loc(i))%x(3)=radius*sin(pi*real(2*i-1)/dummy_inject_size)+rand2
          else if (injection_direction==2) then
            f(inject_loc(i))%x(2)=rand3*(-box_size/2.1)
            f(inject_loc(i))%x(3)=rand3*radius*cos(pi*real(2*i-1)/dummy_inject_size)+rand1
            f(inject_loc(i))%x(1)=radius*sin(pi*real(2*i-1)/dummy_inject_size)+rand2
          else if (injection_direction==3) then
            f(inject_loc(i))%x(3)=rand3*(-box_size/2.1)
            f(inject_loc(i))%x(1)=rand3*radius*cos(pi*real(2*i-1)/dummy_inject_size)+rand1
            f(inject_loc(i))%x(2)=radius*sin(pi*real(2*i-1)/dummy_inject_size)+rand2
          end if
          if (i==1) then
            f(inject_loc(i))%behind=inject_loc(dummy_inject_size)
            f(inject_loc(i))%infront=inject_loc(i+1)
          else if (i==dummy_inject_size) then 
            f(inject_loc(i))%behind=inject_loc(i-1)
            f(inject_loc(i))%infront=inject_loc(1)
          else
            f(inject_loc(i))%behind=inject_loc(i-1)
            f(inject_loc(i))%infront=inject_loc(i+1)
          end if
          !zero the stored velocities
          f(inject_loc(i))%u1=0. ; f(inject_loc(i))%u2=0.
        end do
      case('rand-yz-loop') !as yz-loop but with a slight pertubation in positions
        radius=(0.75*dummy_inject_size*delta)/(2*pi) !75% of potential size
        !loop over particles setting spatial and 'loop' position
        call random_number(rand1)
        call random_number(rand2)
        rand1=box_size*(rand1*2.-1.)/20.
        rand2=box_size*(rand2*2.-1.)/20.
        do i=1, dummy_inject_size
          f(inject_loc(i))%x(1)=-box_size/2.
          f(inject_loc(i))%x(2)=radius*cos(pi*real(2*i-1)/dummy_inject_size)+rand1
          f(inject_loc(i))%x(3)=radius*sin(pi*real(2*i-1)/dummy_inject_size)+rand2
          if (i==1) then
            f(inject_loc(i))%behind=inject_loc(dummy_inject_size)
            f(inject_loc(i))%infront=inject_loc(i+1)
          else if (i==dummy_inject_size) then 
            f(inject_loc(i))%behind=inject_loc(i-1) 
            f(inject_loc(i))%infront=inject_loc(1)
          else
            f(inject_loc(i))%behind=inject_loc(i-1) 
            f(inject_loc(i))%infront=inject_loc(i+1)
          end if
          !zero the stored velocities
          f(inject_loc(i))%u1=0. ; f(inject_loc(i))%u2=0.
        end do
      case('rand-yz-loop-rot') !pertubation in both position and angle
        radius=(0.75*dummy_inject_size*delta)/(2*pi) !75% of potential size
        !loop over particles setting spatial and 'loop' position
        call random_number(rand1)
        call random_number(rand2)
        call random_number(anglex)
        call random_number(angley)
        call random_number(anglez)
        anglex=(2.*anglex-1.)*2*pi*rotation_factor
        angley=(2.*angley-1.)*2*pi*rotation_factor
        anglez=(2.*anglez-1.)*2*pi*rotation_factor
        rand1=box_size*(rand1*2.-1.)*lattice_ratio
        rand2=box_size*(rand2*2.-1.)*lattice_ratio
        do i=1, dummy_inject_size
          dummy_xp_1(1)=0.
          dummy_xp_1(2)=radius*cos(pi*real(2*i-1)/dummy_inject_size)
          dummy_xp_1(3)=radius*sin(pi*real(2*i-1)/dummy_inject_size)

          dummy_xp_2(1)=dummy_xp_1(1)
          dummy_xp_2(2)=dummy_xp_1(2)*cos(anglex)+dummy_xp_1(3)*sin(anglex)
          dummy_xp_2(3)=dummy_xp_1(2)*(-sin(anglex))+dummy_xp_1(3)*cos(anglex)

          dummy_xp_3(1)=dummy_xp_2(1)*cos(angley)+dummy_xp_2(3)*(-sin(angley))
          dummy_xp_3(2)=dummy_xp_2(2)
          dummy_xp_3(3)=dummy_xp_2(1)*sin(angley)+dummy_xp_2(3)*cos(angley)

          dummy_xp_4(1)=dummy_xp_3(1)*cos(anglez)+dummy_xp_3(2)*sin(anglez)
          dummy_xp_4(2)=dummy_xp_3(1)*(-sin(anglez))+dummy_xp_3(2)*cos(anglez)
          dummy_xp_4(3)=dummy_xp_3(3)
    
          f(inject_loc(i))%x(1)=dummy_xp_4(1)-box_size/2.2
          f(inject_loc(i))%x(2)=dummy_xp_4(2)+rand1
          f(inject_loc(i))%x(3)=dummy_xp_4(3)+rand2
          if (i==1) then
            f(inject_loc(i))%behind=inject_loc(dummy_inject_size)
            f(inject_loc(i))%infront=inject_loc(i+1)
          else if (i==dummy_inject_size) then 
            f(inject_loc(i))%behind=inject_loc(i-1)
            f(inject_loc(i))%infront=inject_loc(1)
          else
            f(inject_loc(i))%behind=inject_loc(i-1)
            f(inject_loc(i))%infront=inject_loc(i+1)
          end if
          !zero the stored velocities
          f(inject_loc(i))%u1=0. ; f(inject_loc(i))%u2=0.
        end do
      case('focus-yz-loop')
        radius=(0.75*dummy_inject_size*delta)/(2*pi) !75% of potential size
        !loop over particles setting spatial and 'loop' position
        call random_number(rand1)
        call random_number(rand2)
        call random_number(anglex)
        call random_number(angley)
        call random_number(anglez)
        anglex=(2.*anglex-1.)*2*pi/10.
        angley=(2.*angley-1.)*2*pi/10.
        anglez=(2.*anglez-1.)*2*pi/10.
        rand1=box_size*(rand1*2.-1.)/20.
        rand2=box_size*(rand2*2.-1.)/20.
        do i=1, dummy_inject_size
          dummy_xp_1(1)=-box_size/2.
          dummy_xp_1(2)=radius*cos(pi*real(2*i-1)/dummy_inject_size)+rand1
          dummy_xp_1(3)=radius*sin(pi*real(2*i-1)/dummy_inject_size)+rand2

          dummy_xp_2(1)=dummy_xp_1(1)
          dummy_xp_2(2)=dummy_xp_1(2)*cos(anglex)+dummy_xp_1(3)*sin(anglex)
          dummy_xp_2(3)=dummy_xp_1(2)*(-sin(anglex))+dummy_xp_1(3)*cos(anglex)

          dummy_xp_3(1)=dummy_xp_2(1)*cos(angley)+dummy_xp_2(3)*(-sin(angley))
          dummy_xp_3(2)=dummy_xp_2(2)
          dummy_xp_3(3)=dummy_xp_2(1)*sin(angley)+dummy_xp_2(3)*cos(angley)

          dummy_xp_4(1)=dummy_xp_3(1)*cos(anglez)+dummy_xp_3(2)*sin(anglez)
          dummy_xp_4(2)=dummy_xp_3(1)*(-sin(anglez))+dummy_xp_3(2)*cos(anglez)
          dummy_xp_4(3)=dummy_xp_3(3)
    
          f(inject_loc(i))%x(1)=dummy_xp_4(1)
          f(inject_loc(i))%x(2)=dummy_xp_4(2)
          f(inject_loc(i))%x(3)=dummy_xp_4(3)
          if (i==1) then
            f(inject_loc(i))%behind=inject_loc(dummy_inject_size)
            f(inject_loc(i))%infront=inject_loc(i+1)
          else if (i==dummy_inject_size) then 
            f(inject_loc(i))%behind=inject_loc(i-1)
            f(inject_loc(i))%infront=inject_loc(1)
          else
            f(inject_loc(i))%behind=inject_loc(i-1)
            f(inject_loc(i))%infront=inject_loc(i+1)
          end if
          !zero the stored velocities
          f(inject_loc(i))%u1=0. ; f(inject_loc(i))%u2=0.
        end do
      case('rand-xyz-loop') !completetly random loops
        radius=(0.75*dummy_inject_size*delta)/(2*pi) !75% of potential size
        !loop over particles setting spatial and 'loop' position
        call random_number(rand1)
        call random_number(rand2)
        call random_number(rand3)
        call random_number(anglex)
        call random_number(angley)
        call random_number(anglez)
        anglex=(2.*anglex-1.)*2*pi*rotation_factor
        angley=(2.*angley-1.)*2*pi*rotation_factor
        anglez=(2.*anglez-1.)*2*pi*rotation_factor
        rand1=box_size*(rand1*2.-1.)/2.
        rand2=box_size*(rand2*2.-1.)/2.
        rand3=box_size*(rand3*2.-1.)/2.
        do i=1, dummy_inject_size
          dummy_xp_1(1)=0.
          dummy_xp_1(2)=radius*cos(pi*real(2*i-1)/dummy_inject_size)
          dummy_xp_1(3)=radius*sin(pi*real(2*i-1)/dummy_inject_size)

          dummy_xp_2(1)=dummy_xp_1(1)
          dummy_xp_2(2)=dummy_xp_1(2)*cos(anglex)+dummy_xp_1(3)*sin(anglex)
          dummy_xp_2(3)=dummy_xp_1(2)*(-sin(anglex))+dummy_xp_1(3)*cos(anglex)

          dummy_xp_3(1)=dummy_xp_2(1)*cos(angley)+dummy_xp_2(3)*(-sin(angley))
          dummy_xp_3(2)=dummy_xp_2(2)
          dummy_xp_3(3)=dummy_xp_2(1)*sin(angley)+dummy_xp_2(3)*cos(angley)

          dummy_xp_4(1)=dummy_xp_3(1)*cos(anglez)+dummy_xp_3(2)*sin(anglez)
          dummy_xp_4(2)=dummy_xp_3(1)*(-sin(anglez))+dummy_xp_3(2)*cos(anglez)
          dummy_xp_4(3)=dummy_xp_3(3)
    
          f(inject_loc(i))%x(1)=dummy_xp_4(1)+rand1
          f(inject_loc(i))%x(2)=dummy_xp_4(2)+rand2
          f(inject_loc(i))%x(3)=dummy_xp_4(3)+rand3
          if (i==1) then
            f(inject_loc(i))%behind=inject_loc(dummy_inject_size)
            f(inject_loc(i))%infront=inject_loc(i+1)
          else if (i==dummy_inject_size) then 
            f(inject_loc(i))%behind=inject_loc(i-1)
            f(inject_loc(i))%infront=inject_loc(1)
          else
            f(inject_loc(i))%behind=inject_loc(i-1)
            f(inject_loc(i))%infront=inject_loc(i+1)
          end if
          !zero the stored velocities
          f(inject_loc(i))%u1=0. ; f(inject_loc(i))%u2=0.
        end do
      case default
        call fatal_error('inject.mod','inject_type not set to useable value')    
    end select
    !finally see if we need to switch the injection direction
    if (mod(itime,inject_freq)==0) then
      !OK we can change injection direction
      if (injection_direction==1) then
        injection_direction=2
      else if (injection_direction==2) then
        injection_direction=3
      else if (injection_direction==3) then
        injection_direction=1
      end if      
    end if
    deallocate(inject_loc) !deallocate the location array
  end subroutine
end module
