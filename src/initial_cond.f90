!> Module which contains initial conditions for filament (set in run.in) for more 
!!information on the options see \ref INIT
module initial_cond
  use cdata
  use general
  use periodic
  contains
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
    !*************************************************************************
  !>set up an ellips in the x-y plane, it's size is dictated by the initial 
  !>number of particles set and the size of \f$\delta\f$ 
  !>\todo put eccentricity in run.in
  subroutine setup_ellipse
    implicit none
    real :: velocity
    real :: radius! not really a radius just a helper
    real, parameter :: elip_a=1.5, elip_b=0.5
    integer :: i 
    radius=(0.5*pcount*delta)/(2*pi) !50% of potential size 
    velocity=(quant_circ/(4*pi*radius))*log(8E8*radius)
    write(*,*) 'initf: ellipse, parameters a/b:', elip_a, elip_b
    write(*,*) 'eccentricity is:', sqrt((elip_a**2-elip_b**2)/elip_a**2)
    
    !loop over particles setting spatial and 'loop' position
    do i=1, pcount
      f(i)%x(1)=elip_a*radius*sin(pi*real(2*i-1)/pcount)
      f(i)%x(2)=elip_b*radius*cos(pi*real(2*i-1)/pcount)
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
      f(i)%x(1)=0. !box_size/3
      f(i)%x(2)=0. !-box_size/3
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
  !>draw hypotrochoid (http://en.wikipedia.org/wiki/Hypotrochoid) in the x-y plane
  !> with sinusoidal variation in z
  subroutine setup_hypotrochoid
    implicit none
    real :: velocity
    real :: radius, innerr
    real :: theta
    integer :: i 
    radius=(0.75*pcount*delta)/(4*pi) !75% of potential size 
    innerr=3*radius/5. !radius of inner circle
    write(*,*) 'initf: hypotrochoid'
    !loop over particles setting spatial and 'loop' position
    do i=1, pcount
      theta=3*pi*(2.*i-1.)/pcount
      f(i)%x(1)=(radius-innerr)*cos(theta)+radius*cos(theta*(radius-innerr)/innerr)
      f(i)%x(2)=(radius-innerr)*sin(theta)-radius*sin(theta*(radius-innerr)/innerr)
      f(i)%x(3)=10*delta*sin(4*theta)
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
  !>lines from the lop of the box to the bottom arranged in a lattice, the number of
  !>lines should be a square number
  subroutine setup_lattice
    implicit none
    real :: xpos, ypos
    integer :: pcount_required
    integer :: line_size, line_position
    integer :: counter=0
    integer :: i, j, k
    !test run.in parameters, if wrong program will exit
    if (line_count==0) then
      call fatal_error('init.mod:setup_lattice', &
      'you have not set a value for line_count in run.in')
    end if
    if (periodic_bc) then
      !work out the number of particles required for single line
      !given the box size specified in run.i
      pcount_required=line_count*nint(box_size/(.75*delta))
      write(*,*) 'changing size of pcount to fit with box_length and delta'
      write(*,'(a,i7.1)') ' pcount is now:', pcount_required
      deallocate(f) ; pcount=pcount_required ; allocate(f(pcount))
    else
      call fatal_error('init.mod:setup_lattice', &
      'periodic boundary conditions required')
    end if
    write(*,'(a,i3.1,a)') ' drawing', line_count, ' lines from -z to +z'
    write(*,'(a,f6.3,a)') ' lines in a lattice design occupying ', lattice_ratio, ' of box area'
    line_size=int(pcount/line_count)
    !START THE LOOP
    do i=1, floor(sqrt(real(line_count))) ; do k=1, floor(sqrt(real(line_count)))
      xpos=(-box_size/2.+box_size*((2.*i-1.)/(2*sqrt(real(line_count)))))*lattice_ratio
      ypos=(-box_size/2.+box_size*((2.*k-1.)/(2*sqrt(real(line_count)))))*lattice_ratio
      do j=1, line_size
        line_position=j+counter*line_size
        f(line_position)%x(1)=xpos
        f(line_position)%x(2)=ypos
        f(line_position)%x(3)=-box_size/2.+box_size*real(2*j-1)/(2.*line_size)
        if(j==1) then
          f(line_position)%behind=(counter+1)*line_size
          f(line_position)%infront=line_position+1
        else if (j==line_size) then
          f(line_position)%behind=line_position-1
          f(line_position)%infront=counter*line_size+1
        else
          f(line_position)%behind=line_position-1
          f(line_position)%infront=line_position+1
        end if
        f(line_position)%u1=0. ; f(line_position)%u2=0.
      end do
      counter=counter+1
    end do ; end do
    if (counter/=line_count) then
      call fatal_error('init.mod:setup_lattice', &
      'line_count must be a square number')  
    end if
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
    write(*,'(a,i5.1,a)') ' drawing ', line_count, ' random loops in the box'
    write(*,'(a,i5.1,a)') ' each loop contains ', loop_size, ' particles'
    write(*,'(a,f7.4)') ' radius of each loop: ', loop_radius
    if (rotation_factor<1.) then
      write(*,'(a,f7.4,a)') ' rotation applied to loops: ', rotation_factor, '*2\pi'
    end if
    if (lattice_ratio<1.) then
      write(*,'(a,f7.4,a)') ' loops occupy ', 100*lattice_ratio, '% of box volume'
    end if
    do i=1, line_count
      call random_number(anglex)
      call random_number(angley)
      call random_number(anglez)
      call random_number(translate)
      anglex=anglex*2*pi*rotation_factor
      angley=angley*2*pi*rotation_factor
      anglez=anglez*2*pi*rotation_factor
      translate=lattice_ratio*((box_size*translate-box_size/2.)-loop_radius)
        
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
!*************************************************************************
  !>set up a single loop in the x-y plane, it's size is dictated by the initial 
  !>number of particles set and the size of \f$\delta\f$, this matches with the 
  !>initial SPH condition SPH_loop
  subroutine setup_SPH_loop
    implicit none
    real :: radius
    integer :: i 
    !check that the velocity field is SPH
    select case(velocity)
      case('SPH')
        !we are OK
      case default
        call fatal_error('init.mod','must have SPH velocity field set for this initial condition')
    end select
    !check that the initial SPH condition is a loop
    select case(SPH_init)
      case('loop')
        !we are OK
      case default
        call fatal_error('init.mod','must have SPH initial condition set to loop')
    end select
    !finally check that the initial number of particles are the same
    if (pcount/=SPH_count) then
      call fatal_error('init.mod','must have pcount=SPH_count')
    end if
    radius=(0.75*pcount*delta)/(2*pi) !75% of potential size 
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
      !set the sph particle
      f(i)%sph=i
    end do   
  end subroutine
end module
!******************************************************************
!>\page INIT Initial condition for filament
!!Initial condition for filament is set in run.in through the parameter initf\n
!!Options are:\n
!!- \p single_loop - Single loop in the x-y plane (z=0). The size of the loop is dictated by the number of intial particles set.
!!\image html single_loop_thumb.png
!!- \p single_line - Single line from -z to +z (x=y=0). The number of particles needed is set by the box size and delta, pcount is automatically adjusted to fit this.
!!\image html single_line_thumb.png
!!- \p random_loops - multiple random loops (number set by the parameter \p line_count) distributed throughout the box. The size of the loops is dictated by the number of loops and the number of intial particles set. You must ensure \p pcount is a multiple of the number of loops
!!\image html random_loops_thumb.png
!!- \p crow - Two lines from -z to +z (x=y=0), anti-parallel to trigger the
!!crow instability. The number of particles needed is set by the box size and
!!delta, pcount is automatically adjusted to fit this.
!!\image html crow_thumb.png
!!- \p leap-frog - two loops in the x-y plane at z=-delta, z=delta this sets
!!the loops leapfrogging (not with LIA). 
!!\image html leap-frog_thumb.png
!!- \p criss-cross - random lines criss-crossing the box
!!\image html criss-cross_thumb.png
!!- \p tangle - curved lines from -z to +z which will drive a tangle, 
!!superceeded in many ways by criss-cross but kept for posterity.
!!\image html tangle_thumb.png
!!- \p lattice - lines from -z to +z in a lattice desing, 
!!line_count must be a square number, the area of the box the lattice occupies
!!can be changed with the lattice_ratio variable in run.in.
!!\image html lattice_thumb.png
!!- \p wave_line -  lines from -z to z with wave pertubations added, these can
!!be either helical or planar with an imposed spectra and wave count. Relevant
!!parameters in run.in are wave_count, wave_slope (spectral slope), wave_amp 
!!(amplitude of pertubations) and wave_type (planar or helical). Note there is
!!also the wave_loop intial condition which is very similar.
!!\image html wave_line_thumb.png

