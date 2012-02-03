!> Module which contains initial conditions (loops) for filament (set in run.in) for more 
!!information on the options see \ref INIT
module initial_loop
  use cdata
  use general
  use periodic
  contains
  !*************************************************************************
  !>set up a single loop in the x-y plane, it's size is dictated by the initial 
  !>number of particles set and the size of \f$\delta\f$
  subroutine setup_single_loop
    implicit none
    real :: velocity, ring_energy
    real :: radius
    integer :: i 
    radius=(0.75*pcount*delta)/(2*pi) !75% of potential size 
    velocity=(quant_circ/(4*pi*radius))*(log(8*radius/corea)-.5)
    ring_energy=0.5*(quant_circ**2)*radius*(log(8*radius/corea)-2.)
    write(*,*) 'initf: single loop, radius of loop:', radius
    write(*,*) 'velocity should be:', velocity
    write(*,*) 'energy should be:', ring_energy
    !loop over particles setting spatial and 'loop' position
    do i=1, pcount
      f(i)%x(1)=radius*sin(pi*real(2*i-1)/pcount)
      f(i)%x(2)=radius*cos(pi*real(2*i-1)/pcount) !+ box_size/2.
      f(i)%x(3)=0.
      if (i==1) then
        f(i)%behind=pcount ; f(i)%infront=i+1
      else if (i==pcount) then 
        f(i)%behind=i-1 ; f(i)%infront=1
      else
        f(i)%behind=i-1 ; f(i)%infront=i+1
      end if
      !zero the stored velocities
      f(i)%u1=0. ; f(i)%u2=0. ; f(i)%u3=0.
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
      f(i)%u1=0. ; f(i)%u2=0. ; f(i)%u3=0.
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
      f(i)%u1=0. ; f(i)%u2=0. ; f(i)%u3=0.
    end do   
  end subroutine
  !*************************************************************************
  !>draw torus knot (http://en.wikipedia.org/wiki/Torus_knot) in the x-y plane
  !> asymptotic form from Ricca, Chaos [1993]. In cylindrical polars \f$(r,\theta,z)\f$ 
  !>\f[ r=r_0+\epsilon \sin(w\phi)\f] \f[ \theta=\phi+\frac{\epsilon}{wr_0} \cos(w\phi)\f] \f[ z=\epsilon (1+w^{-2})^{1/2} cos(w\phi)\f] 
  !> where \f$ w=q/p\f$ is the winding number, p, q and \f$\epsilon \ll 1\f$
  !!are set in run.in
  subroutine setup_torus_knot
    implicit none
    integer :: i 
    real :: phi, w, radius
    real :: rho, theta !cylindrical co-ordinates
    radius=(0.75*pcount*delta)/(6*pi) !75% of potential size 
    write(*,'(a,i4.2,a,i4.2)') 'initf: torus knot p= ', torus_p, ' q= ', torus_q
    w=real(torus_q)/real(torus_p)
    write(*,'(a,f4.3,a,f6.4)') 'winding number ', w, ' loop radius ', radius 
    do i=1, pcount
      phi=torus_p*2.*pi*(2.*i-1.)/(2.*pcount) !phase \in [0,2*pi]
      rho=radius+torus_epsilon*sin(w*phi) !radius
      theta=phi+(torus_epsilon/(w*radius))*cos(w*phi)
      f(i)%x(1)=rho*cos(theta)
      f(i)%x(2)=rho*sin(theta)
      f(i)%x(3)=torus_epsilon*((1+1./w)**(0.5))*cos(w*phi)
      if (i==1) then
        f(i)%behind=pcount ; f(i)%infront=i+1
      else if (i==pcount) then 
        f(i)%behind=i-1 ; f(i)%infront=1
      else
        f(i)%behind=i-1 ; f(i)%infront=i+1
      end if
      !zero the stored velocities
      f(i)%u1=0. ; f(i)%u2=0. ; f(i)%u3=0.
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
      f(i)%u1=0. ; f(i)%u2=0. ; f(i)%u3=0.
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
      f(i)%u1=0. ; f(i)%u2=0. ; f(i)%u3=0.
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
      f(i)%u1=0. ; f(i)%u2=0. ; f(i)%u3=0.
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
      f(i)%u1=0. ; f(i)%u2=0. ; f(i)%u3=0.
    end do
    !second loop
    do i=pcount/2+1, pcount
      f(i)%x(1)=radius*sin(pi*real(2.*(i-pcount/2)-1)/(pcount/2))+radius/0.55
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
      f(i)%u1=0. ; f(i)%u2=0. ; f(i)%u3=0.
    end do    
  end subroutine
  !*************************************************************************
  !>two unlinked loops which drive a reconnection via crow instability
  subroutine setup_crow_loop
    !two linked loops 
    implicit none
    real :: radius
    integer :: i 
    if (mod(pcount,2)/=0) then
      call fatal_error('init.mod:setup_crow_loop', &
      'pcount is not a multiple of 2-aborting')
    end if
    radius=(0.75*pcount*delta)/(4*pi) !75% of potential size
    write(*,'(a,f9.6)') ' initf: linked filaments, radius of loops:', radius 
    !loop over particles setting spatial and 'loop' position
    do i=1, pcount/2
      f(i)%x(1)=radius*sin(pi*real(2*i-1)/(pcount/2))
      f(i)%x(2)=0.
      f(i)%x(3)=radius*cos(pi*real(2*i-1)/(pcount/2))
      if (i==1) then
        f(i)%behind=pcount/2 ; f(i)%infront=i+1
      else if (i==pcount/2) then 
        f(i)%behind=i-1 ; f(i)%infront=1
      else
        f(i)%behind=i-1 ; f(i)%infront=i+1
      end if
      !zero the stored velocities
      f(i)%u1=0. ; f(i)%u2=0. ; f(i)%u3=0.
    end do
    !second loop
    do i=pcount/2+1, pcount
      f(i)%x(1)=radius*sin(pi*real(2.*(i-pcount/2)-1)/(pcount/2))+2.1*radius
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
      f(i)%u1=0. ; f(i)%u2=0. ; f(i)%u3=0.
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
      f(i)%u1=0. ; f(i)%u2=0. ; f(i)%u3=0.
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
      f(i)%u1=0. ; f(i)%u2=0. ; f(i)%u3=0.
    end do    
  end subroutine
  !*************************************************************************
  !>set up a single loop in the y-z plane, it's size is dictated by the initial 
  !>number of particles set and the size of \f$\delta\f$
  subroutine setup_single_loop_zy
    implicit none
    real :: velocity, ring_energy
    real :: radius
    integer :: i
    radius=(0.75*pcount*delta)/(2*pi) !75% of potential size 
    velocity=(quant_circ/(4*pi*radius))*(log(8*radius/corea)-.5)
    ring_energy=0.5*(quant_circ**2)*radius*(log(8*radius/corea)-2.)
    write(*,*) 'initf: single loop zy, radius of loop:', radius
    write(*,*) 'velocity should be:', velocity
    write(*,*) 'energy should be:', ring_energy
    !loop over particles setting spatial and 'loop' position
    do i=1, pcount
      f(i)%x(1)=-box_size/2.
      f(i)%x(2)=radius*cos(pi*real(2*i-1)/pcount)
      f(i)%x(3)=radius*sin(pi*real(2*i-1)/pcount)
      if (i==1) then
        f(i)%behind=pcount ; f(i)%infront=i+1
      else if (i==pcount) then
        f(i)%behind=i-1 ; f(i)%infront=1
      else
        f(i)%behind=i-1 ; f(i)%infront=i+1
      end if
      !zero the stored velocities
      f(i)%u1=0. ; f(i)%u2=0. ; f(i)%u3=0.
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
      f(i)%u1=0. ; f(i)%u2=0. ; f(i)%u3=0.
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
      f(i)%u1=0. ; f(i)%u2=0. ; f(i)%u3=0.
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
      f(i)%u1=0. ; f(i)%u2=0. ; f(i)%u3=0.
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
      f(i)%u1=0. ; f(i)%u2=0. ; f(i)%u3=0.
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
    write(*,'(a,f9.5,a)') ' loops scattered by ', lattice_ratio, ' in z direction'
    write(*,'(a,i3.1)') ' starting wavenumber ', wave_start
    write(*,'(a,i3.1)') ' wavenumber separation', wave_skip
    write(*,'(i4.1,a,a,a,f9.5)') wave_count, ' ',trim(wave_type),' wave pertubations, with spectral slope:', wave_slope
    loop_size=int(pcount/line_count)
    loop_radius=loop_size*(0.75*delta)/(2*pi) !75% of potential size
    write(*,'(a,f9.3)') ' loop radius: ', loop_radius
    !START THE LOOP
    do i=1, line_count
      call random_number(zpos)
      zpos=(zpos-.5)*box_size*lattice_ratio !lattice ratio sets scattering
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
        f(loop_position)%u1=0. ; f(loop_position)%u2=0. ; f(loop_position)%u3=0.
      end do
      !we have now drawn the basic loop, now we add the wave pertubations
      prefactor=wave_amp/(wave_start**wave_slope) !our starting wavenumber is wave_start
      if (i==1) then !only write this once
        write(*,*)'wave information recorded in ./data/wave_info.log'
        open(unit=34,file='./data/wave_info.log')
        write(34,*) '%------------------WAVE INFORMATION-------------------'
      end if
      do k=1, wave_count !wave_count set in run.in
        call ghostp !we must call this routine at the start of adding every wave 
                    !the routines normalf, binormalf rely on correct ghostpoints
        !on a loop so wavenumbers must be an integer
        wave_number=wave_start+(k-1)*wave_skip
        amp=prefactor*(wave_number**wave_slope)
        call random_number(random_shift) !help things along with a 
        random_shift=random_shift*2*pi   !random shift \in (0,2\pi)
        !random_shift=0.
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
  !>a macroscopic ring created by many quantised loops with a hexagonal pattern
  subroutine setup_macro_ring
    implicit none
    real :: zpos
    real, allocatable :: rad_shift(:), x_shift(:)
    real, allocatable :: rad_shift_r(:), rad_shift_theta(:)
    real :: average_loop_sep
    integer :: shift_help, shift_counter
    integer :: actual_line_count, pcount_required
    integer :: loop_size, loop_position
    integer :: i, j
    !test run.in parameters, if wrong program will exit
    if (line_count==0) then
      call fatal_error('init.mod:setup_macro_ring', &
      'you have not set a value for line_count in run.in')
    end if
    if ((macro_ring_R<epsilon(0.)).or.(macro_ring_a<epsilon(0.))) then
      call fatal_error('init.mod:setup_macro_ring', &
      'you have not set a value for either macro_ring_a or R in run.in')
    end if
    if (macro_ring_a>macro_ring_R) then
      call warning_message('init.mod','major axis of torus is less than minor axis?')
    end if
    !find the hexagonal number
    actual_line_count=3*line_count*(line_count-1)+1
    write(*,'(a,i3.1,a)') ' creating a macro ring with ', line_count, ' layer(s)'
    write(*,'(a,f6.3,a,f6.3)') ' major radius ', macro_ring_R, ' minor radius ', macro_ring_a
    write(*,'(a,i3.1,a)') ' drawing', actual_line_count, ' loops in the z-y plane'
    !how big does a loop need to be?
    loop_size=nint(2*pi*macro_ring_R/(0.75*delta))
    pcount_required=loop_size*actual_line_count !# particles needed
    print*, 'changing size of pcount to fit with major_R, delta and line_count'
    print*, 'pcount is now', pcount_required
    average_loop_sep=sqrt(pi)*macro_ring_a/sqrt(real(actual_line_count))
    if  (average_loop_sep< delta/2.) then
      call fatal_error('init.mod','separation of loops is less than resolution')
    end if
    deallocate(f) ; pcount=pcount_required ; allocate(f(pcount))
    !now we need to find shift values to create the torus
    allocate(rad_shift(actual_line_count),x_shift(actual_line_count))
    allocate(rad_shift_r(actual_line_count),rad_shift_theta(actual_line_count))
    rad_shift_r(1)=0. ; rad_shift_theta(1)=0.
    shift_counter=1
    if (line_count>1) then !we may not need this if statement
      do i=2,line_count
        do j=1,((i-1)*6)
          shift_counter=shift_counter+1
          rad_shift_r(shift_counter)=macro_ring_a*(real(i-1)/line_count)
          rad_shift_theta(shift_counter)=pi*(2.*j-1.)/((i-1)*6)
        end do
      end do
    end if
    open(unit=37,file='./data/macro_loop_profile.log')
    do i=1, actual_line_count
      rad_shift(i)=rad_shift_r(i)*cos(rad_shift_theta(i))
      x_shift(i)=rad_shift_r(i)*sin(rad_shift_theta(i))
      write(37,*) i, rad_shift(i), x_shift(i)
    end do
   close(37)
    !START THE LOOP
    do i=1, actual_line_count
      do j=1, loop_size
        loop_position=j+(i-1)*loop_size
        f(loop_position)%x(1)=-box_size/2.+1.5*macro_ring_a+x_shift(i)
        f(loop_position)%x(2)=(macro_ring_R+rad_shift(i))*cos(pi*real(2*j-1)/loop_size)
        f(loop_position)%x(3)=(macro_ring_R+rad_shift(i))*sin(pi*real(2*j-1)/loop_size)
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
        f(loop_position)%u1=0. ; f(loop_position)%u2=0. ; f(loop_position)%u3=0.
      end do
    end do 
    deallocate(rad_shift,x_shift,rad_shift_r,rad_shift_theta)
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
          f(loop_position)%u1=0. ; f(loop_position)%u2=0. ; f(loop_position)%u3=0.
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
          f(loop_position)%u1=0. ; f(loop_position)%u2=0. ; f(loop_position)%u3=0.
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
  subroutine setup_random_loops
    implicit none
    real :: loop_radius
    real :: anglex,angley,anglez
    real,dimension(3)::translate
    real,dimension(3)::dummy_xp_1, dummy_xp_2, dummy_xp_3, dummy_xp_4
    integer :: mean_loop_size, loop_size
    integer:: loop_position
    integer :: used_pcount, total_loop_count=0
    integer :: i,j
    !test run.in parameters, if wrong program will exit
    if (line_count==0) then
      call fatal_error('init.mod:setup_random_loops', &
      'you have not set a value for line_count in run.in')
    end if
    !!\todo may be able to remove this constraint
    if (mod(pcount,line_count)/=0) then
      call fatal_error('init.mod:setup_random_loops', &
      'pcount/line_count is not an integer')
    end if
    mean_loop_size=int(pcount/line_count)
    loop_radius=mean_loop_size*(0.75*delta)/(2*pi) !75% of potential size
    if (line_sigma>0.) then
      write(*,'(a,a,a)') ' loop sizes follow a ', trim(initial_distribution), ' distribution'
      write(*,'(a,i5.1,a)') ' mean particle count per loop ', mean_loop_size, ' particles'
      write(*,'(a,f8.4)') ' standard deviation of loop sizes: ', line_sigma
    else
      write(*,'(a,i5.1,a)') ' each loop contains ', mean_loop_size, ' particles'
      write(*,'(a,f7.4)') ' radius of each loop: ', loop_radius
    end if
    if (rotation_factor<1.) then
      write(*,'(a,f7.4,a)') ' rotation applied to loops: ', rotation_factor, '*2\pi'
    end if
    if (maxval(loop_translate)>1.) then
      call fatal_error('init.mod','max of loop_translate >1')
    end if
    if (minval(loop_translate)<1.) then
      write(*,'(a,3f7.2,a)') ' loops translated by (x,y,z)', 100*loop_translate, ' % of box volume'
    end if
    used_pcount=0 !zero this before we start do loop
    open(unit=79,file='data/random_loop_sizes.log',status='replace')
    do while (used_pcount<pcount)
      call random_number(anglex)
      call random_number(angley)
      call random_number(anglez)
      call random_number(translate)
      anglex=anglex*2*pi*rotation_factor
      angley=angley*2*pi*rotation_factor
      anglez=anglez*2*pi*rotation_factor
      translate=loop_translate*((box_size*translate-box_size/2.)-loop_radius)
      !set the loop size using specified probability distribution
      select case(initial_distribution)
        case('uniform')
          loop_size=nint(runif(real(mean_loop_size)-sqrt(3.)*line_sigma,real(mean_loop_size)+sqrt(3.)*line_sigma))
        case('gaussian')
          loop_size=nint(rnorm(real(mean_loop_size),line_sigma**2))
        case('laplace')  
          loop_size=nint(rlaplace(real(mean_loop_size),line_sigma))
        case('gamma')  
           loop_size=nint(rgamma((real(mean_loop_size)/line_sigma)**2,(line_sigma**2)/real(mean_loop_size)))
        case('exp')  
          loop_size=nint(rexp(1./real(mean_loop_size)))
        case default
          call fatal_error('init_cond.mod','initial_distribution set to incorrect argument')  
      end select
      if (loop_size<0) loop_size=0 !avoid negative loop sizes
      if (used_pcount+loop_size>pcount) then 
        loop_size=pcount-used_pcount
      end if
      loop_radius=loop_size*(0.75*delta)/(2*pi) !75% of potential size
      write(79,*) loop_size,loop_radius
      do j=1, loop_size

        loop_position=j+used_pcount

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
          f(loop_position)%behind=used_pcount+loop_size	
          f(loop_position)%infront=loop_position+1
        else if (j==loop_size) then
          f(loop_position)%behind=loop_position-1
          f(loop_position)%infront=used_pcount+1
        else
          f(loop_position)%behind=loop_position-1
          f(loop_position)%infront=loop_position+1
        end if
        f(loop_position)%u1=0. ; f(loop_position)%u2=0. ; f(loop_position)%u3=0.
        !test the loop size if its too small then we simply kill the loop here
        if (loop_size<5) then
          call clear_particle(loop_position)
        end if
      end do
      used_pcount=used_pcount+loop_size
      total_loop_count=total_loop_count+1
    end do
    close(79)
    write(*,'(i5.1,a)') total_loop_count, ' random loops have been created'
  end subroutine
  !*************************************************************************
  subroutine setup_loop_train
    implicit none
    integer :: pcount_required
    integer :: loop_size, loop_position
    real :: loop_radius, z_pos
    integer :: i,j
    !test run.in parameters, if wrong program will exit
    if (line_count==0) then
      call fatal_error('init.mod:setup_loop_train', &
      'you have not set a value for line_count in run.in')
    end if
    if (periodic_bc) then
      !work out the number of particles required for our lines
      !given the box size specified in run.in
      loop_radius=0.45*box_size
      loop_size=nint(loop_radius*2*pi/(0.75*delta)) !75%
      pcount_required=line_count*loop_size
      write(*,*) 'changing size of pcount to fit with box_length and delta'
      write(*,*) 'pcount is now', pcount_required
      deallocate(f) ; pcount=pcount_required ; allocate(f(pcount))
    else
      call fatal_error('init.mod:setup_loop_train',&
                       'periodic boundary conditions required for this initial condition') !cdata.mod
    end if
    do i=1, line_count
      z_pos=(2*i-1)/(2.*line_count)*box_size-box_size/2.
      do j=1, loop_size
        loop_position=j+(i-1)*loop_size
        if (mod(i,2)==0) then
          f(loop_position)%x(1)=loop_radius*sin(pi*real(2*j-1)/loop_size)
          f(loop_position)%x(2)=loop_radius*cos(pi*real(2*j-1)/loop_size)+box_size/2.
        else
          f(loop_position)%x(1)=loop_radius*sin(pi*real(2*j-1)/loop_size)+box_size/2.
          f(loop_position)%x(2)=loop_radius*cos(pi*real(2*j-1)/loop_size)
        end if
        f(loop_position)%x(3)=z_pos
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
        f(loop_position)%u1=0. ; f(loop_position)%u2=0. ; f(loop_position)%u3=0.
      end do
    end do
  end subroutine
  !*************************************************************************
  subroutine setup_loop_stream
    implicit none
    integer :: pcount_required
    integer :: loop_size, loop_position
    real :: loop_radius, z_pos
    integer :: i,j
    !test run.in parameters, if wrong program will exit
    if (line_count==0) then
      call fatal_error('init.mod:setup_loop_train', &
      'you have not set a value for line_count in run.in')
    end if
    if (periodic_bc) then
      !work out the number of particles required for our lines
      !given the box size specified in run.in
      loop_radius=0.245*box_size
      loop_size=nint(loop_radius*2*pi/(0.75*delta)) !75%
      pcount_required=line_count*loop_size
      write(*,*) 'changing size of pcount to fit with box_length and delta'
      write(*,*) 'pcount is now', pcount_required
      deallocate(f) ; pcount=pcount_required ; allocate(f(pcount))
    else
      call fatal_error('init.mod:setup_loop_train',&
                       'periodic boundary conditions required for this initial condition') !cdata.mod
    end if
    do i=1, line_count
      z_pos=(2*i-1)/(2.*line_count)*box_size-box_size/2.
      do j=1, loop_size
        loop_position=j+(i-1)*loop_size
        if (mod(i,2)==0) then
          f(loop_position)%x(1)=loop_radius*sin(pi*real(2*j-1)/loop_size)+box_size/5
          f(loop_position)%x(2)=loop_radius*cos(pi*real(2*j-1)/loop_size)+box_size/5.
        else
          f(loop_position)%x(1)=loop_radius*sin(pi*real(2*j-1)/loop_size)-box_size/5.
          f(loop_position)%x(2)=-loop_radius*cos(pi*real(2*j-1)/loop_size)-box_size/5.
        end if
        f(loop_position)%x(3)=z_pos
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
        f(loop_position)%u1=0. ; f(loop_position)%u2=0. ; f(loop_position)%u3=0.
      end do
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
!!- \p torus_knot - A torus knot in the xy-plane the winding number w=p/q can
!!be set (via torus_p and torus_q) in run.in.
!!\image html knot_thumb.png
