!>The particles aspect of the code. At present this evolves either fluid or
!>inertial particles in the code. The type of particle is set by the
!>particle_type parameter in run.in Fluid/Inertial particles use the 
!>normal fluid velocity. The initial setup of the particles is set throught the initp parameter. See \ref PARTICLES for more details
module particles
  use cdata
  use general
  use timestep
  use output
  use biofluid
  real :: part_maxu=0., part_maxdu=0., part_urms=0. !velocity information
  real :: part_sep=0. !average particle separation
  contains
  !************************************************************
  !>setup the particles in as set by initp
  subroutine setup_particles
    implicit none
    integer :: i, j, k
    integer :: counter
    real :: rand1, rand2, rand3, sphere_theta, sphere_phi
    allocate(p(part_count))
    write(*,*) 'setting up particles in ', trim(initp),' configuration'
    write(*,'(a,a,a,i4.1)') ' particle type: ', trim(particle_type), ', particle count: ', part_count
    select case(initp)
      !**************************FLUID/INERITAL-PARTICLES*************************
      !put in all the fluid/inertial particle initial conditions in the code here
      case('one-side')
        !particles all on one side of the box (x-axis) 
        do i=1, part_count
          p(i)%x(1)=-box_size/2.
          call random_number(p(i)%x(2)) ; call random_number(p(i)%x(3))
          p(i)%x(2)=box_size*p(i)%x(2)-box_size/2.
          p(i)%x(3)=box_size*p(i)%x(3)-box_size/2.
          p(i)%u=0. ; p(i)%u1=0. ; p(i)%u2=0. !0 the velocity fields
        end do
      case('random')
        !particles in random positions
        do i=1, part_count
          call random_number(p(i)%x(1))
          call random_number(p(i)%x(2)) ; call random_number(p(i)%x(3))
          if (box_size>0.) then
            p(i)%x(1)=box_size*p(i)%x(1)-box_size/2.
            p(i)%x(2)=box_size*p(i)%x(2)-box_size/2.
            p(i)%x(3)=box_size*p(i)%x(3)-box_size/2.
          end if
          p(i)%u=0. ; p(i)%u1=0. ; p(i)%u2=0. !0 the velocity fields
        end do
      case('ring')
        !particles in random positions
        do i=1, part_count
          p(i)%x(1)=0.75*box_size/2.*cos(2*pi*(2.*i-1.)/(2.*part_count))
          p(i)%x(2)=0.75*box_size/2.*sin(2*pi*(2.*i-1.)/(2.*part_count))
          !p(i)%x(3)=-box_size/2.
          p(i)%x(3)=0.
          p(i)%u=0. ; p(i)%u1=0. ; p(i)%u2=0. !0 the velocity fields
        end do
      case('tube')
        write(*,'(a,f6.4,a)') ' tube radius: ', part_tube_ratio, ' of box size'
        !particles in random positions
        do i=1, part_count
          rand1=runif(0.,box_size*part_tube_ratio)
          rand2=runif(0.,2.*pi)
          rand3=runif(-box_size/2.,box_size/2.)
          p(i)%x(1)=rand1*cos(rand2)
          p(i)%x(2)=rand1*sin(rand2)
          p(i)%x(3)=rand3
          p(i)%u=0. ; p(i)%u1=0. ; p(i)%u2=0. !0 the velocity/Users/abaggaley/Work/code/qvort fields
        end do
      case('crow_tube')
        write(*,'(a,f6.4,a)') ' tube radius: ', part_tube_ratio, ' of box size'
        do i=1, part_count/2
          rand1=runif(0.,box_size*part_tube_ratio)
          rand2=runif(0.,2.*pi)
          p(i)%x(1)=rand1*cos(rand2)
          p(i)%x(3)=-box_size/2.+box_size*real(2*i-1)/(part_count)
          p(i)%x(2)=delta-(delta/16.)*sin(pi*(box_size/2.+p(i)%x(3))/box_size)+rand1*sin(rand2)
          p(i)%u=0. ; p(i)%u1=0. ; p(i)%u2=0. !0 the velocity fields
        end do
        do i=part_count/2+1, part_count
          rand1=runif(0.,box_size*part_tube_ratio)
          rand2=runif(0.,2.*pi)
          p(i)%x(1)=rand1*cos(rand2)
          p(i)%x(3)=box_size/2.-box_size*real(2*(i-part_count/2)-1)/(part_count)
          p(i)%x(2)=-delta+(delta/16.)*sin(pi*(box_size/2.-p(i)%x(3))/box_size)+rand1*sin(rand2)
          p(i)%u=0. ; p(i)%u1=0. ; p(i)%u2=0. !0 the velocity fields
        end do
      case('orthog_tube')
        write(*,'(a,f6.4,a)') ' tube radius: ', part_tube_ratio, ' of box size'
        do i=1, part_count/2
          rand1=runif(0.,box_size*part_tube_ratio)
          rand2=runif(0.,2.*pi)
          p(i)%x(1)=rand1*cos(rand2)
          p(i)%x(3)=-box_size/2.+box_size*real(2*i-1)/(part_count)
          p(i)%x(2)=0.5*delta+rand1*sin(rand2)
          p(i)%u=0. ; p(i)%u1=0. ; p(i)%u2=0. !0 the velocity fields
        end do
        do i=part_count/2+1, part_count
          rand1=runif(0.,box_size*part_tube_ratio)
          rand2=runif(0.,2.*pi)
          p(i)%x(1)=box_size/2.-box_size*real(2*(i-pcount/2)-1)/(part_count)
          p(i)%x(2)=-0.5*delta+rand1*sin(rand2)
          p(i)%x(3)=rand1*cos(rand2)
          p(i)%u=0. ; p(i)%u1=0. ; p(i)%u2=0. !0 the velocity fields
        end do
      case('central_sphere')
        write(*,'(a,f6.4,a)') ' sphere radius: ', part_sphere_radius, ' of box size'
        do i=1, part_count
          rand1=runif(0.,part_sphere_radius*box_size) !r - distance from centre
          rand2=runif(0.,1.)
          rand3=runif(0.,1.)
          sphere_theta=2.*pi*rand2
          sphere_phi=acos(2.*rand3-1.)
          p(i)%x(1)=rand1*cos(sphere_theta)*sin(sphere_phi)
          p(i)%x(2)=rand1*sin(sphere_theta)*sin(sphere_phi)
          p(i)%x(3)=rand1*cos(sphere_phi)
          p(i)%u=0. ; p(i)%u1=0. ; p(i)%u2=0. !0 the velocity fields
        end do
      case('pairs')
        !particles in pairs in a random position
        !check that the number of particles is a multiple of 2
        if (mod(part_count,2)/=0) then
          call fatal_error('setup_particles','part_count must be a multiple of 2')
        end if
        do i=1,  part_count, 2
          call random_number(p(i)%x(1))
          call random_number(p(i)%x(2)) ; call random_number(p(i)%x(3))
          p(i)%x(1)=box_size*p(i)%x(1)-box_size/2.
          p(i)%x(2)=box_size*p(i)%x(2)-box_size/2.
          p(i)%x(3)=box_size*p(i)%x(3)-box_size/2.
          p(i)%u=0. ; p(i)%u1=0. ; p(i)%u2=0. !0 the velocity fields
          !finally sort out the position of its neighbour
          p(i+1)%x=p(i)%x ; p(i+1)%x(3)=p(i)%x(3)+0.001
          p(i+1)%u=0. ; p(i+1)%u1=0. ; p(i+1)%u2=0. !0 the velocity fields
        end do
      case('lattice')
        counter=0
        do i=1, floor(sqrt(real(part_count))) ; do k=1, floor(sqrt(real(part_count)))
          p((i-1)*floor(sqrt(real(part_count)))+k)%x(1)=(-box_size/2.+box_size*((2.*i-1.)/(2*sqrt(real(part_count)))))*lattice_ratio
          p((i-1)*floor(sqrt(real(part_count)))+k)%x(2)=(-box_size/2.+box_size*((2.*k-1.)/(2*sqrt(real(part_count)))))*lattice_ratio
          p((i-1)*floor(sqrt(real(part_count)))+k)%x(3)=0.
          counter=counter+1
        end do ; end do
        if (counter/=part_count) then
          call fatal_error('init.mod:setup_particles', &
          'part_count must be a square number')  
        end if
      case default
        call fatal_error('setup_particles','initp not set to available value')
    end select
    !check the particles type
    select case(particle_type)
      case('fluid')
        !do nothing
      case('inertial')
        if (part_stokes<epsilon(0.)) call fatal_error('setup_particles',&
        'part_stokes is 0, change particle_type if you want tracers')
        write(*,'(a,f10.4)') 'stokes number is ', part_stokes
      case('swimmers')
        write(*,*) 'swimming particles, G= ', bio_G, ' V= ', bio_V
        !sort out particle directions for swimmers
        do i=1,part_count
          call random_unit_vector(p(i)%p) !general.f90
        end do
      case default
        call fatal_error('setup_particles','particle type incorrect')
    end select
  end subroutine
  !************************************************************
  !>evolve the particles in the code and print to file/perform 
  !>diagnostic tests
  subroutine particles_evolution
    implicit none
    integer :: i
    !timestep the particles
    call timestep_particles
    !particle diagnostics
    if (mod(itime,shots)==0) then
      call diagnostics_particles
    end if
    !print the particles to file
    if (mod(itime,shots)==0) then
      call printp(itime/shots) !output.mod
    end if
    !any other business in here
  end subroutine  
  !***********************************************************
  !>timestep the particles 
  subroutine timestep_particles
    implicit none
    real :: u(3)
    integer :: i
    !$omp parallel do
    do i=1, part_count
      select case(particle_type)
        case('fluid')
          call velocity_fluidp(i,u)
          p(i)%u=u
        case('inertial')
          call velocity_fluidp(i,u)
          !an euler step to get the velocity
          p(i)%u=p(i)%u1+dt*(u-p(i)%u1)*part_stokes
          !adjust the velocity due to stokes drag
        case('swimmers')
          call bio_swimming(i) !biofluid.mod
      end select
      if (maxval(abs(p(i)%u1))==0) then
        p(i)%x(:)=p(i)%x(:)+dt*p(i)%u(:) !euler
      else if (maxval(abs(p(i)%u2))==0) then
        !first order adams-bashforth
        p(i)%x(:)=p(i)%x(:)+three_twos*dt*p(i)%u(:)-one_half*dt*p(i)%u1(:)
      else 
        !second order adams-bashforth
        p(i)%x(:)=p(i)%x(:)+twenty_three_twelve*dt*p(i)%u(:)-four_thirds*dt*p(i)%u1(:)+five_twelths*dt*p(i)%u2(:)
      end if
      !store the old velocities
      p(i)%u2=p(i)%u1
      p(i)%u1=p(i)%u
      !enforce periodicity
      if (periodic_bc) then
        !if a particle leaves one side of the box
        !-------------x------------------     
        if (p(i)%x(1)>(box_size/2.)) then
          p(i)%x(1)=p(i)%x(1)-box_size
        else if (p(i)%x(1)<(-box_size/2.)) then
          p(i)%x(1)=p(i)%x(1)+box_size
        end if
        !-------------y------------------
        if (p(i)%x(2)>(box_size/2.)) then
          p(i)%x(2)=p(i)%x(2)-box_size
        else if (p(i)%x(2)<(-box_size/2.)) then
          p(i)%x(2)=p(i)%x(2)+box_size
        end if
        !-------------z------------------
        if (p(i)%x(3)>(box_size/2.)) then
          p(i)%x(3)=p(i)%x(3)-box_size
        else if (p(i)%x(3)<(-box_size/2.)) then
          p(i)%x(3)=p(i)%x(3)+box_size
        end if
        !--------------------------------
        !put it back in the other side....
      end if
    end do
    !$omp end parallel do
  end subroutine
  !************************************************************
  !> calculate the velocity at each particle
  !> note we can set this in here to be the velocity induced by
  !> the superfluid vortices but this is unphysical
  !> really this should be false by default
  subroutine velocity_fluidp(i,u)
    implicit none 
    integer, intent(IN) :: i
    real, intent(OUT) :: u(3)
    real :: u_sup(3), u_norm(3)
    integer :: peri, perj, perk
    u_sup=0. !must be zeroed for intially
    if (particle_super_velocity) then
      select case(velocity)
        case('LIA','BS')
          call biot_savart_general(p(i)%x,u_sup) !timestep.mod
        case('Tree')
          call tree_walk_general(p(i)%x,vtree,(/0.,0.,0./),u_sup)
          if (periodic_bc) then
            !we must shift the mesh in all 3 directions, all 26 permutations needed!
            do peri=-1,1 ; do perj=-1,1 ; do perk=-1,1
              if (peri==0.and.perj==0.and.perk==0) cycle
              call tree_walk_general(p(i)%x,vtree, &
                   (/peri*box_size,perj*box_size,perk*box_size/),u_sup) !tree.mod
            end do ; end do ;end do
          end if
      end select
    end if
    !normal fluid velocity
    call get_normal_velocity(p(i)%x,u_norm) !normal_fluid.mod
    u=u_sup+u_norm
  end subroutine
  !************************************************************
  !>diagnostics for fluid/inertial particles
  subroutine diagnostics_particles
    implicit none
    real :: uinfo(part_count,2)
    integer :: i
    uinfo(:,1)=sqrt(p(:)%u(1)**2+p(:)%u(2)**2+p(:)%u(3)**2)
    uinfo(:,2)=sqrt((p(:)%u1(1)-p(:)%u2(1))**2+&
                    (p(:)%u1(2)-p(:)%u2(2))**2+&
                    (p(:)%u1(3)-p(:)%u2(3))**2)
    !in here we need to determine
    part_maxu=maxval(uinfo(:,1)) ; part_maxdu=maxval(uinfo(:,2))
    part_urms=sqrt(sum(uinfo(:,1)**2)/part_count)
    !2 particle separation
    select case(initp)
      case('pairs')
        part_sep=0.
        do i=1, part_count, 2
          part_sep=part_sep+sqrt((p(i)%x(1)-p(i+1)%x(1))**2+&
                                 (p(i)%x(2)-p(i+1)%x(2))**2+&
                                 (p(i)%x(3)-p(i+1)%x(3))**2)
        end do
        part_sep=part_sep/(part_count/2)
    end select
    open(unit=78,file='data/par_ts.log',position='append')
    if (itime==shots) then
      write(78,*) '%-var--t-----maxu---maxdu----urms---part_sep--'
      if (particles_only) then
        !printing to screen
        write(*,*) '%-var--t-----maxu---maxdu----urms---part_sep--'
      end if
    end if
    write(78,'(i5.3,f6.2,f8.5,f8.5,f8.5,f8.5)') &
    itime/shots,t,part_maxu,part_maxdu,part_urms,part_sep
    if (particles_only) then
      write(*,'(i5.3,f6.2,f8.5,f8.5,f8.5,f8.5)') &
      itime/shots,t,part_maxu,part_maxdu,part_urms,part_sep
    end if
    close(78)
  end subroutine
end module 
!>\page PARTICLES Fluid/Inertial particles
!!Initial condition for particles is set in run.in through the parameter initp.\n
!!Particle count is set in run.in through the parameter part_count.\n
!!Options are:\n
!!- \p one_side - All particles on one side of the box
!!\image html par_one_side_thumb.png
!!- \p random - random positions within the box
!!\image html par_random_thumb.png
!!- \p pairs - Particles in pairs randomly placed within the box, can be used to look at two-particle dispersion.
!!\image html par_pairs_thumb.png
!!
!!We can choose between particle type using the particle_count parameter in run.in, options available are:\n
!!
!!- \p fluid - Default option, massless particles \f${d\mathbf{x}_i}/{dt}=u(\mathbf{x}_i,t)\f$, the following plot shows the trajectories of 100 particles in the ABC flow:
!!\image html ABC_stokes_0_thumb.png 
!!- \p inertial - Particles have a mass therefore feel a stokes drag: 
!!\f[\frac{d\mathbf{x}_i}{dt}=\mathbf{u}_i, \qquad \frac{d\mathbf{u}_i}{dt}=\frac{u(\mathbf{x}_i,t)-\mathbf{u}_i}{\tau},\f] where \f$\tau\f$ is the stokes number. Note we define the reciprocal of the stokes number in run.in via the parameter part_stokes. The following plot shows the trajectories of 100 particles in the ABC flow with a stokes number 1:
!!\image html ABC_stokes_1_thumb.png
