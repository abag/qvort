!>The particles aspect of the code. At present this evolves either fluid or
!>inertial particles in the code. The type of particle is set by the
!>particle_type parameter in run.in Fluid/Inertial particles use the 
!>normal fluid velocity. The initial setup of the particles is set throught the initp parameter.
module particles
  use cdata
  use general
  use timestep
  use output
  real :: part_maxu=0., part_maxdu=0., part_urms=0. !velocity information
  real :: part_sep=0. !average particle separation
  contains
  !************************************************************
  !>setup the particles in as set by initp
  subroutine setup_particles
    implicit none
    integer :: i
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
    do i=1, part_count
      call velocity_fluidp(i,u)
      select case(particle_type)
        case('fluid')
          p(i)%u=u
        case('inertial')
          !an euler step to get the velocity
          p(i)%u=p(i)%u1+dt*(u-p(i)%u1)*part_stokes
          !adjust the velocity due to stokes drag
      end select
      !euler step the particles
      p(i)%x(:)=p(i)%x(:)+dt*p(i)%u
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
  end subroutine
  !************************************************************
  !> calculate the velocity at each particle
  !> note we can set this in here to be the velocity induced by
  !> the superfluid vortices but this is unphysical
  !> really this should be false by default - consider putting in run.in
  subroutine velocity_fluidp(i,u)
    implicit none 
    integer, intent(IN) :: i
    real, intent(OUT) :: u(3)
    real :: u_sup(3), u_norm(3)
    integer :: peri, perj, perk
    logical :: superfluid_on=.false. !set to false to only use normal fluid
    u_sup=0. !must be zeroed for intially
    if (superfluid_on) then
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
        do i=1, quasi_pcount, 2
          part_sep=part_sep+sqrt((p(i)%x(1)-p(i+1)%x(1))**2+&
                                 (p(i)%x(2)-p(i+1)%x(2))**2+&
                                 (p(i)%x(3)-p(1+1)%x(3))**2)
        end do
        part_sep=part_sep/(part_count/2)
    end select
    open(unit=78,file='data/par_ts.log',position='append')
    if (itime==shots) then
      write(78,*) '%-var--t-----maxu---maxdu----urms---part_sep--'
      if (particles_only) write(*,*) '%-var--t-----maxu---maxdu----urms---part_sep--'
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
