module quasip
  use cdata
  use general
  use timestep
  use output
  real :: qp_maxu=0., qp_maxdu=0., qp_urms=0. !velocity information
  real :: qp_sep=0. !avergae particle separation
  contains
  !************************************************************
  subroutine setup_quasip
    !SET UP THE QUASI PARTICLES INTIAL POSITION (VELOCITY)
    implicit none
    integer :: i
    allocate(g(quasi_pcount))
    write(*,*) 'setting up particles in ', trim(initg),' configuration'
    write(*,'(a,a,a,i4.1)') ' particle type: ', trim(particle_type), ', particle count: ', quasi_pcount
    select case(initg)
      case('one-side')
        !particles all on one side of the box (x-axis) 
        do i=1, quasi_pcount
          if (box_size<epsilon(0.))  then
            call fatal_error('setup_quasip','particles need periodic bcs')
          end if
          g(i)%x(1)=-box_size/2.
          call random_number(g(i)%x(2)) ; call random_number(g(i)%x(3))
          g(i)%x(2)=box_size*g(i)%x(2)-box_size/2.
          g(i)%x(3)=box_size*g(i)%x(3)-box_size/2.
          g(i)%u=0. ; g(i)%u1=0. ; g(i)%u2=0. !0 the velocity fields
          g(i)%u1(1)=1. ; g(i)%u2(1)=1. !now give an intial (x) velocity
        end do
      case('random')
        !particles in random positions
        do i=1, quasi_pcount
          call random_number(g(i)%x(1))
          call random_number(g(i)%x(2)) ; call random_number(g(i)%x(3))
          if (box_size>0.) then
            g(i)%x(1)=box_size*g(i)%x(1)-box_size/2.
            g(i)%x(2)=box_size*g(i)%x(2)-box_size/2.
            g(i)%x(3)=box_size*g(i)%x(3)-box_size/2.
          end if
          g(i)%u=0. ; g(i)%u1=0. ; g(i)%u2=0. !0 the velocity fields
        end do
      case('pairs')
        !particles in pairs in a random position
        !check that the number of particles is a multiple of 2
        if (mod(quasi_pcount,2)/=0) then
          call fatal_error('setup_quasip','quasi_pcount must be a multiple of 2')
        end if
        do i=1,  quasi_pcount, 2
          call random_number(g(i)%x(1))
          call random_number(g(i)%x(2)) ; call random_number(g(i)%x(3))
          if (box_size>0.) then
            g(i)%x(1)=box_size*g(i)%x(1)-box_size/2.
            g(i)%x(2)=box_size*g(i)%x(2)-box_size/2.
            g(i)%x(3)=box_size*g(i)%x(3)-box_size/2.
          end if
          g(i)%u=0. ; g(i)%u1=0. ; g(i)%u2=0. !0 the velocity fields
          !finally sort out the position of its neighbour
          g(i+1)%x=g(i)%x ; g(i+1)%x(3)=g(i)%x(3)+0.001
          g(i+1)%u=0. ; g(i+1)%u1=0. ; g(i+1)%u2=0. !0 the velocity fields
        end do
      case default
        call fatal_error('setup_quasip','initg not set to available value')
    end select
    select case(particle_type)
      case('quasi')
        call fatal_error('setup_qausi', 'quasi particles not ready yet')
        !Anything extra that the quasi particles setup needs can go in here
    end select
  end subroutine
  !************************************************************
  subroutine quasip_evolution
    !PUT TOGETHER ALL THE QUASI PARTICLE ROUTINES TO EVOLVE
    !PARTICLES IN THE CODE AND THEN PERFORM ANALYSIS
    implicit none
    integer :: i
    !timestep the particles
    call timestep_quasip
    !particle diagnostics
    if (mod(itime,shots)==0) then
      call diagnostics_quasip
    end if
    !print the particles to file
    if (mod(itime,shots)==0) then
      call printg(itime/shots) !output.mod
    end if
    !any other business in here
  end subroutine
  !************************************************************
  subroutine timestep_quasip
    !TIMESTEP THE QUASI PARTICLE I USING SELECTED SCHEME
    implicit none
    real, parameter :: nstokes=0. !stokes number
    real :: u(3)
    integer :: i
    do i=1, quasi_pcount

      call velocity_quasip(i,u)
      select case(particle_type)
        case('fluid')
          g(i)%u=u
        case('inertial')
          !an euler step to get the velocity
          g(i)%u=g(i)%u1+dt*(u-g(i)%u1)*nstokes
          !adjust the velocity due to stokes drag
      end select
      !euler step the particles
      g(i)%x(:)=g(i)%x(:)+dt*g(i)%u
      !store the old velocities
      g(i)%u2=g(i)%u1
      g(i)%u1=g(i)%u
      !enforce periodicity
      if (periodic_bc) then
        !if a particle leaves one side of the box
        !-------------x------------------     
        if (g(i)%x(1)>(box_size/2.)) then
          g(i)%x(1)=g(i)%x(1)-box_size
        else if (g(i)%x(1)<(-box_size/2.)) then
          g(i)%x(1)=g(i)%x(1)+box_size
        end if
        !-------------y------------------
        if (g(i)%x(2)>(box_size/2.)) then
          g(i)%x(2)=g(i)%x(2)-box_size
        else if (g(i)%x(2)<(-box_size/2.)) then
          g(i)%x(2)=g(i)%x(2)+box_size
        end if
        !-------------z------------------
        if (g(i)%x(3)>(box_size/2.)) then
          g(i)%x(3)=g(i)%x(3)-box_size
        else if (g(i)%x(3)<(-box_size/2.)) then
          g(i)%x(3)=g(i)%x(3)+box_size
        end if
        !--------------------------------
        !put it back in the other side....
      end if
    end do
  end subroutine
  !************************************************************
  subroutine velocity_quasip(i,u)
    !CALCULATE THE VELOCITY INDUCED BY THE VORTICES
    implicit none 
    integer, intent(IN) :: i
    real, intent(OUT) :: u(3)
    real :: u_sup(3), u_norm(3)
    !get the superfluid velocity
    call biot_savart_general(g(i)%x,u_sup) !timestep.mod
    !we can put a tree biot-savart in here eventually
    !normal fluid velocity
    call get_normal_velocity(g(i)%x,u_norm) !normal_fluid.mod
    u=u_sup+u_norm
  end subroutine
  !************************************************************
  subroutine diagnostics_quasip
    implicit none 
    real :: uinfo(quasi_pcount,2)
    integer :: i
    uinfo(:,1)=sqrt(g(:)%u(1)**2+g(:)%u(2)**2+g(:)%u(3)**2)
    uinfo(:,2)=sqrt((g(:)%u1(1)-g(:)%u2(1))**2+&
                    (g(:)%u1(2)-g(:)%u2(2))**2+&
                    (g(:)%u1(3)-g(:)%u2(3))**2)
    maxu=maxval(uinfo(:,1)) ; maxdu=maxval(uinfo(:,2))
    !in here we need to determine
    qp_maxu=maxval(uinfo(:,1)) ; qp_maxdu=maxval(uinfo(:,2))
    qp_urms=sqrt(sum(uinfo(:,1)**2)/quasi_pcount)
    !2 particle separation
    select case(initg)
      case('pairs')
        qp_sep=0.
        do i=1, quasi_pcount, 2
          qp_sep=qp_sep+sqrt((g(i)%x(1)-g(i+1)%x(1))**2+&
                             (g(i)%x(2)-g(i+1)%x(2))**2+&
                             (g(i)%x(3)-g(1+1)%x(3))**2)
        end do
        qp_sep=qp_sep/(quasi_pcount/2)
    end select
    open(unit=78,file='data/par_ts.log',position='append')
    if (itime==shots) then
      write(78,*) '%-var--t-----maxu---maxdu----urms---qp_sep--'
    end if
    write(78,'(i5.3,f6.2,f8.5,f8.5,f8.5,f8.5)') &
itime/shots,t,qp_maxu,qp_maxdu,qp_urms,qp_sep
    close(78)
  end subroutine
  !************************************************************
end module
