module quasip
  !QUASI PARTICLES - AT PRESENT THIS EVOLVES FLUID/INERTIAL/QUASI PARTICLES IN THE CODE
  !IN THE INTEREST OF CALCULATING PARTICLE STATISTICS ONE CAN IGNORE THE SUPERFLUID VELOCITY
  !AND USE ONLY THE NORMAL FLUID 
  use cdata
  use general
  use timestep
  use output
  use tree
  use stiff_solver
  use hamiltonian
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
      !**************************FLUID/INERITAL-PARTICLES*************************
      !put in all the fluid/inertial particle initial conditions in the code here
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
      !****************************QUASI-PARTICLES***********************
      !put in all the quasi particle initial conditions in the code here
      !they all should have the check for particle type quasi
      case('quasi1')
        !set up the stiff ode solver
        select case(initf)
          case('single_line')
          case default
            call warning_message('setup_quasip',&
            'ideally we want initf=single_line for this config')
        end select
        select case(particle_type)
          case('quasi')
            !single particle at one side of the box
            if (quasi_pcount>1) call fatal_error('setup_quasip','only one particle for this initial conditon only')
            g(1)%x(1)=-0.49*box_size
            g(1)%x(2)=-delta/2. ; g(1)%x(3)=0.
            g(1)%p(1)=1.0001*pfermi
            g(1)%p(2)=0. ; g(1)%p(3)=0.
            g(1)%xold(1,:)=g(1)%x ; g(1)%pold(1,:)=g(1)%p
         case default
            call fatal_error('setup_quasip','particle type must be quasi')
      end select
      case('quasi2')
        select case(particle_type)
          case('quasi')
            !multiple particles at one side of the box
            do i=1, quasi_pcount
              g(i)%x(1)=-0.49*box_size
              call random_number(g(i)%x(2))
              call random_number(g(i)%x(3))
              g(i)%x(2)=5.*(1.-2.*g(i)%x(2))*delta ; g(i)%x(3)=5.*(1.-2.*g(i)%x(3))*delta
              g(i)%p(1)=1.0001*pfermi
              g(i)%p(2)=0. ; g(i)%p(3)=0.
              g(i)%xold=0. ; g(i)%pold=0.
              g(i)%xold(1,:)=g(i)%x ; g(i)%pold(1,:)=g(i)%p
            end do
         case default
            call fatal_error('setup_quasip','particle type must be quasi')
      end select
      case default
        call fatal_error('setup_quasip','initg not set to available value')
    end select
    select case(particle_type)
      case('quasi')
        !finally intialise the backwards difference coefficients array
        call set_BDF_coeff !stiff_solver.mod
    end select
  end subroutine
  !************************************************************
  subroutine quasip_evolution
    !PUT TOGETHER ALL THE NORMAL/QUASI PARTICLE ROUTINES TO EVOLVE
    !PARTICLES IN THE CODE AND THEN PERFORM ANALYSIS
    implicit none
    integer :: i
    !timestep the particles
     select case(particle_type)
       case('quasi')
         call timestep_quasip
       case default
         call timestep_fluid
     end select
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
!******************QUASI-PARTICLES*****************************
  subroutine timestep_quasip
    !TIMESTEP THE QUASI-PARTICLES USING VARIOUS SCHEMES
    implicit none
    integer, parameter :: order=4 !order of backwards difference scheme
    integer :: i !used to loop over particles
    do i=1, quasi_pcount
      !move the particle - backwards difference
      call BDF(i,order) !stiff_solver.mod

      !enforce periodicity
      if (periodic_bc) then
        !if a particle leaves one side of the box
        !-------------x------------------     
        if (g(i)%x(1)>(box_size/2.)) then
          select case(initg)
            !for certain initial conditions end the run once the particles have left the box
            case('quasi1','quasi2')
              call fatal_error('quasip_timestep','particles have left the box!')
          end select
          g(i)%x(1)=g(i)%x(1)-box_size
        else if (g(i)%x(1)<(-box_size/2.)) then
          select case(initg)
            !for certain initial conditions end the run once the particles have left the box
            case('quasi1','quasi2')
              call fatal_error('quasip_timestep','particles have left the box!')
          end select
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
  !******************NORMAL PARTICLES*****************************
  subroutine timestep_fluid
    !TIMESTEP THE FLUID/INERTIAL PARTICLE USING SELECTED SCHEME
    implicit none
    real, parameter :: nstokes=0. !stokes number
    real :: u(3)
    integer :: i
    do i=1, quasi_pcount
      call velocity_fluidp(i,u)
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
  subroutine velocity_fluidp(i,u)
    !CALCULATE THE VELOCITY INDUCED BY THE VORTICES
    implicit none 
    integer, intent(IN) :: i
    real, intent(OUT) :: u(3)
    real :: u_sup(3), u_norm(3)
    integer :: peri, perj, perk
    logical :: superfluid_on=.true. !set to false to only use normal fluid
      u_sup=0. !must be zeroed for intially
    if (superfluid_on) then
      select case(velocity)
        case('LIA','BS')
          call biot_savart_general(g(i)%x,u_sup) !timestep.mod
        case('Tree')
          call tree_walk_general(g(i)%x,vtree,(/0.,0.,0./),u_sup)
          if (periodic_bc) then
            !we must shift the mesh in all 3 directions, all 26 permutations needed!
            do peri=-1,1 ; do perj=-1,1 ; do perk=-1,1
              if (peri==0.and.perj==0.and.perk==0) cycle
              call tree_walk_general(g(i)%x,vtree, &
                   (/peri*box_size,perj*box_size,perk*box_size/),u_sup) !tree.mod
            end do ; end do ;end do
          end if
      end select
    end if
    !normal fluid velocity
    call get_normal_velocity(g(i)%x,u_norm) !normal_fluid.mod
    u=u_sup+u_norm
  end subroutine
  !************************************************************
  subroutine diagnostics_quasip
    !DIAGNOSTICS FOR ALL PARTICLES
    implicit none
    real :: energy !quasi_particles
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
    select case(particle_type)
      case('quasi')
        !tidy this up, add qp_energy to the top and write energy into the ts structure
        !may end up with a fluid/inertial particle energy routine
        energy=0.
        do i=1, quasi_pcount
          energy=energy+g(i)%energy
        end do
        open(unit=73,file='data/qp_energy.log',position='append')
          write(73,*) t, energy
        close(73)
    end select
  end subroutine
  !************************************************************
end module
