module quasip
  !QUASI PARTICLES - AT PRESENT THIS EVOLVES FLUID/INERTIAL PARTICLES IN THE CODE
  !IN THE INTEREST OF CALCULATING PARTICLE STATISTICS ONE CAN IGNORE THE SUPERFLUID VELOCITY
  !AND USE ONLY THE NORMAL FLUID 
  use cdata
  use general
  use timestep
  use output
  use tree
  real :: qp_maxu=0., qp_maxdu=0., qp_urms=0. !velocity information
  real :: qp_sep=0. !avergae particle separation
  !*******CONSTANTS ASSOCIATED WITH QUASI-PARTICLES****************
  real, private, parameter :: mbare=5.017e-24
  real, private, parameter :: mfermi=1.51e-23 
  real, private, parameter :: pfermi=8.281e-20
  real, private, parameter :: efermi=2.27e-16
  real, private, parameter :: vfermi=5.48e3
  real, private, parameter :: delta_gap=2.43e-19
  real, private, parameter :: unitdt=2.9e-6 !unit time for matlab  ode solver 
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
            end do
         case default
            call fatal_error('setup_quasip','particle type must be quasi')
      end select
      case default
        call fatal_error('setup_quasip','initg not set to available value')
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
    !TIMESTEP THE QUASI-PARTICLES USING EULER SCHEME
    implicit none
    real :: rdot(3), pdot(3)
    integer :: i
    do i=1, quasi_pcount
      call velocity_quasip(i,rdot,pdot)
      !euler step both the position and momentum
      g(i)%x(:)=g(i)%x(:)+dt*rdot(:)
      g(i)%p(:)=g(i)%p(:)+dt*pdot(:)
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
  subroutine velocity_quasip(i,rdot,pdot)
  !for each particle we 
  ! 1. define the hamiltonian
       !in order to do this we require the velocity (vortex) derivatives in all 3 dimension (7 evaluations)
  ! 2. intially update the particles using an euler step
  ! 3. build up to solving the stiff problem 
    implicit none
    integer, intent(IN) :: i
    real, intent(OUT) :: rdot(3), pdot(3)
    real :: u_sup(3), du_sup_x(3), du_sup_y(3), du_sup_z(3) !superfluid velocity and derivatives
    real :: u_plusx(3), u_minusx(3) !superfluid velocity at plus minus x
    real :: u_plusy(3), u_minusy(3) !superfluid velocity at plus minus y
    real :: u_plusz(3), u_minusz(3) !superfluid velocity at plus minus z
    real :: epsilonp !used in the hamiltonian
    real :: dist, min_dist !need the minimum distance between the particle and the vortex 
    integer :: peri, perj, perk !used to loop in periodic cases
    integer :: j !used for loops
    !Get the superfluid velocity at the position of the particle
    !account for all possible velocity fields
    u_sup=0. !must be zeroed for intially
    u_plusx=0. ; u_minusx=0. ; u_plusy=0. ; u_minusy=0.
    u_plusz=0. ; u_minusz=0.
    !first we need to find the minimum distance between the particle and the vortice's
    select case(velocity)
      case('LIA','BS')
        min_dist=100. !arbitrarily large
        do j=1,pcount !loop over all vortex points
          dist=dist_gen(g(i)%x,f(j)%x)
          if (dist<min_dist) min_dist=dist
        end do
      case('Tree')
        !we need a tree minimum distance algorithm
        call fatal_error('velocity_quasip','tree adaptation not ready yet')
    end select
    !now we check wether the particle is within the vortex core or outside
    if (min_dist<corea) then
      !empty core so 0 velocity
      u_sup=0. ; u_plusx=0. ; u_minusx=0.
      u_plusy=0. ; u_minusy=0.
      u_plusz=0. ; u_minusz=0.
    else
      select case(velocity)
        case('LIA','BS')
          call biot_savart_general(g(i)%x,u_sup) !timestep.mod
          call biot_savart_general(g(i)%x+(/delta,0.,0./),u_plusx) !timestep.mod
          call biot_savart_general(g(i)%x+(/-delta,0.,0./),u_minusx) !timestep.mod
          call biot_savart_general(g(i)%x+(/0.,delta,0./),u_plusy) !timestep.mod
          call biot_savart_general(g(i)%x+(/0.,-delta,0./),u_minusy) !timestep.mod
          call biot_savart_general(g(i)%x+(/0.,0.,delta/),u_plusz) !timestep.mod 
          call biot_savart_general(g(i)%x+(/0.,0.,-delta/),u_minusz) !timestep.mod
          if (periodic_bc) then
            !we must shift the mesh in all 3 directions, all 26 permutations needed!
            !needs adpating for periodic boundary conditons - see above
            do peri=-1,1 ; do perj=-1,1 ; do perk=-1,1
              if (peri==0.and.perj==0.and.perk==0) cycle
              call biot_savart_general_shift(g(i)%x,u_sup, &
                   (/peri*box_size,perj*box_size,perk*box_size/))!timestep.mod
            end do ; end do ;end do
          end if
        case('Tree')
          call fatal_error('velocity_quasip','tree adaptation not ready yet')
      end select
    end if
    !now we must find both rdot and pdot using the hamiltonian for the quasi-particle
    epsilonp=(dot_product(g(i)%p,g(i)%p)/(2.*mfermi))-efermi
    rdot(:)=(epsilonp/sqrt(epsilonp**2+delta_gap**2))*g(i)%p(:)/mfermi+u_sup
    !let us now calculate the derivative of the superfluid velocity field
    du_sup_x(:)=((u_plusx-u_minusx)/(2.*delta))
    du_sup_y(:)=((u_plusy-u_minusy)/(2.*delta))
    du_sup_z(:)=((u_plusz-u_minusz)/(2.*delta))
    !calculate pdot
    pdot(1)=-dot_product(g(i)%p,du_sup_x)
    pdot(2)=-dot_product(g(i)%p,du_sup_y)
    pdot(3)=-dot_product(g(i)%p,du_sup_z)
    !finally store the energy for output
    g(i)%energy=sqrt(epsilonp**2+delta_gap**2)+dot_product(g(i)%p,u_sup)
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
