!>The quasi particles aspect of the code. The superfluid
!>velocity if used to evolve the particle using the stiff_solver and hamiltonian
!>module. The initial setup of the particles is set throught the initg parameter.
module quasip
  use cdata
  use general
  use timestep
  use output
  use tree
  use stiff_solver
  use hamiltonian
  real :: qp_max_u=0., qp_max_pdot=0., qp_urms=0. !velocity/momenta information
  contains
  !************************************************************
  !>setup the quasi particles initial positions and momenta
  subroutine setup_quasip
    implicit none
    integer :: i
    allocate(g(quasi_pcount))
    write(*,*) 'setting up qausi particles in ', trim(initg),' configuration'
    select case(initg)
      !****************************QUASI-PARTICLES***********************
      !put in all the quasi particle initial conditions in the code here
      case('quasi1')
        !set up the stiff ode solver
        select case(initf)
          case('single_line')
          case default
            call warning_message('setup_quasip',&
            'ideally we want initf=single_line for this config')
        end select
        !single particle at one side of the box
        if (quasi_pcount>1) call fatal_error('setup_quasip','only one particle for this initial conditon only')
        g(1)%x(1)=-0.49*box_size
        g(1)%x(2)=-delta/2. ; g(1)%x(3)=0.
        g(1)%p(1)=1.0001*pfermi
        g(1)%p(2)=0. ; g(1)%p(3)=0.
        g(1)%xold=0. ; g(1)%pold=0.
        g(1)%xold(1,:)=g(1)%x ; g(1)%pold(1,:)=g(1)%p
      case('quasi2')
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
        call fatal_error('setup_quasip','initg not set to available value')
    end select
    !finally intialise the backwards difference coefficients array
    call set_BDF_coeff !stiff_solver.mod
  end subroutine
  !************************************************************
  !>evolve the quasi particles and print to file
  subroutine quasip_evolution
    implicit none
    integer :: i
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
  !******************QUASI-PARTICLES*****************************
  !>timestep the quasi particles by calling the backwards difference scheme in 
  !>stiff_solver module
  subroutine timestep_quasip
    implicit none
    integer, parameter :: order=6 !order of backwards difference scheme
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
  !************************************************************
  !>quasi particle diagnostics, maxium velocity and rate of change of
  !>momentum (rdot, pdot)
  subroutine diagnostics_quasip
    implicit none
    real :: energy !quasi_particles energy
    real :: qpinfo(quasi_pcount,2) ! a helper array
    integer :: i
    qpinfo(:,1)=sqrt(g(:)%rdot(1)**2+g(:)%rdot(2)**2+g(:)%rdot(3)**2)
    qpinfo(:,2)=sqrt(g(:)%pdot(1)**2+g(:)%pdot(2)**2+g(:)%pdot(3)**2)
    qp_max_u=maxval(qpinfo(:,1)) ; qp_max_pdot=maxval(qpinfo(:,2))
    qp_urms=sqrt(sum(qpinfo(:,1)**2)/quasi_pcount)
    open(unit=78,file='data/qp_ts.log',position='append')
    if (itime==shots) then
      write(78,*) '%-var-----t-----------max_u---------max_p--------u_rms-----'
    end if
    write(78,'(i5.3,e14.7,e14.7,e14.7,e14.7,e14.7)') &
    itime/shots,t,qp_max_u,qp_max_pdot,qp_urms
    close(78)
    energy=0.
    do i=1, quasi_pcount
      energy=energy+g(i)%energy
    end do
    open(unit=73,file='data/qp_energy.log',position='append')
      write(73,*) t, energy
    close(73)
  end subroutine
  !************************************************************
end module
