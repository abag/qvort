module timestep
  !ANYTHING RELATED TO VELOCITY SHOULD ENTER HERE
  use cdata
  use general
  use normal_fluid
  use forcing
  contains
  !*************************************************
  subroutine pmotion()
    !implement adams bashforth time-stepping scheme to move particles
    implicit none
    real :: u(3) !dummy variable used to store velocities
    integer :: i
    do i=1, pcount
      if (f(i)%infront==0) cycle !check for 'empty' particles
      call calc_velocity(u,i)
      f(i)%u(:)=u(:) !store the velocity for time-step
    end do
    do i=1, pcount
      if (f(i)%infront==0) cycle !check for 'empty' particles
      if (maxval(abs(f(i)%u1))==0) then
        f(i)%x(:)=f(i)%x(:)+dt*f(i)%u(:) !euler
      else if (maxval(abs(f(i)%u2))==0) then
        !first order adams-bashforth
        f(i)%x(:)=f(i)%x(:)+three_twos*dt*f(i)%u(:)-one_half*dt*f(i)%u1(:)
      else
        !second order adams-bashforth
        f(i)%x(:)=f(i)%x(:)+twenty_three_twelve*dt*f(i)%u(:)-four_thirds*dt*f(i)%u1(:)+five_twelths*dt*f(i)%u2(:)
      end if
      f(i)%u2(:)=f(i)%u1(:) !store our old velocities 
      f(i)%u1(:)=f(i)%u(:)
    end do
    t=t+dt !time needs incrementing 
  end subroutine
  !*************************************************
  subroutine calc_velocity(u,i)
    !get the velocity of each particle subject to the superfluid velocity
    !plus any normal fluid/forcing
    implicit none
    integer, intent(IN) :: i
    real :: u(3), u_norm(3), u_force(3)
    real :: curv, beta
    real :: f_dot(3), f_ddot(3) !first and second derivs
    !what scheme are we using? (LIA/BS)
    select case(velocity)
      case('LIA')
        !use the local induction approximation
        call get_deriv_1(i,f_dot) !general.mod
        call get_deriv_2(i,f_ddot) !general.mod
        !calculate the curvature
        curv=sqrt(dot_product(f_ddot,f_ddot))
        if (curv<epsilon(0.)) then
          curv=epsilon(0.) !we must check for zero curvature
        end if
        !caluculate beta based on the curvature
        beta=(9.97E-4/(4.*pi))*log(1E8/curv)
        u=beta*cross_product(f_dot,f_ddot) !general.mod
      case('BS')
        !full (nonlocal) biot savart
        print*, 'not yet ready'
        call fatal_error('timestep.mod:calc_velocity', &
                         'BS not ready yet!') !cdata.mod
      case default
        print*, 'correct value for velocity in run.in has not been set'
        print*, 'options are: LIA, BS'
        call fatal_error('timestep.mod:calc_velocity', & 
        'correct value for "velocity" in run.in has not been set') !cdata.mod
    end select
    !now account for mutual friction - test if alpha's are 0
    if (sum(alpha)>epsilon(0.)) then
      call get_normal_velocity(f(i)%x,u_norm) !normal_fluid.mod
      u=u+alpha(1)*cross_product(f_dot,(u_norm-u))- &
          alpha(2)*cross_product(f_dot,cross_product(f_dot,(u_norm-u)))
    end if
    !forcing?
    call get_forcing(i,u_force)
    u=u+u_force
  end subroutine
  !*************************************************
  subroutine mesh_velocity
    !get the velocity at each point on the mesh for spectra etc.
    implicit none
    integer :: i,j,k
    do k=1, mesh_size
      do j=1, mesh_size
        do i=1, mesh_size
          !superfluid velocity
          !select case(velocity)
          !  case('BS')
            !operations would go in here
          !end select
          !normal fluid
          call get_normal_velocity(mesh(k,j,i)%x,mesh(k,j,i)%u_norm)
        end do
      end do
    end do
  end subroutine
end module
