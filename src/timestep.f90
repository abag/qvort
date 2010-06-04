module timestep
  use cdata
  use general
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
    implicit none
    integer, intent(IN) :: i
    real :: u(3)
    real :: f_dot(3), f_ddot(3) !first and second derivs
    !what scheme are we using? (LIA/BS)
    select case(velocity)
      case('LIA')
        !use the local induction approximation
        call get_deriv_1(i,f_dot) !general.mod
        call get_deriv_2(i,f_ddot) !general.mod
        u=5E-4*cross_product(f_dot,f_ddot)
      case('BS')
        !full (nonlocal) biot savart
        print*, 'not yet ready'
        call fatal_error !cdata.mod
      case default
        print*, 'correct value for velocity in run.in has not been set'
        print*, 'options are: LIA, BS'
        call fatal_error !cdata.mod
    end select
  end subroutine
end module
