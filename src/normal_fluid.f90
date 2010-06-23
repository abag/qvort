module normal_fluid
  !NORMAL FLUID COMPONENT IN THE EQUATION OF MOTION
  use cdata
  use general
    !parameters used in get_normal_fluid
    real, parameter, private :: vel_xflow=1.
    real, parameter, private :: abc_A=1., abc_B=1., abc_C=1.
    real, private :: abc_k
    contains 
    !************************************************************
    subroutine setup_normal_fluid
      implicit none
      abc_k=2.*pi/box_size
      write(*,*) 'normal fluid velocity field is: ', trim(normal_velocity)
      select case(normal_velocity)
        case('xflow')
          write(*,'(a,f6.3)') ' u(x)=', vel_xflow
        case('ABC')
          write(*,*) 'A=', abc_A, ' B=', abc_B, ' C=', abc_C
      end select
      write(*,*) 'mutual friction coefficients alpha=',alpha(1),' alpha`=',alpha(2) 
    end subroutine
    !************************************************************
    subroutine get_normal_velocity(x,u)
      implicit none
      real, intent(IN) :: x(3) !position of the particle
      real, intent(OUT) :: u(3) !velocity at x
      select case(normal_velocity)
        case('zero')
          u=0.
        case('xflow')
          u=0. ; u(1)=vel_xflow
        case('ABC')
          u(1)=abc_B*cos(abc_k*x(2))+abc_C*sin(abc_k*x(3))
          u(2)=abc_C*cos(abc_k*x(3))+abc_A*sin(abc_k*x(1))
          u(3)=abc_A*cos(abc_k*x(1))+abc_B*sin(abc_k*x(2))
        case('KS')
          call fatal_error('normal_fluid', 'KS not ready yet')
        case default
          call fatal_error('normal_fluid.mod:get_normal_fluid', &
          'correct parameter for normal_veloctity not set')
      end select
    end subroutine
end module
