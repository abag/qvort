module forcing
  !THIS MODULE CONTAINS ALL THE ROUTINES TO EMPLOY FORCING IN THE CODE
  use cdata
  use general
  contains
  subroutine setup_forcing()
    !essentially check all the necessary conditions to use forcing are set in run.in
    implicit none
    select case(force)
      case('off')
        write(*,*) 'no forcing employed'
      case('top_boundary')
        !check initial conditions
         select case(initf)
           case('single_line', 'line_motion')
             write(*,*) 'forcing top boundary with amplitude', force_amp, &
             'frequency', force_freq
           case default
             call fatal_error('cdata.mod:init_setup', &
             'incorrect initf for forcing to be applied') !cdata.mod
         end select
      case default
        call fatal_error('cdata.mod:init_setup', &
        'incorrect forcing parameter set in run.in') !cdata.mod
      end select
  end subroutine
  !***********************************************
  subroutine get_forcing(i,u)
    !force the particles - an additional velocity
    implicit none
    integer, intent(IN) :: i
    real, intent(OUT) :: u(3)
    select case(force)
      case('off')
        u=0.
      case('top_boundary')
        if ((f(i)%x(3)-box_size/2.)>-1.5*delta) then
          !particle is sufficiently close to top boundary to force
          u=0.
          !sinusoidal forcing in the x direction
          u(1)=delta*force_amp*sin(force_freq*t/(2*pi))
        end if
    end select
  end subroutine
end module
