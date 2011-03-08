!> all the routines to employ forcing in the code note this is different to the normal velocity for a vortex system
!> normal fluid enters as \f$\mathbf{u}=\alpha \mathbf{s}' \times (\mathbf{u}_\mathrm{n}-\mathbf{u}_\mathrm{s})+
!>\alpha' \mathbf{s}' \times [\mathbf{s}' \times (\mathbf{u}_\mathrm{n}-\mathbf{u}_\mathrm{s})]\f$ whereas forcing is directly
!>added to vortex velocity i.e. \f$\mathbf{u}=\mathbf{u}_\mathrm{force}\f$
module forcing
  use cdata
  !> @param force_direction a helper vector to force a particle in 3 spatial dimensions
  real, private :: force_direction(3)=0. 
  contains
  !>check all the necessary conditions to use forcing are set in run.in
  subroutine setup_forcing()
    implicit none
    select case(force)
      case('off')
        write(*,*) 'no forcing employed'
      case('top_boundary')
        !check initial conditions
         select case(initf)
           case('single_line', 'line_motion')
             write(*,'(a,f6.3,a,f6.3)') 'forcing top boundary with amplitude', force_amp, &
             'frequency', force_freq
           case default
             call fatal_error('forcing.mod:setup_forcing', &
             'incorrect initf for forcing to be applied') !cdata.mod
         end select
      case('box_shake')
        write(*,'(a,f5.3,a,i4.3,a)') 'box shaking forcing with amplitude ', force_amp, &
        ' forcing direction changed every ', ceiling(force_freq), 'timesteps'
      case('delta_corr')
        write(*,'(a,f5.3)') 'delta correlated (time/space) forcing with amplitude ', force_amp
      case default
        call fatal_error('forcing.mod:setup_forcing', &
        'incorrect forcing parameter set in run.in') !cdata.mod
      end select
  end subroutine
  !***********************************************
  !>force the particles - an additional velocity at position
  !>of f(i)%x
  subroutine get_forcing(i,u)
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
      case('box_shake')
        if (i==1) then
          !do we need to generate a new forcing direction?
          if (mod(itime,ceiling(force_freq))==0) then
            call random_number(force_direction)
            force_direction=force_direction*2.-1.
            !normalise
            force_direction=force_direction/sqrt(dot_product(force_direction,force_direction))
          end if
        end if
        u=force_direction*force_amp
      case('delta_corr')
        !at each time and position generate a new random vector
        call random_number(force_direction)
        force_direction=force_direction*2.-1.
        !normalise it
        force_direction=force_direction/sqrt(dot_product(force_direction,force_direction))
        !use this as the forcing
        u=force_direction*force_amp
    end select
  end subroutine
  !***********************************************
  !>a general forcing routine which accepts a position
  subroutine get_forcing_gen(x,u)
    implicit none
    real, intent(IN) :: x(3)
    real, intent(OUT) :: u(3)
    select case(force)
      case('off')
        u=0.
      case('top_boundary')
        if ((x(3)-box_size/2.)>-1.5*delta) then
          !particle is sufficiently close to top boundary to force
          u=0.
          !sinusoidal forcing in the x direction
          u(1)=delta*force_amp*sin(force_freq*t/(2*pi))
        end if
      case('delta_corr')
        !at each time and position generate a new random vector
        call random_number(force_direction)
        force_direction=force_direction*2.-1.
        !normalise it
        force_direction=force_direction/sqrt(dot_product(force_direction,force_direction))
        !use this as the forcing
        u=force_direction*force_amp
      case default
        call fatal_error('forcing.mod','this forcing function is not compatable with a general call')
    end select
  end subroutine
end module
