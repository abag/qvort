!> all the routines to employ forcing in the code note this is different to the normal velocity for a vortex system
!> normal fluid enters as \f$\mathbf{u}=\alpha \mathbf{s}' \times (\mathbf{u}_\mathrm{n}-\mathbf{u}_\mathrm{s})+
!>\alpha' \mathbf{s}' \times [\mathbf{s}' \times (\mathbf{u}_\mathrm{n}-\mathbf{u}_\mathrm{s})]\f$ whereas forcing is directly
!>added to vortex velocity i.e.
!>\f$\mathbf{u}=\mathbf{u}_\textrm{vortex}+\mathbf{u}_\mathrm{force}\f$.
!>See \ref FORCE for more details.
!>\todo forcing shutoff after certain time
!>\todo forcing frequency in LS_force
module forcing
  use cdata
  use general
  !> @param force_direction a helper vector to force a particle in 3 spatial dimensions
  real, private :: force_direction(3)=0. 
  !> @param LS_k Large scale forcing wavenumber
  !> @param LS_A Large scale forcing - vector perp to k
  !> @param LS_B Large scale forcing - vector perp to k
  real, private :: LS_k(3), LS_A(3), LS_B(3)
  !> @param LS_k2 Large scale forcing squared wavenumber
  real :: LS_k2
  contains
  !>check all the necessary conditions to use forcing are set in run.in
  subroutine setup_forcing()
    implicit none
    select case(force)
      case('off')
        write(*,*) 'no forcing employed'
      case('xflow')
        write(*,*) 'adding an imposed flow in the x direction'
        write(*,'(a,f6.3)') ' u(x)=', force_amp
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
      case('LS_force')
        write(*,'(a,f5.3)') 'large scale forcing (similar to KS) with amplitude ', force_amp
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
      case('xflow')
        u=0.
        u(1)=force_amp
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
      case('LS_force')
        !this has it's own subroutine
        if (i==1) then
          !only do this once per time-step
          call get_LS_kAB!forcing.mod
          if (mod(itime,shots)==0) then
            !print to file
            open(unit=47,file='./data/LS_forcing.log',position='append')
              write(47,*) LS_k, LS_A, LS_B
            close(47)
          end if 
        end if
        u=cross_product(LS_A,LS_k)*sin(dot_product(LS_k,f(i)%x))/LS_k2+cross_product(LS_B,LS_k)*cos(dot_product(LS_k,f(i)%x))/LS_k2
        !use this as the forcing-multiply by amplitude
        u=u*force_amp
    end select
  end subroutine
  !***********************************************
  !>Large scale forcing routine, very similar to KS with one
  !>random large scale mode, this sets up k, A and B for use in
  !>main forcing routine
  subroutine get_LS_kAB
    real,dimension(3) ::  newa, newb, unitk !helper 
    integer :: i
    !create random vectors
    call random_number(newa) ; call random_number(newb)
    newa=2.*newa-1. ; newa2=2.*newa2-1.
    newa=newa/(sqrt(sum(newa**2))) ; newb=newb/(sqrt(sum(newb**2)))
    !generate random wavevector
    LS_k=.0
    do while (sum(LS_k)<epsilon(0.))
      do i=1,3
        call random_number(LS_k(i))
        if (LS_k(i)<0.33333) then
          LS_k(i)=1.
        else if (LS_k(i)<0.666) then
          LS_k(i)=-1
        else
          LS_k(i)=0.
        end if
      end do
    end do
    LS_k=LS_k*2*pi/box_size !scale by 2\pi for periodicity - box size must enter here
    LS_k2=sqrt(sum(LS_k**2))
    unitk=LS_k/LS_k2 !unit vector
    !create A and B by taking cross product of newa/b with unit k
    LS_A=cross_product(newa,unitk) ; LS_B=cross_product(newb,unitk)
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
        call fatal_error('forcing.mod','this forcing function is not compatible with a general call')
    end select
  end subroutine
end module
!> \page FORCE Forcing
!!Forcing of filaments is set in run.in throught the parameter
!!force this accepts the following arguements\n
!!
!!\p zero - no forcing, this is the default \n
!!\p xflow - imposed flow in the x direction, velocity set by force_amp \n
!!\p box_shake - correlated forcing moving box randomly with frequency set by
!!\p force_feq in run.in \n
!!\p delta_corr - delta correlated forcing in time and space with amplitude set
!!by \p force_amp parameter \n
!!\p top_boundary - sinusoidal forcing in the x direction at the top of the box
!!with a frequency and amplitude set in run.in (force_freq, force_amp)\n
!!\p LS_force - please document me\n

