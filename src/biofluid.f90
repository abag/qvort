module biofluid
  use cdata
  use general
  use timestep
  use output
  use normal_fluid
  real, parameter :: bio_k(3)=(/0,0,1/) !gravity
  contains
  !************************************************************
  subroutine bio_swimming(i)
    implicit none
    integer, intent(IN) :: i
    real :: pdot(3), u(3), omega(3)
    !update the particles swimming direction
    !get the velocity field at p(i)%x
    call get_normal_velocity(p(i)%x,u)
    call get_normal_vorticity(p(i)%x,omega)
    pdot=(1./(2*bio_G))*(bio_k-dot_product(bio_k,p(i)%p)*p(i)%p)+&
         0.5*cross_product(omega,p(i)%p)
    if (maxval(abs(p(i)%pdot1))==0) then
      p(i)%p(:)=p(i)%p(:)+dt*pdot !euler
    else if (maxval(abs(p(i)%pdot2))==0) then
      !first order adams-bashforth
      p(i)%p(:)=p(i)%p(:)+three_twos*dt*pdot-one_half*dt*p(i)%pdot1(:)
    else
      !second order adams-bashforth
      p(i)%p(:)=p(i)%p(:)+twenty_three_twelve*dt*pdot-four_thirds*dt*p(i)%pdot1(:)+five_twelths*dt*p(i)%pdot2(:)
    end if
    !store old pdots
    p(i)%pdot2(:)=p(i)%pdot1(:)
    p(i)%pdot1(:)=pdot
    !now get the particle velocity
    p(i)%u=bio_V*p(i)%p+u
    !print*, i, omega, p(i)%p
  end subroutine
end module 
