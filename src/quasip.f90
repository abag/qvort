module quasip
  use Cdata
  use General
  integer :: quasi_pcount
  type quasi !quasi particle structure
    real :: x(3) !position
    real :: u(3), u1(3), u2(3) !stored velocities (adam bash)
  end type
  type(quasi), allocatable :: g(:) !vector
  contains
  !************************************************************
  subroutine initiate_quasip
    !SET UP THE QUASI PARTICLES INTIAL POSITION (VELOCITY)
    implicit none
    integer :: i
    quasi_pcount=100
    allocate(g(quasi_pcount))
    !particles all on one side of the box (x-axis) 
    do i=1, quasi_pcount
      g(i)%x(1)=-box_size/2.
      call random_number(g(i)%x(2)) ; call random_number(g(i)%x(3))
      g(i)%x(2)=box_size*g(i)%x(2)-box_size/2.
      g(i)%x(3)=box_size*g(i)%x(3)-box_size/2.
      g(i)%u=0. ; g(i)%u1=0. ; g(i)%u2=0. !0 the velocity fields
      g(i)%u1(1)=1. ; g(i)%u2(1)=1. !now give an intial (x) velocity
    end do
  end subroutine
  !************************************************************
  subroutine timestep_quasip
    !TIMESTEP THE QUASI PARTICLES USING ADAMS-BASH
    implicit none
  end subroutine
  !************************************************************
  subroutine velocity_quasip
    !CALCULATE THE VELOCITY INDUCED BY THE VORTICES
    implicit none 
  end subroutine
  !************************************************************
  subroutine diagnostics_quasip
    implicit none 
  end subroutine
  !************************************************************
end module
