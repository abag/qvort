!>the routines contained within this module will inject loops into the code every nth timestep
!>(n is selected in run.in via inject_skip), the size of the loop is set using inject_size and the 
!>topology of the injected loop is set via loop_type, again in run.in.
module inject
  use cdata
  use general
  use periodic
  contains
  subroutine vortex_inject
    implicit none    
    type(qvort), allocatable, dimension(:) :: tmp
    real :: radius !used for loops
    real :: rand1, rand2 !random numbers
    integer :: old_pcount
    integer :: i, j
    !check that inject_size is not 0
    if (inject_size==0) return !leave routine if it is
    !check if inject is set to off
    select case(inject_type)
      case('off')
        return !leave routine
    end select
    !first store the old point count
    old_pcount=pcount
    !now we must resize the array 
    !increment by inject_size using a dummy array tmp
    allocate(tmp(size(f)+inject_size)) ; tmp(:)%infront=0 !0 the infront array
    !copy accross information
    tmp(1:size(f)) = f
    !now deallocate tmp and transfer particles back to f
    !move_alloc is new intrinsic function (fortran 2003)
    call move_alloc(from=tmp,to=f)
    pcount=pcount+inject_size !increase the particle count
    !---------------now check what the inject_type is--------------
    select case(inject_type)
      case('xy-loop')
        radius=(0.75*inject_size*delta)/(2*pi) !75% of potential size
        !loop over particles setting spatial and 'loop' position
        do i=old_pcount+1, pcount
          f(i)%x(1)=radius*sin(pi*real(2*i-1)/inject_size)
          f(i)%x(2)=radius*cos(pi*real(2*i-1)/inject_size)
          f(i)%x(3)=-box_size/2.
          if (i==old_pcount+1) then
            f(i)%behind=pcount ; f(i)%infront=i+1
          else if (i==pcount) then 
            f(i)%behind=i-1 ; f(i)%infront=old_pcount+1
          else
            f(i)%behind=i-1 ; f(i)%infront=i+1
          end if
          !zero the stored velocities
          f(i)%u1=0. ; f(i)%u2=0.
        end do
      case('yz-loop')
        radius=(0.75*inject_size*delta)/(2*pi) !75% of potential size
        !loop over particles setting spatial and 'loop' position
        do i=old_pcount+1, pcount
          f(i)%x(1)=-box_size/2.
          f(i)%x(2)=radius*cos(pi*real(2*i-1)/inject_size)
          f(i)%x(3)=radius*sin(pi*real(2*i-1)/inject_size)
          if (i==old_pcount+1) then
            f(i)%behind=pcount ; f(i)%infront=i+1
          else if (i==pcount) then 
            f(i)%behind=i-1 ; f(i)%infront=old_pcount+1
          else
            f(i)%behind=i-1 ; f(i)%infront=i+1
          end if
          !zero the stored velocities
          f(i)%u1=0. ; f(i)%u2=0.
        end do
      case('rand-yz-loop')
        radius=(0.75*inject_size*delta)/(2*pi) !75% of potential size
        !loop over particles setting spatial and 'loop' position
        call random_number(rand1)
        call random_number(rand2)
        rand1=box_size*(rand1*2.-1.)/20.
        rand2=box_size*(rand2*2.-1.)/20.
        do i=old_pcount+1, pcount
          f(i)%x(1)=-box_size/2.
          f(i)%x(2)=radius*cos(pi*real(2*i-1)/inject_size)+rand1
          f(i)%x(3)=radius*sin(pi*real(2*i-1)/inject_size)+rand2
          if (i==old_pcount+1) then
            f(i)%behind=pcount ; f(i)%infront=i+1
          else if (i==pcount) then 
            f(i)%behind=i-1 ; f(i)%infront=old_pcount+1
          else
            f(i)%behind=i-1 ; f(i)%infront=i+1
          end if
          !zero the stored velocities
          f(i)%u1=0. ; f(i)%u2=0.
        end do
      case default
        call fatal_error('inject.mod','inject_type not set to useable value')    
    end select
  end subroutine
end module
