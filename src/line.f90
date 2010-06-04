module line
  use cdata
  use general
  use periodic
  contains
  !**************************************************************
  subroutine pinsert
    implicit none    
    type(qvort), allocatable, dimension(:) :: tmp
    real :: disti, curv, f_ddot(3)
    integer :: par_new
    integer :: old_pcount, i
    old_pcount=pcount
    do i=1, old_pcount
      if (i>size(f)) print*, 'I think there is a problem' ; exit
      if (f(i)%infront==0) cycle !empty particle
      !get the distance between the particle and the one infront
      disti=dist_gen(f(i)%x,f(i)%ghosti)
      if (disti>delta) then               
        print*, 'in here' ; stop
        !we need a new particle 
        !the first step is assess where to put the particle?
        !1. is there an empty slot in out array?
        if (minval(f(:)%infront)==0) then
          !there is an empty slot in the array
          !now find its location minloc is fortran intrinsic function
          par_new=minloc(f(:)%infront,1)
        else
          !2. we must resize the array - 
          !increment by 1 (speed) - use a dummy array tmp
          allocate(tmp(size(f)+1)) ; tmp(:)%infront=0 !0 the infront array
          !copy accross information
          tmp(1:size(f)) = f
          !now deallocate tmp and transfer particles back to f
          !move_alloc is new intrinsic function (fortran 2003)
          call move_alloc(from=tmp,to=f)
          !pcount+1 must be free - so put the particle there
          par_new=pcount+1 ; pcount=pcount+1 !increase the particle count
        end if
        !insert a new particle between i and infront using curvature
        !get second derivative at i
        call get_deriv_2(i,f_ddot) !general.mod
        curv=sqrt(dot_product(f_ddot,f_ddot)) !get the curvature
        f(par_new)%x=0.5*(f(i)%x+f(i)%ghosti)+&
                     (sqrt(curv**2-0.25*disti**2)-curv)*curv*f_ddot
        !0 the velocities
        f(par_new)%u=0. ; f(par_new)%u1=0. ; f(par_new)%u2=0.

        !set correct infront and behinds & ghostzones
        f(par_new)%behind=i ; f(par_new)%infront=f(i)%infront
        call get_ghost_p(par_new,f(par_new)%ghosti, f(par_new)%ghostb) !periodic.mod
        f(f(i)%infront)%behind=par_new
        call get_ghost_p(f(i)%infront,f(f(i)%infront)%ghosti, f(f(i)%infront)%ghostb) !periodic.mod         
        f(i)%infront=par_new
        call get_ghost_p(i,f(i)%ghosti, f(i)%ghostb) !periodic.mod         
      end if
    end do
  end subroutine
  !******************************************************************
  subroutine premove
    !maximum curvature cutoff based on an equilateral triangle with side
    !of length delta is \sqrt(3)/delta. This routine will remove particles
    !whose curvature is too high
    implicit none
    real :: distii
    integer :: infront, tinfront
    integer :: i    
    do i=1, pcount
      if (f(i)%infront==0) cycle !empty particle
      !get the distance between the particle and the one twice infront
      distii=distf(i,f(f(i)%infront)%infront)
      if (distii<delta) then
        infront=f(i)%infront ; tinfront=f(f(i)%infront)%infront
        !remove the particle at infront
        f(tinfront)%behind=i ; f(i)%infront=tinfront
        f(infront)%x(:)=0. ; f(infront)%infront=0
        f(infront)%behind=0
      end if
    end do
  end subroutine
end module
