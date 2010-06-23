module periodic
  use cdata
  contains 
  !******************************************************************
  subroutine ghostp
    implicit none
    integer :: i
    do i=1, pcount
      if (f(i)%infront==0) cycle !empty particles
      call get_ghost_p(i,f(i)%ghosti, f(i)%ghostb)
    end do
  end subroutine
  !******************************************************************
  subroutine get_ghost_p(i,ginfront,gbehind)
    implicit none
    integer, intent(IN) :: i
    real :: ginfront(3), gbehind(3)
    ginfront(:)=f(f(i)%infront)%x(:)
    gbehind(:)=f(f(i)%behind)%x(:)
    !if periodic then must do more in here
    if (periodic_bc) then
      !we must ensure that ginfront/gbehind is not on the other side of the box
      !---------------------x------------------------------
      if ((f(i)%x(1)-ginfront(1))>(box_size/2.)) then
        ginfront(1)=ginfront(1)+box_size
      elseif ((f(i)%x(1)-ginfront(1))<(-box_size/2.)) then
        ginfront(1)=ginfront(1)-box_size
      end if
      if ((f(i)%x(1)-gbehind(1))>(box_size/2.)) then
        gbehind(1)=gbehind(1)+box_size
      elseif ((f(i)%x(1)-gbehind(1))<(-box_size/2.)) then
        gbehind(1)=gbehind(1)-box_size
      end if
      !---------------------y------------------------------
      if ((f(i)%x(2)-ginfront(2))>(box_size/2.)) then
        ginfront(2)=ginfront(2)+box_size
      elseif ((f(i)%x(2)-ginfront(2))<(-box_size/2.)) then
        ginfront(2)=ginfront(2)-box_size
      end if
      if ((f(i)%x(2)-gbehind(2))>(box_size/2.)) then
        gbehind(2)=gbehind(2)+box_size
      elseif ((f(i)%x(2)-gbehind(2))<(-box_size/2.)) then
        gbehind(2)=gbehind(2)-box_size
      end if
      !---------------------z------------------------------
      if ((f(i)%x(3)-ginfront(3))>(box_size/2.)) then
        ginfront(3)=ginfront(3)+box_size
      elseif ((f(i)%x(3)-ginfront(3))<(-box_size/2.)) then
        ginfront(3)=ginfront(3)-box_size
      end if
      if ((f(i)%x(3)-gbehind(3))>(box_size/2.)) then
        gbehind(3)=gbehind(3)+box_size
      elseif ((f(i)%x(3)-gbehind(3))<(-box_size/2.)) then
        gbehind(3)=gbehind(3)-box_size
      end if
    end if
  end subroutine
  !******************************************************************
  subroutine enforce_periodic()
    implicit none
    integer :: i
    do i=1, pcount
      if (f(i)%infront==0) cycle !empty particle
      !-------------x------------------     
      if (f(i)%x(1)>(box_size/2.)) then
        f(i)%x(1)=f(i)%x(1)-box_size
      else if (f(i)%x(1)<(-box_size/2.)) then
        f(i)%x(1)=f(i)%x(1)+box_size
      end if
      !-------------y------------------
      if (f(i)%x(2)>(box_size/2.)) then
        f(i)%x(2)=f(i)%x(2)-box_size
      else if (f(i)%x(2)<(-box_size/2.)) then
        f(i)%x(2)=f(i)%x(2)+box_size
      end if
      !-------------z------------------
      if (f(i)%x(3)>(box_size/2.)) then
        f(i)%x(3)=f(i)%x(3)-box_size
      else if (f(i)%x(3)<(-box_size/2.)) then
        f(i)%x(3)=f(i)%x(3)+box_size
      end if
    end do
  end subroutine
  
!this will contain all the periodic routines
!1. move particles from one side of a box to another
!2. a ghost particle routine
!   if not periodic then just copy all of the particles infront and behind
!   to a part of the qvort structure
!  type qvort !our main structure
!    real :: x(3) !position
!    real :: u(3), u1(3), u2(3) !stored velocities (adam bash)
!    real :: ghosti(3), ghostb(3) !ghost particles TO BE ADDED!!!!!
!    integer :: infront, behind !to form a line/loop
!  end type
!
!  if in periodic then have a routine which calculates the distance between particles infront 
!  and behind and if larger than the box size then recalculates the ghosti ghostb data
!  THINK ABOUT THIS ANDREW, PERHAPS THE OPTIMAL WAY IS TO DO THIS ONLY WHEN PARTICLES MOVE FROM 
!  ONE SIDE OF THE BOX TO THE OTHER
!  we then use these ghost particles in all our routines
end module
