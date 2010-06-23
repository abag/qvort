module periodic
  !ALL THE ROUTINES NEEDED TO KEEP THE CODE PERIODIC
  use cdata
  contains 
  !******************************************************************
  subroutine ghostp
    !dummy routine, calls get_ghost_p below
    implicit none
    integer :: i
    do i=1, pcount
      if (f(i)%infront==0) cycle !empty particles
      call get_ghost_p(i,f(i)%ghosti, f(i)%ghostb)
    end do
  end subroutine
  !******************************************************************
  subroutine get_ghost_p(i,ginfront,gbehind)
    !set the ghost particles, essentially these are the positions of the
    !particles infront/behind, if they are at the other side of the box 
    !due to periodic b.c. this position must be adjusted
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
    !if a particle leaves one side of the box, reinsert it on the opposite side
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
end module
