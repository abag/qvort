!>all the routines used to keep the code periodic
module periodic
  use cdata
  use general
  contains 
  !******************************************************************
  !>dummy routine, calls get_ghost_p below
  subroutine ghostp
    implicit none
    integer :: i
    do i=1, pcount
      if (f(i)%infront==0) cycle !empty particles
      call get_ghost_p(i,f(i)%ghosti, f(i)%ghostb,f(i)%ghostii, f(i)%ghostbb)
    end do
  end subroutine
  !******************************************************************
  !>set the ghost particles, essentially these are the positions of the
  !!particles infront/behind and twice infront/behind
  !!if they are at the other side of the box 
  !!due to periodic b.c. this position must be adjusted
  subroutine get_ghost_p(i,ginfront,gbehind,giinfront,gbbehind)
    implicit none
    integer, intent(IN) :: i
    real :: ginfront(3), gbehind(3)
    real :: giinfront(3), gbbehind(3)
    ginfront(:)=f(f(i)%infront)%x(:)
    gbehind(:)=f(f(i)%behind)%x(:)
    giinfront(:)=f(f(f(i)%infront)%infront)%x(:)
    gbbehind(:)=f(f(f(i)%behind)%behind)%x(:)
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
      if ((f(i)%x(1)-giinfront(1))>(box_size/2.)) then
        giinfront(1)=giinfront(1)+box_size
      elseif ((f(i)%x(1)-giinfront(1))<(-box_size/2.)) then
        giinfront(1)=giinfront(1)-box_size
      end if
      if ((f(i)%x(1)-gbbehind(1))>(box_size/2.)) then
        gbbehind(1)=gbbehind(1)+box_size
      elseif ((f(i)%x(1)-gbbehind(1))<(-box_size/2.)) then
        gbbehind(1)=gbbehind(1)-box_size
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
      if ((f(i)%x(2)-giinfront(2))>(box_size/2.)) then
        giinfront(2)=giinfront(2)+box_size
      elseif ((f(i)%x(2)-giinfront(2))<(-box_size/2.)) then
        giinfront(2)=giinfront(2)-box_size
      end if
      if ((f(i)%x(2)-gbbehind(2))>(box_size/2.)) then
        gbbehind(2)=gbbehind(2)+box_size
      elseif ((f(i)%x(2)-gbbehind(2))<(-box_size/2.)) then
        gbbehind(2)=gbbehind(2)-box_size
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
      if ((f(i)%x(3)-giinfront(3))>(box_size/2.)) then
        giinfront(3)=giinfront(3)+box_size
      elseif ((f(i)%x(3)-giinfront(3))<(-box_size/2.)) then
        giinfront(3)=giinfront(3)-box_size
      end if
      if ((f(i)%x(3)-gbbehind(3))>(box_size/2.)) then
        gbbehind(3)=gbbehind(3)+box_size
      elseif ((f(i)%x(3)-gbbehind(3))<(-box_size/2.)) then
        gbbehind(3)=gbbehind(3)-box_size
      end if
    else if (periodic_bc_notx) then
      !this could be neater!!!
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
      if ((f(i)%x(2)-giinfront(2))>(box_size/2.)) then
        giinfront(2)=giinfront(2)+box_size
      elseif ((f(i)%x(2)-giinfront(2))<(-box_size/2.)) then
        giinfront(2)=giinfront(2)-box_size
      end if
      if ((f(i)%x(2)-gbbehind(2))>(box_size/2.)) then
        gbbehind(2)=gbbehind(2)+box_size
      elseif ((f(i)%x(2)-gbbehind(2))<(-box_size/2.)) then
        gbbehind(2)=gbbehind(2)-box_size
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
      if ((f(i)%x(3)-giinfront(3))>(box_size/2.)) then
        giinfront(3)=giinfront(3)+box_size
      elseif ((f(i)%x(3)-giinfront(3))<(-box_size/2.)) then
        giinfront(3)=giinfront(3)-box_size
      end if
      if ((f(i)%x(3)-gbbehind(3))>(box_size/2.)) then
        gbbehind(3)=gbbehind(3)+box_size
      elseif ((f(i)%x(3)-gbbehind(3))<(-box_size/2.)) then
        gbbehind(3)=gbbehind(3)-box_size
      end if
    end if
    if (mirror_bc) then
      if (f(i)%pinnedi) then
        if (f(i)%x(1)>(box_size/2.-delta))  ginfront(1)=box_size-ginfront(1)
        if (f(i)%x(1)<-(box_size/2.-delta)) ginfront(1)=-box_size-ginfront(1)
        if (f(i)%x(2)>(box_size/2.-delta))  ginfront(2)=box_size-ginfront(2)
        if (f(i)%x(2)<-(box_size/2.-delta)) ginfront(2)=-box_size-ginfront(2)
        if (f(i)%x(3)>(box_size/2.-delta))  ginfront(3)=box_size-ginfront(3)
        if (f(i)%x(3)<-(box_size/2.-delta)) ginfront(3)=-box_size-ginfront(3)
      end if
      if (f(i)%pinnedb) then
        if (f(i)%x(1)>(box_size/2.-delta))  gbehind(1)=box_size-gbehind(1)
        if (f(i)%x(1)<-(box_size/2.-delta)) gbehind(1)=-box_size-gbehind(1)
        if (f(i)%x(2)>(box_size/2.-delta))  gbehind(2)=box_size-gbehind(2)
        if (f(i)%x(2)<-(box_size/2.-delta)) gbehind(2)=-box_size-gbehind(2)
        if (f(i)%x(3)>(box_size/2.-delta))  gbehind(3)=box_size-gbehind(3)
        if (f(i)%x(3)<-(box_size/2.-delta)) gbehind(3)=-box_size-gbehind(3)
      end if
    end if
  end subroutine
  !******************************************************************
  !>if a point/particle leaves one side of the box, 
  !>reinsert it on the opposite side
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
  !******************************************************************
  !>if a point/particle leaves one side of the box, 
  !>reinsert it on the opposite side - do not do the x direction
  !> in the x direction we simple remove the full loop
  subroutine enforce_periodic_yz()
    implicit none
    integer :: i
    do i=1, pcount
      if (f(i)%infront==0) cycle !empty particle
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
      !-------------x------------------     
      if (f(i)%x(1)>(box_size/2.)) then
        call boundary_loop_remove(i)
      else if (f(i)%x(1)<(-box_size/2.)) then
        call boundary_loop_remove(i)
      end if
    end do
  end subroutine
  !******************************************************************
  !> remove loops which hit the boundaries
  subroutine enforce_open_removal()
    implicit none
    integer :: i
    do i=1, pcount
      if (f(i)%infront==0) cycle !empty particle
      !-------------x------------------     
      if (f(i)%x(1)>(box_size/2.)) then
        call boundary_loop_remove(i)
      else if (f(i)%x(1)<(-box_size/2.)) then
        call boundary_loop_remove(i)
      end if
      !-------------y------------------
      if (f(i)%x(2)>(box_size/(2.*xdim_scaling_factor))) then
        call boundary_loop_remove(i)
      else if (f(i)%x(2)<(-box_size/(2.*xdim_scaling_factor))) then
        call boundary_loop_remove(i)
      end if
      !-------------z------------------
      if (f(i)%x(3)>(box_size/(2.*xdim_scaling_factor))) then
        call boundary_loop_remove(i)
      else if (f(i)%x(3)<(-box_size/(2.*xdim_scaling_factor))) then
        call boundary_loop_remove(i)
      end if
    end do
  end subroutine
  !**************************************************
  !>remove loops that have left the box in \pm x directions
  !>this is very similar to loop_killer in line.f90
  !!a better way would be have a separate loop count
  !!and loop removal routine in general.mod
  subroutine boundary_loop_remove(particle)
    implicit none
    integer :: particle, next
    integer :: store_next
    integer :: i, counter
    next=particle 
    do i=1, pcount
      store_next=f(next)%infront
      call clear_particle(next) !general.mod
      next=store_next
      if (next==particle) then
        exit  
      end if
    end do
  end subroutine
end module
