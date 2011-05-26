module mirror
  !ALL THE ROUTINES REQUIRED TO IMPOSE CLOSED BOUNDARIES USING THE METHOD OF IMAGES
  !THIS IS DONE BY MIRRORING THE VORTICES AT EACH OF THE 6 EDGES OF THE BOX
  use cdata
  use general
  type(qvort), allocatable, private :: m(:,:) !the mirror array
  contains
  !**************************************************************************
  subroutine mirror_init
    !set up the mirror array
    implicit none
    integer :: i
    !note everything is in mirror_biot_savart is done using m(i)%infront
    !so copying beind to infront is OK, behind is never used
    !begin by allocating the mirror array
    allocate(m(6,pcount)) !6 faces of the cube
    !now we need to fill the array
    do i=1,6 !again 6 faces
      m(i,:)=f(:)
      if ((i==1).or.(i==2)) then
        if (i==1) then
          m(i,:)%x(1)=box_size-m(i,:)%x(1)
          m(i,:)%infront=m(i,:)%behind
          !reflect in positive x
        else
          m(i,:)%x(1)=-box_size-m(i,:)%x(1)
          m(i,:)%infront=m(i,:)%behind
          !reflect in negative x
        end if
      else if ((i==3).or.(i==4)) then
        if (i==3) then
          m(i,:)%x(2)=box_size-m(i,:)%x(2)
          m(i,:)%infront=m(i,:)%behind
          !reflect in positive y
        else
          m(i,:)%x(2)=-box_size-m(i,:)%x(2)
          m(i,:)%infront=m(i,:)%behind
          !reflect in negative y
        end if
      else !i=5 or 6
        if (i==5) then
          m(i,:)%x(3)=box_size-m(i,:)%x(3)
          m(i,:)%infront=m(i,:)%behind
          !reflect in positive z
        else
          m(i,:)%x(3)=-box_size-m(i,:)%x(3)
          m(i,:)%infront=m(i,:)%behind
          !reflect in negative z
        end if
      end if
    end do
  end subroutine
  !**************************************************************************
  subroutine mirror_close
    !deallocate the mirror array and print if asked
    implicit none
    integer, parameter :: mnum=6
    integer :: i
    character (len=40) :: print_file
    !if we want to print it goes in here
    if (mirror_print.and.(mod(itime, shots)==0)) then
      write(unit=print_file,fmt="(a,i4.4,a)")"./data/mirror",itime/shots,".dat"
      open(unit=98,file=print_file,status='replace',form='unformatted',access='stream')
        write(98) t
        write(98) pcount*mnum
        do i=1,mnum
          write(98) m(i,:)%x(1)
        end do
        do i=1,mnum
          write(98) m(i,:)%x(2)
        end do
        do i=1,mnum
          write(98) m(i,:)%x(3)
        end do
        do i=1,mnum
          write(98) m(i,:)%infront+int((i-1)*pcount)
        end do
      close(98)
    end if
    deallocate(m) !deallocate the array to avoid memory leak
  end subroutine
  !**************************************************************************
  subroutine biot_savart_mirror(i,u)
    !get the induced velocity from the method of images
    implicit none
    integer, intent(IN) :: i !the particle we want the velocity at
    real :: u(3) !the velocity at i (non-local mirror field)
    real :: u_mir(3) !helper array  
    real :: a_bs, b_bs, c_bs !helper variables
    integer :: j, k
    do k=6, 6 !all 6 faces
      do j=1, pcount
        !check that the particle is not empty
        if ((m(k,j)%infront==0)) cycle
        a_bs=dist_gen_sq(m(k,j)%x,f(i)%x) !distance squared between i and j on mirror
        b_bs=2.*dot_product((m(k,j)%x-f(i)%x),(m(k,m(k,j)%infront)%x-m(k,j)%x))
        c_bs=dist_gen_sq(m(k,m(k,j)%infront)%x,m(k,j)%x) !distance sqd between j, j+1
        !add non local contribution to velocity vector
        if (4*a_bs*c_bs-b_bs**2==0) cycle !avoid 1/0
        u_mir=cross_product((m(k,j)%x-f(i)%x),(m(k,m(k,j)%infront)%x-m(k,j)%x))
        u_mir=u_mir*quant_circ/((2*pi)*(4*a_bs*c_bs-b_bs**2))
        u_mir=u_mir*((2*c_bs+b_bs)/sqrt(a_bs+b_bs+c_bs)-(b_bs/sqrt(a_bs)))
        u=u+u_mir !add on the non-local contribution of j
      end do
    end do
  end subroutine
  !**************************************************************************
  subroutine mirror_flux_check(i,u)
    !ensure that the normal velocity at the boundaries is 0
    implicit none
    integer, intent(IN) :: i !the particle we want the velocity at
    real :: u(3) !the velocity at i
    real :: threshold !force zero flux at the boundaries
    threshold=box_size/2.-delta/2.
    if (f(i)%x(1)>threshold) u(1)=0.
    if (f(i)%x(1)<-threshold) u(1)=0.
    if (f(i)%x(2)>threshold) u(2)=0.
    if (f(i)%x(2)<-threshold) u(2)=0.
    if (f(i)%x(3)>threshold) u(3)=0.
    if (f(i)%x(3)<-threshold) u(3)=0.
  end subroutine
end module
