module diagnostic
  use cdata
  use general
  contains
  !*************************************************
  subroutine velocity_info()
    implicit none
    real, allocatable :: uinfo(:,:)
    integer :: i
    allocate(uinfo(pcount,2))
    uinfo(:,1)=sqrt(f(:)%u(1)**2+f(:)%u(1)**2+f(:)%u(1)**2)
    uinfo(:,2)=sqrt((f(:)%u1(1)-f(i)%u2(:))**2+&
                    (f(:)%u1(2)-f(i)%u2(:))**2+&
                    (f(:)%u1(3)-f(i)%u2(:))**2)
    maxu=maxval(uinfo(:,1)) ; maxdu=maxval(uinfo(:,2))
    deallocate(uinfo)
  end subroutine
  !*************************************************
  subroutine energy_info()
    implicit none
    real :: sdot(3)
    integer :: i
    energy=0.
    do i=1, pcount
      if (f(i)%infront==0) cycle !check for 'empty' particles
      call get_deriv_1(i,sdot)
      energy=energy+distf(i,f(i)%infront)*dot_product(f(i)%u,cross_product((f(i)%x-dt*f(i)%u),sdot))
    end do
  end subroutine
end module
