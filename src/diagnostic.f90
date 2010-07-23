module diagnostic
  !DIAGNOSTICS - MAINLY TO PRINT TO SCREEN
  use cdata
  use general
  contains
  !*************************************************
  subroutine velocity_info()
    implicit none
    real, allocatable :: uinfo(:,:)
    allocate(uinfo(pcount,2))
    uinfo(:,1)=sqrt(f(:)%u(1)**2+f(:)%u(2)**2+f(:)%u(3)**2)
    uinfo(:,2)=sqrt((f(:)%u1(1)-f(:)%u2(1))**2+&
                    (f(:)%u1(2)-f(:)%u2(2))**2+&
                    (f(:)%u1(3)-f(:)%u2(3))**2)
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
      !energy=energy+dot_product(f(i)%u,cross_product(f(i)%x,f(i)%ghosti-f(i)%x))
      energy=energy+dot_product(f(i)%u,cross_product(f(i)%x,sdot))*dist_gen(f(i)%x,f(i)%ghosti)
    end do
  end subroutine
  !*************************************************
  subroutine mean_curv()
    !caculate the mean curvature of the vortex system
    implicit none
    integer :: i
    kappa_bar=0.
    do i=1, pcount
      if (f(i)%infront==0) cycle !check for 'empty' particles
      kappa_bar=kappa_bar+curvature(i)
    end do
    !average this quantity
    kappa_bar=kappa_bar/count(mask=f(:)%infront>0)
  end subroutine
end module
