module general
  !general routines e.g. derivatives, vector products, distance evaluation
  use cdata
  contains
  subroutine get_deriv_1(i,s_dot)
    implicit none
    integer, intent(IN) :: i
    real :: s_dot(3), disti, distb
    disti=dist_gen(f(i)%x,f(i)%ghosti) ; distb=dist_gen(f(i)%x,f(i)%ghostb)
    s_dot(1:3)=distb*f(i)%ghosti(1:3)+ &
               (disti-distb)*f(i)%x(1:3)- &
               disti*f(i)%ghostb(1:3)
    s_dot=s_dot/(2.*distb*disti)
  end subroutine
  !*********************************************************************
  subroutine get_deriv_2(i,s_ddot)
    implicit none
    integer, intent(IN) :: i
    real :: s_ddot(3), disti, distb
    disti=dist_gen(f(i)%x,f(i)%ghosti) ; distb=dist_gen(f(i)%x,f(i)%ghostb)
    s_ddot(1:3)=2.*(f(i)%ghosti(1:3)/(disti*(disti+distb))- &
                   f(i)%x(1:3)/(disti*distb)+ &
                   f(i)%ghostb(1:3)/(distb*(disti+distb)))
  end subroutine
  !*********************************************************************
  subroutine clear_particle(i)
    !empty all a particles information
    implicit none
    integer, intent(IN) :: i
    f(i)%x=0.
    f(i)%u=0. ; f(i)%u1=0. ; f(i)%u2=0.
    f(i)%ghosti=0. ; f(i)%ghostb=0.
    f(i)%infront=0 ; f(i)%behind=0
    f(i)%closest=0 ; f(i)%closestd=0.
  end subroutine
  !*********************************************************************
  real function distf(i,j)
    !calculate the distance between particles in the f vector
    use Cdata
    implicit none
    integer, intent(IN) :: i, j
    distf=sqrt((f(i)%x(1)-f(j)%x(1))**2+&
               (f(i)%x(2)-f(j)%x(2))**2+&
               (f(i)%x(3)-f(j)%x(3))**2)
  end function
  !*********************************************************************
  real function distfsq(i,j)
    !calculate the squared distance between particles in the f vector
    use Cdata
    implicit none
    integer, intent(IN) :: i, j
    distfsq=(f(i)%x(1)-f(j)%x(1))**2+&
            (f(i)%x(2)-f(j)%x(2))**2+&
            (f(i)%x(3)-f(j)%x(3))**2
  end function
  !*********************************************************************
  real function dist_gen(a,b)
    !calculate the distance between points a(1:3), b(1:3)
    use Cdata
    implicit none
    real, dimension(3), intent(IN) :: a, b
    dist_gen=sqrt((a(1)-b(1))**2+&
                  (a(2)-b(2))**2+&
                  (a(3)-b(3))**2)
  end function
  !*********************************************************************
  real function dist_gen_sq(a,b)
    !calculate the squared distance between points a(1:3), b(1:3)
    use Cdata
    implicit none
    real, dimension(3), intent(IN) :: a, b
    dist_gen_sq=(a(1)-b(1))**2+&
                (a(2)-b(2))**2+&
                (a(3)-b(3))**2
  end function
  !*********************************************************************
  real function tangentf(i)
    !calculate the tangent vector at point i
    use Cdata
    implicit none
    real, dimension(3) :: tangentf
    integer, intent(IN) :: i
    real :: dist
    !get distance between particle and ghost particle (infront)
    dist=dist_gen(f(i)%x,f(i)%ghosti)
    !now determine the vector (low order finite diff.)
    tangentf(:)=(f(i)%ghosti(:)-f(i)%x(:))/dist
  end function
  !*********************************************************************
  real function norm_tanf(i)
    !calculate the normalised tangent vector at point i
    use Cdata
    implicit none
    real, dimension(3) :: norm_tanf
    integer, intent(IN) :: i
    real :: length !length of vector
    !now determine the vector (low order finite diff.)
    norm_tanf(:)=(f(i)%ghosti(:)-f(i)%x(:))
    !calculate length of vector
    length=sqrt(norm_tanf(1)**2+norm_tanf(2)**2+norm_tanf(3)**2)
    norm_tanf(:)=norm_tanf(:)/length !normalise
  end function
  !*********************************************************************
  real function cross_product(a,b)
    !calculate the cross product of a and b
    implicit none
    real, dimension(3) :: cross_product
    real, dimension(3), intent(IN) :: a, b
    cross_product(1)=a(2)*b(3)-a(3)*b(2)
    cross_product(2)=a(3)*b(1)-a(1)*b(3)
    cross_product(3)=a(1)*b(2)-a(2)*b(1)
  end function
  !*********************************************************************
  subroutine zero_finder(location)
    !only used for code testing
    implicit none
    integer :: i, zcount=0
    character(len=*) :: location
    do i=1, pcount
      if (f(i)%infront==0) zcount=zcount+1
    end do
    write(*,*) 'zero finder called at', trim(location)
    write(*,*) 'zero count= ', zcount
    write(*,*) 'ending run...' ; stop 
  end subroutine
end module
