!>general routines e.g. derivatives, vector products, distance evaluation
!>this is used by all other modules (excluding cdata)
module general
  use cdata
  contains
  !> get the first derivative at position i using a second order adaptive mesh
  !>finite difference scheme:
  !!\f[
  !!\frac{d \mathbf{s}_i}{d \xi}=\frac{\ell_{i-1}\mathbf{s}_{i+1}+(\ell_{i+1}-\ell_{i-1})\mathbf{s}_i+\ell_{i+1}\mathbf{s}_{i-1}}
  !!  {2\ell_{i+1}\ell_{i-1}}+{\cal O}(\ell^2)
  !! \f]
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
  !>get the second derivative at position i using a second order adaptive 
  !>mesh finite difference scheme:
  !!\f[
  !!\frac{d^2 \mathbf{s}_i}{d \xi^2}=\frac{2\mathbf{s}_{i+1}}{\ell_{i+1}(\ell_{i+1}+\ell_{i-  1})}-
  !!\frac{2\mathbf{s}_i}{\ell_{i+1}\ell_{i-1}}+\frac{2\mathbf{s}_{i-1}}{\ell_{i-1}(\ell_{i+  1}+\ell_{i-1})}+{\cal O}(\ell^2)
  !! \f]
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
  !>empty all a points information - used if a point is removed
  !>(e.g. contraction of filament)
  subroutine clear_particle(i)
    implicit none
    integer, intent(IN) :: i
    f(i)%x=0.
    f(i)%u=0. ; f(i)%u1=0. ; f(i)%u2=0.
    f(i)%ghosti=0. ; f(i)%ghostb=0.
    f(i)%infront=0 ; f(i)%behind=0
    f(i)%closest=0 ; f(i)%closestd=0.
  end subroutine
  !*********************************************************************
  !>calculate the distance between points in the f vector
  !!The distance between \f$(x_1,y_1,z_1)\f$ and \f$(x_2,y_2,z_2)\f$ is 
  !!\f$\sqrt{(x_2-x_1)^2+(y_2-y_1)^2+(z_2-z_1)^2}\f$.
  real function distf(i,j)
    use Cdata
    implicit none
    integer, intent(IN) :: i, j
    distf=sqrt((f(i)%x(1)-f(j)%x(1))**2+&
               (f(i)%x(2)-f(j)%x(2))**2+&
               (f(i)%x(3)-f(j)%x(3))**2)
  end function
  !*********************************************************************
  !>calculate the squared distance between particles in the f vector
  !!\f${(x_2-x_1)^2+(y_2-y_1)^2+(z_2-z_1)^2}\f$.
  real function distfsq(i,j)
    use Cdata
    implicit none
    integer, intent(IN) :: i, j
    distfsq=(f(i)%x(1)-f(j)%x(1))**2+&
            (f(i)%x(2)-f(j)%x(2))**2+&
            (f(i)%x(3)-f(j)%x(3))**2
  end function
  !*********************************************************************
  !>calculate the curvature at the particle i: \f$|\mathbf{s}''|\f$
  real function curvature(i)
    use Cdata
    implicit none
    integer, intent(IN) :: i
    real :: fddot(3)
    call get_deriv_2(i,fddot)
    curvature=sqrt(dot_product(fddot,fddot))
  end function
  !*********************************************************************
  !>calculate the distance between points \f$\mathbf{a}\f$ and \f$\mathbf{b}\f$
  !>again just pythagoras
  real function dist_gen(a,b)
    use Cdata
    implicit none
    real, dimension(3), intent(IN) :: a, b
    dist_gen=sqrt((a(1)-b(1))**2+&
                  (a(2)-b(2))**2+&
                  (a(3)-b(3))**2)
  end function
  !*********************************************************************
  !>calculate the squared distance between points \f$\mathbf{a}\f$ and \f$\mathbf{b}\f$
  real function dist_gen_sq(a,b)
    use Cdata
    implicit none
    real, dimension(3), intent(IN) :: a, b
    dist_gen_sq=(a(1)-b(1))**2+&
                (a(2)-b(2))**2+&
                (a(3)-b(3))**2
  end function
  !*********************************************************************
  !>calculate the angle between vectors \f$\mathbf{a}\f$ and  \f$\mathbf{b}\f$
  !!\f[ \theta=\cos^{-1} \frac{ \mathbf{a} \cdot \mathbf{b}}{ab} \f]
  real function vector_angle(a,b)
    use Cdata
    implicit none
    real, dimension(3), intent(IN) :: a, b
    vector_angle=acos(dot_product(a,b)/ &
                 (sqrt(dot_product(a,a))*sqrt(dot_product(b,b))))
  end function
  !*********************************************************************
  !>low order finite difference to calculate the tangent vector at point i
  !!\f[
  !!\frac{d \mathbf{s}_i}{d \xi}=\frac{\mathbf{s}_{i+1}-\mathbf{s}_{i}}{\ell_i}
  !! \f]
  real function tangentf(i)
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
  !>calculate the normalised tangent vector at point i low order fininte diff.
  real function norm_tanf(i)    
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
  !>calculate the normalised normal vector at point i:
  !>\f$\hat{\mathbf{s}''}\f$
  real function normalf(i)
    use Cdata
    implicit none
    real, dimension(3) :: normalf
    integer, intent(IN) :: i
    real :: length !length of vector
    !get the second derivative
    call get_deriv_2(i,normalf)
    !calculate length of vector
    length=sqrt(normalf(1)**2+normalf(2)**2+normalf(3)**2)
    normalf(:)=normalf(:)/length !normalise
  end function
  !*********************************************************************
  !>calculate the normalised binormal vector at point i:
  !>\f$\hat{\mathbf{s}'} \times \hat{\mathbf{s}''}\f$
  real function binormalf(i)
    use Cdata
    implicit none
    real, dimension(3) :: s_dot, s_ddot, binormalf
    integer, intent(IN) :: i
    real :: length !length of vector
    !get the first/second derivative
    call get_deriv_1(i,s_dot)
    call get_deriv_2(i,s_ddot)
    binormalf=cross_product(s_dot,s_ddot)
    !calculate length of vector
    length=sqrt(binormalf(1)**2+binormalf(2)**2+binormalf(3)**2)
    binormalf(:)=binormalf(:)/length !normalise
  end function
  !*********************************************************************
  !>calculate the cross product of \f$\mathbf{a}\f$ and \f$\mathbf{b}\f$
  real function cross_product(a,b)
    implicit none
    real, dimension(3) :: cross_product
    real, dimension(3), intent(IN) :: a, b
    cross_product(1)=a(2)*b(3)-a(3)*b(2)
    cross_product(2)=a(3)*b(1)-a(1)*b(3)
    cross_product(3)=a(1)*b(2)-a(2)*b(1)
  end function
  !*********************************************************************
  !> a routine purely used for code testing finds empty particles 
  subroutine zero_finder(location)
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
