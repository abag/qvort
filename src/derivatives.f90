!>routines to calculate spatial derivatives
module derivatives
  use cdata
  contains
  !> get the first derivative at position i using a second/fourth order adaptive mesh
  !>finite difference scheme:
  !!\f[
  !!\frac{d \mathbf{s}_i}{d \xi}=\frac{\ell_{i-1}\mathbf{s}_{i+1}+(\ell_{i+1}-\ell_{i-1})\mathbf{s}_i+\ell_{i+1}\mathbf{s}_{i-1}}
  !!  {2\ell_{i+1}\ell_{i-1}}+{\cal O}(\ell^2)
  !! \f]
  subroutine get_deriv_1(i,s_dot)
    implicit none
    integer, intent(IN) :: i
    real :: s_dot(3), disti, distb
    real :: disti_ii, distb_bb
    real :: A, B, C, D ,E !helper variables
    select case(deriv_order)
      case('second')
      disti=dist_gen(f(i)%x,f(i)%ghosti) ; distb=dist_gen(f(i)%x,f(i)%ghostb)
      s_dot(1:3)=distb*f(i)%ghosti(1:3)+ &
                (disti-distb)*f(i)%x(1:3)- &
                 disti*f(i)%ghostb(1:3)
      s_dot=s_dot/(2.*distb*disti)
      case('fourth')
      disti=dist_gen(f(i)%x,f(i)%ghosti) ; distb=dist_gen(f(i)%x,f(i)%ghostb)
      disti_ii=dist_gen(f(i)%ghostii,f(i)%ghosti) 
      distb_bb=dist_gen(f(i)%ghostbb,f(i)%ghostb)
      A=distb*(disti**2)+distb*disti*disti_ii
      A=A/(distb_bb*(distb_bb+distb)*(distb_bb+distb+disti)*(distb_bb+distb+disti+disti_ii))
      B=-distb_bb*(disti**2)-distb*(disti**2)-distb_bb*disti*disti_ii-distb*disti*disti_ii
      B=B/(distb_bb*distb*(distb+disti)*(distb+disti+disti_ii))
      D=distb_bb*distb*disti+(distb**2)*disti+distb_bb*distb*disti_ii+(distb**2)*disti_ii
      D=D/(disti*disti_ii*(distb+disti)*(distb_bb+distb+disti))
      E=-disti*(distb**2)-distb_bb*distb*disti
      E=E/(disti_ii*(disti+disti_ii)*(distb+disti+disti_ii)*(distb_bb+distb+disti+disti_ii))
      C=-(A+B+D+E)
      s_dot(1:3)=A*f(i)%ghostbb+B*f(i)%ghostb+C*f(i)%x+D*f(i)%ghosti+E*f(i)%ghostii
    end select
  end subroutine
  !*********************************************************************
  !>get the second derivative at position i using a second/fourth order adaptive 
  !>mesh finite difference scheme:
  !!\f[
  !!\frac{d^2 \mathbf{s}_i}{d \xi^2}=\frac{2\mathbf{s}_{i+1}}{\ell_{i+1}(\ell_{i+1}+\ell_{i-  1})}-
  !!\frac{2\mathbf{s}_i}{\ell_{i+1}\ell_{i-1}}+\frac{2\mathbf{s}_{i-1}}{\ell_{i-1}(\ell_{i+  1}+\ell_{i-1})}+{\cal O}(\ell^2)
  !! \f]
  subroutine get_deriv_2(i,s_ddot)
    implicit none
    integer, intent(IN) :: i
    real :: s_ddot(3), disti, distb
    real :: disti_ii, distb_bb
    real :: A, B, C, D ,E !helper variables
    select case(deriv_order)
      case('second')
      disti=dist_gen(f(i)%x,f(i)%ghosti) ; distb=dist_gen(f(i)%x,f(i)%ghostb)
      s_ddot(1:3)=2.*(f(i)%ghosti(1:3)/(disti*(disti+distb))- &
                      f(i)%x(1:3)/(disti*distb)+ &
                      f(i)%ghostb(1:3)/(distb*(disti+distb)))
      case('fourth')
      disti=dist_gen(f(i)%x,f(i)%ghosti) ; distb=dist_gen(f(i)%x,f(i)%ghostb)
      disti_ii=dist_gen(f(i)%ghostii,f(i)%ghosti) 
      distb_bb=dist_gen(f(i)%ghostbb,f(i)%ghostb)
      A=2.*(-2.*distb*disti+disti**2-distb*disti_ii+disti*disti_ii)
      A=A/(distb_bb*(distb_bb+distb)*(distb_bb+distb+disti)*(distb_bb+distb+disti+disti_ii))
      B=2.*(2.*distb_bb*disti+2.*distb*disti-disti**2+distb_bb*disti_ii+distb*disti_ii-disti*disti_ii)
      B=B/(distb_bb*distb*(distb+disti)*(distb+disti+disti_ii))
      D=2.*(-distb_bb*distb-distb**2+distb_bb*disti+2.*distb*disti+distb_bb*disti_ii+2.*distb*disti_ii)
      D=D/(disti*disti_ii*(distb+disti)*(distb_bb+distb+disti))
      E=2.*(distb_bb*distb+distb**2-distb_bb*disti-2.*distb*disti)
      E=E/(disti_ii*(disti+disti_ii)*(distb+disti+disti_ii)*(distb_bb+distb+disti+disti_ii))
      C=-(A+B+D+E)
      !print*, A, B, C, D, E
      s_ddot(1:3)=A*f(i)%ghostbb+B*f(i)%ghostb+C*f(i)%x+D*f(i)%ghosti+E*f(i)%ghostii    
    end select
  end subroutine
  !*********************************************************************
  !>calculate the distance between points \f$\mathbf{a}\f$ and \f$\mathbf{b}\f$
  !>just pythagoras
  real function dist_gen(a,b)
    use Cdata
    implicit none
    real, dimension(3), intent(IN) :: a, b
    dist_gen=sqrt((a(1)-b(1))**2+&
                  (a(2)-b(2))**2+&
                  (a(3)-b(3))**2)
  end function
end module

