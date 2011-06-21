!> contains a subroutine to calculate kernel density estimation
!> an improvement on a histogram. Also contains routines to perform
!> a quick sorting of the data which is required by the algorithm
module kernel_density
  contains
  !> subroutine for kernel density estimation with iterative plug-in bandwidth selection,
  !> only uses the gauss kernel and estimates only the density itself and not its derivatives.
  !> @param x(n)   input|dp - sorted data  
  !> @param n  input|integer - length of x
  !> @param z(m)  input|dp - output grid (sorted)
  !> @param m  input|integer - length of z
  !> @param f(m)  output|dp - estimated density
  !> @param h  output|dp - estimated iterative plugin bandwidth
  subroutine get_kernel_density(x,n,z,m,f,h)
    implicit none
    integer, parameter  :: dp = selected_real_kind(12, 60)
    
    real (dp), intent(in)   :: x(:)
    integer, intent(in)     :: n
    real (dp), intent(in)   :: z(:)
    integer, intent(in)     :: m
    real (dp), intent(out)  :: f(:)
    real (dp), intent(out)  :: h
    
    ! local variables
    integer    :: i, it, iter, j, jbegin, jend
    real (dp)  :: a, co1, const, d2, d3, e2, e3, h2, h3, pi, r2, rhat2, rhat3, &
                  rt2pi, s, s2, s3, t, xiqr, xn
    
    xn = n
    xiqr = x(nint(.75*xn)) - x(nint(.25*xn)+1)
    pi = 3.14159265358979324_dp
    rt2pi = sqrt(2._dp*pi)
    iter = 5
    !-------  estimate inflation constant c
    h2 = (.920*xiqr) / (xn**(1._dp/7._dp))
    h3 = (.912*xiqr) / (xn**(1._dp/9._dp))
    s2 = 0._dp
    s3 = 0._dp
    loop20:  do  i = 1, n - 1
      do  j = i + 1, n
        d2 = ((x(i)-x(j))/h2) ** 2
        d3 = ((x(i)-x(j))/h3) ** 2
        if (d2 > 50 .and. d3 > 60) cycle loop20
        e2 = exp(-0.5_dp*d2)
        e3 = exp(-0.5_dp*d3)
        s2 = s2 + (d2**2 - 6*d2 + 3._dp) * e2
        s3 = s3 + (d3**3 - 15*d3**2 + 45*d3 - 15) * e3
      end do
    end do loop20
    rhat2 = (2.0_dp*s2) / ((xn**2)*(h2**5)*rt2pi)
    rhat2 = rhat2 + 3._dp / (rt2pi*xn*(h2**5))
    rhat3 = (-2.0_dp*s3) / ((xn**2)*(h3**7)*rt2pi)
    rhat3 = rhat3 + 15.0_dp / (rt2pi*xn*(h3**7))
    co1 = 1.357_dp * (rhat2/rhat3) ** (1.0_dp/7.0_dp)
    !-
    !-------  compute constant of asymptotic formula
    !-
    const = 1._dp / (2._dp*sqrt(pi))
    a = 1.132795764 / rhat3 ** (1._dp/7._dp) * xn ** (-1.0_dp/2.0_dp)
    !-
    !------  loop over iterations
    !-
    do  it = 1, iter
    !-
    !-------  estimate functional
      
      s = 0._dp
      loop40:  do  i = 1, n - 1
        do  j = i + 1, n
          d2 = ((x(i) - x(j))/a) ** 2
          if (d2 > 50) cycle loop40
          e2 = exp(-0.5_dp*d2)
          s = s + (d2**2 - 6._dp*d2 + 3._dp) * e2
        end do
      end do loop40
      r2 = (2.0_dp*s) / ((xn**2)*(a**5)*rt2pi)
      r2 = r2 + 3._dp / (rt2pi*xn*(a**5))
    !-
    !-------  estimate bandwidth by asymptotic formula
    !-
      h = (const/(r2*xn)) ** (0.2_dp)
      a = co1 * h ** (5.0_dp/7.0_dp)
    end do
    !-
    !------- estimate density with plugin bandwidth
    !-
    jbegin = 1
    jend = 1
    do  i = 1, m
      s = 0._dp
      do  j = jbegin, jend
        t = (z(i) - x(j)) / h
        if (t > 5.0 .and. jbegin < n) then
          jbegin = jbegin + 1
          cycle
        end if
        s = s + exp(-t*t/2.)
      end do
      do  jend = j, n
        t = (z(i) - x(jend)) / h
        if (t < -5.0_dp) exit
        s = s + exp(-t*t/2.)
      end do
      f(i) = s / (xn*h*rt2pi)
      jend = jend - 1
    !-
    end do  
    return
  end subroutine
  !********************************************************** 
  recursive subroutine qsort(A)
    real, intent(in out), dimension(:) :: A
    integer :: iq
    if(size(A) > 1) then
       call partition(A, iq)
       call qsort(A(:iq-1))
       call qsort(A(iq:))
    endif
  end subroutine
  !**********************************************************
  subroutine partition(A, marker)
    real, intent(in out), dimension(:) :: A
    integer, intent(out) :: marker
    integer :: i, j
    real :: temp
    real :: x      ! pivot point
    x = A(1)
    i= 0
    j= size(A) + 1
    do
       j = j-1
       do
          if (A(j) <= x) exit
          j = j-1
       end do
       i = i+1
       do
          if (A(i) >= x) exit
          i = i+1
       end do
       if (i < j) then
          ! exchange A(i) and A(j)
          temp = A(i)
          A(i) = A(j)
          A(j) = temp
       elseif (i == j) then
          marker = i+1
          return
       else
          marker = i
          return
       endif
    end do
  end subroutine
end module
