!>statistical routines for generating random numbers from
!!a variety of statistical distributions, as well as evaluating
!!pdf's
module statistics
  use cdata
  contains
  !*********************************************************************
  !**************ROUTINES TO GENERATE RANDOM VARIABLES******************
  !*********************************************************************
  !>used to draw random variables from a \f$N(0,1)\f$ distribution, 
  !>uses box-muller algorithm
  real function rnorm(mu,sigma2)
    implicit none
    real, intent(IN) :: mu, sigma2
    real :: u1, u2
    call random_number(u1) ;  call random_number(u2)
    rnorm=mu+sqrt(sigma2)*sqrt(-2.*log(u1))*cos(2.*pi*u2)   
  end function
  !*********************************************************************
  !>used to draw random variables from laplace distribution
  !!http://en.wikipedia.org/wiki/Laplace_distributio
  real function rlaplace(mu,b) !NOT TESTED
    implicit none
    real, intent(IN) :: mu, b
    real :: u
    call random_number(u)
    u=u-0.5
    rlaplace=mu-b*sign(1.,u)*log(1.-2.*abs(u))
  end function
  !*********************************************************************
  !>draw random variables from a \f$U(\alpha,\beta)\f$ distribution
  real function runif(alpha,beta)
    implicit none
    real, intent(IN) :: alpha, beta
    real :: u
    call random_number(u)
    runif=alpha+u*(beta-alpha)
  end function
  !*********************************************************************
  !>draw random variables from \f$\Gamma(a,b)\f$ distribution 
  real function rgamma(a,b)
    implicit none
    real, intent(IN) :: a, b
    real :: uvect(floor(a))
    real :: delta, nu
    real :: eta, xi
    real :: u1, u2
    real :: var
    !set up parameters
    delta=a-floor(a)
    nu=exp(1.)/(exp(1.)+delta)
    rgamma=0.
    if (delta>0.) then
      !now need a loop
      do 
        call random_number(u1) ; call random_number(u2)
        if (u1<=nu) then
          xi=(u1/nu)**(1./delta)
          eta=u2*(xi**(delta-1.))
        else
         xi=1.-log((u1-nu)/(1.-nu))
         eta=u2*exp(-xi)
        end if
        if (eta<(xi**(delta-1.))*exp(-xi)) then
          rgamma=xi
          exit
        end if
      end do
    end if
    call random_number(uvect) 
    rgamma=b*(rgamma-sum(log(uvect)))
  end function
  !*********************************************************************
  !>used to draw random variables from \f$\exp(\lambda)\f$ distribution
  real function rexp(lambda)
    implicit none
    real, intent(IN) :: lambda
    real :: u
    call random_number(u)
    rexp=-log(u)/lambda
  end function
  !*********************************************************************
  !*************ROUTINES TO EVALUATE PROB. DENSITY FUNCS.***************
  !*********************************************************************
  !>PDF for normal - \f[ f(x) = \frac{1}{\sqrt{2\pi\sigma^2}} e^{-\frac{(x-\mu)^2}{2\sigma^2}}\f] 
  real function dnorm(x,mu,sigma2)
    implicit none
    real, intent(IN) :: x, mu, sigma2
    dnorm=(1./(2.*pi*(sigma2)))*exp(-(x-mu)**2/(2.*(sigma2)))
  end function
  !*********************************************************************
  !>PDF for gamma - \f[ f(x;a,b) = x^{a-1} \frac{e^{-x/b}}{b^a \, \Gamma(a)}\f] 
  real function dgamma(x,a,b)
      implicit none
      real, intent(IN) :: x,a, b
      if (x<0) then
        dgamma=0.
      else
        dgamma=(x**(a-1))*exp(-x/b)/((b**a)*gamma(a))
      end if
  end function
  !**********************************************************************
  !>PDF for cauchy - \f[ f(x; x_0,\gamma) = \frac{1}{\pi \gamma \left[1 + \left(\frac{x - x_0}{\gamma}\right)^2 \right]} \f] 
  real function dcauchy(x,x0,gamm)
    !PDF for cauchy 
    implicit none
    real, intent(IN) :: x, x0, gamm
    dcauchy=(1./pi)*(gamm/((x-x0)**2+gamm**2))
  end function
  !**********************************************************************
  !>PDF for exponential distribution -
  !!\f[ f(x;\lambda) = \lambda e^{-\lambda x}, \; x > 0 \f]
  real function dexp(x,lambda) 
    implicit none
    real, intent(IN) :: x,lambda
    if (x<0) then
      dexp=0.
    else
      dexp=lambda*exp(-lambda*x)
    end if
  end function
  !**********************************************************************
  !> PDF for \f$ \chi^2 \f$ distribution
  real function dchi(x,k) 
    implicit none
    real, intent(IN) :: x, k
    real :: u
    if (x<0) then
      dchi=0.
    else
      dchi=(1./(gamma(k/2.)*2.**(k/2.)))*(x**((k/2.)-1.))*exp(-x/2.)
    end if
  end function
  !**********************************************************************
  !>PDF for Laplace distribution
  real function dlaplace(x,mu,b) 
    implicit none
    real, intent(IN) :: x, mu, b
    dlaplace=(1./(2*b))*exp(-abs(x-mu)/b)
  end function
end module
