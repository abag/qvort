!>all the SPH kernels needed are contained within this module
module sph_kernel
   use cdata
   use general
  contains
  !************************************************************
  !>the smoothing kernal M4 cubic spline
  real function sph_W(r,h)
    implicit none
    real, intent(IN) :: r, h
    real :: q !ration of distance between particles and smoothing length
    q=r/h
    if (q<1) then
      sph_W=1-1.5*q**2+0.75*q**3
    else if (q<2) then
      sph_W=0.25*(2-q)**3
    else 
      sph_W=0.
    end if
    sph_W=sph_W*(1/(pi*h**3))
  end function
  !************************************************************
  !>derivative of the smoothing kernal wrt h
  real function sph_dWdh(r,h)
    implicit none
    real, intent(IN) :: r, h
    real :: q !ratio of distance between particles and smoothing length
    q=r/h
    if (q<1) then
      sph_dWdh=-3+7.5*q**2-4.5*q**3
    else if (q<2) then
      sph_dWdh=-6+12*q-7.5*q**2+1.5*q**3
    else 
      sph_dWdh=0.
    end if
    sph_dWdh=sph_dWdh*(1/(pi*h**4))
  end function
  !************************************************************
  !>gradient of the smoothing kernal
  function sph_grad_W(xi,xj,h)
    implicit none
    real, dimension(3) :: sph_grad_W
    real, intent(IN) :: xi(3), xj(3) !particle positions
    real, intent(IN) :: h
    real :: r, q !ratio of distance between particles and smoothing length
    logical, parameter :: no_clumping=.false.
    r=dist_gen(xi,xj)
    q=r/h
    if (no_clumping) then
      if (q<2/3) then
        sph_grad_W=(xi-xj)
      else if (q<1) then
        sph_grad_W=(xi-xj)*(3*q-2.25*q**2)
      else if (q<2) then
        sph_grad_W=-(xi-xj)*0.75*(2-q)**2
      else
        sph_grad_W=0.
      end if
    else
      if (q<1) then
        sph_grad_W=(xi-xj)*(3*q-2.25*q**2)
      else if (q<2) then
        sph_grad_W=(xi-xj)*0.75*(2-q)**2
      else
        sph_grad_W=0.
      end if
    end if 
    sph_grad_W=sph_grad_W*(-1/(pi*r*h**4))
  end function
  !************************************************************
  !>gravitiational kernal
  function sph_phi(r,h)
    implicit none
    real, dimension(3) :: sph_phi
    real, intent(IN) :: r, h
    real :: q !ratio of distance between particles and smoothing length
    q=r/h
    if (q<1) then
      sph_phi=(4./3.)*q-1.2*q**3+0.5*q**4
    else if (q<2) then
      sph_phi=(8./3.)*q-3*q**2+1.2*q**3-(q**4)/6.-(q**(-2))/15.
    else
      sph_phi=q**(-2)
    end if
    sph_phi=sph_phi/(h**2)
  end function
end module
!>\page SPH Smoothed Particle Hydrodynamics
!>The M4 cubic spline kernel (Monaghan \& Lattanzio 1985) is used in many implementations of SPH, due to its simple form and its compact support.  The M4 kernel is a function of \f$s\equiv r/h\f$ only. For \f$D=1,\;2,\;{\rm and}\;3\,\f$ dimensions, it takes the form
!>\f{eqnarray}{
!>W (s) &=& \frac{\sigma_{\!_D}}{h^D} 
!>\left\{ \;\; \begin{array}{ll} 
!>1 - \frac{3}{2}s^2 + \frac{3}{4}s^3 \;\;\;\; & \;\;0 \leq s \leq 
!>1\,; \\ \frac{1}{4}(2-s)^3 \;\;\;\; & \;\;1 \leq s \leq 2\,; \\ 
!>0 \;\;\;\; & \;\;s > 2 \,.
!>\end{array} \right . 
!>\f}
!>where \f$\sigma_{_1}=2/3\f$, \f$\sigma_{_2}=10/7\pi\f$, and \f$\sigma_{_3}=1/\pi\f$. The first spatial derivative is
!>\f{eqnarray}{
!>\frac{dW}{dr}(s) &=& - \frac{\sigma_{\!_D}}{h^{D+1}}
!>\left\{ \;\; \begin{array}{ll}
!>3s - \frac{9}{4}s^2 \;\;\;\; & \;\;0 \leq s \leq 1\,; \\ 
!>\frac{3}{4}(2-s)^2 \;\;\;\; & \;\;1 \leq s \leq 2\,; \\ 
!>0 \;\;\;\; & \;\;s > 2 \,.
!>\end{array} \right .
!>\f}
!>We should also investigate using the modified derivative proposed by Thomas \& Couchman (1992) to prevent the clumping instability.
!>For `grad-h' SPH, the \f$\Omega\f$ correction kernel function is given by 
!>\f{eqnarray}{
!>\frac{\partial W}{\partial h}(s) = \frac{\sigma_{\!_D}}{h^{D+1}} 
!>\left\{ \;\; \begin{array}{ll} 
!>-D + \frac{3}{2}(D + 2)\,s^2 - \frac{3}{4}(D + 3)\,s^3 
!>\;\;\;\; & \;\;0 \leq s \leq 1\,; \\ 
!>-2\,D + 3\,(D + 1)\,s - \frac{3}{2}(D + 2)\,s^2 + \frac{1}{4}(D + 3)\,s^3 
!> \;\;\;\; & \;\;1 \leq s \leq 2\,; \\ 
!>0 \;\;\;\; & \;\;s > 2 \,.
!>\end{array} \right . 
!>\f}
!> the smoothing length is set to:
!>\f[ h_i  = \eta_{_{\rm SPH}} \left( \frac{m_i }{\rho_i}
!>\right)^{\frac{1}{D}}\ \f]
!>where \f$m_i\f$ is the mass of particle \f$i\f$, \f$\rho_i\f$ is the SPH density at the position of particle \f$i\f$, 
!>\f$D\f$ is the spatial dimensionality, and \f$\eta_{_{\rm SPH}}\f$ is a parameter that controls the mean number of neighbours, 
!>\f$\bar{N}_{_{\rm NEIB}} \simeq 2\,{\cal R}\,\eta_{_{\rm SPH}} \;,\;\; \pi\,({\cal R}\eta_{_{\rm SPH}})^2 \;,\;\;(4\pi /3)({\cal R}\eta_{_{\rm SPH}})^3\,\,\f$ in one, two and three dimensions respectively. \f${\cal R}=2\f$
!>for the kernal choice detailed above. \f$\;\rho_i\f$ is calculated using 
!>\f[ \rho_i = \sum \limits_{j=1}^{N}  m_j W({\bf r}_{ij},h_i), \f]
!>where \f${\bf r}_{ij} \equiv {\bf r}_i - {\bf r}_j\f$, and the summation includes particle \f$i\f$ itself. 
!>Since the smoothing length is needed in order to calculate the density and vice-versa, \f$h_i\f$ and \f$\rho_i\f$ are obtained by iteration.
!>Once \f$h\f$ and \f$\rho\f$ are evaluated for all particles, the terms in the SPH fluid equations can be computed. The momentum equation is 
!>\f{eqnarray}{
!>\frac{d{\bf v}_i }{dt} &=& -
!>\sum \limits_{j=1}^{N}  m_j  \left( 
!>\frac{P_i}{\Omega_i \rho_i^2} \nabla_i W({\bf r}_{ij} ,h_i) + 
!>\frac{P_j}{\Omega_j \rho_j^2} \nabla_i W({\bf r}_{ij} ,h_j) \right)
!>\,,
!>\f}
!>where \f$P_i\f$ is the pressure of particle \f$i\f$, \f$\nabla_i W\f$ is the gradient of the kernel function at the position of particle \f$i\f$, and 
!>\f{eqnarray}
!>\Omega_i  &=& 1 - \frac{\partial h_i }{\partial \rho_i } 
!>\sum \limits_{j=1}^{N}  m_j  \frac{\partial W}{\partial h}
!>({\bf r}_{ij} , h_i )\,.
!>\f}
!>\f$\Omega_i\f$ is a dimensionless quantity that corrects for the spatial variability of \f$h\f$. \f$\partial h_i / \partial \rho_i\f$ is obtained explicitly. \f$\partial W / \partial h\f$ is obtained from the kernel function.
!>We note that the gradient of the kernal function can be calculated in the following manner:
!>\f[ \nabla_i W=\frac{\mathbf{x}_i-\mathbf{x}_j}{{\bf r}_{ij}} 
!>\frac{\partial W}{\partial {\bf r}_{ij}} \f]
!>and that the particle pressure can be determined from the equation of state
