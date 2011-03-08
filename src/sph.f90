!>SPH smooth paricle hydrodynamics.
module sph
  use cdata
  use general
  use output
  use forcing
  contains
  !************************************************************
  !>setup the particles in as set by initp
  subroutine setup_SPH
    implicit none
    integer :: i
    real :: rs, rthet, rphi !for sphere intial conditions
    !check that the mass of the particles has been set in run.in
    if (SPH_mass<epsilon(0.)) call fatal_error('setup_SPH',&
                              'mass of SPH particles not set')
    allocate(s(SPH_count))
    write(*,'(a)') ' --------------------SPH--------------------' 
    write(*,'(a,i4.1,a)') ' initialising ', SPH_count, ', particles '
    write(*,'(a,a,a,e14.7)') ' initial setup ', trim(SPH_init), ', with mass ', SPH_mass
    write(*,'(a,f10.5)') ' adiabatic index, gamma= ', SPH_gamma
    write(*,*) ' forcing will be applied to particles is selected'
    select case(SPH_init)
      case('sphere')
        !particles in random positions
        do i=1, SPH_count
          call random_number(rs)
          call random_number(rthet) ; call random_number(rphi)
          rs=rs*box_size/10. ; rthet=rthet*pi ; rphi=rphi*2*pi
          s(i)%x(1)=rs*sin(rthet)*cos(rphi)
          s(i)%x(2)=rs*sin(rthet)*sin(rphi)
          s(i)%x(3)=rs*cos(rthet)
          s(i)%u=0. ; s(i)%u1=0. ; s(i)%u2=0. !0 the velocity fields
          s(i)%a=0. ; s(i)%a1=0. ; s(i)%a2=0. !0 the acc fields
        end do
      case('random')
        !particles in random positions
        do i=1, SPH_count
          call random_number(s(i)%x(1))
          call random_number(s(i)%x(2)) ; call random_number(s(i)%x(3))
          s(i)%x(1)=box_size*s(i)%x(1)-box_size/2.
          s(i)%x(2)=box_size*s(i)%x(2)-box_size/2.
          s(i)%x(3)=box_size*s(i)%x(3)-box_size/2.
          s(i)%u=0. ; s(i)%u1=0. ; s(i)%u2=0. !0 the velocity fields
          s(i)%a=0. ; s(i)%a1=0. ; s(i)%a2=0. !0 the acc fields
        end do
      case('plane')
        !particles in planes at \pm z
        do i=1, SPH_count
          call random_number(s(i)%x(1)) ; call random_number(s(i)%x(2))
          s(i)%x(1)=box_size*s(i)%x(1)-box_size/2.
          s(i)%x(2)=box_size*s(i)%x(2)-box_size/2.
          if (i<floor(real(SPH_count)/2.)) then
            s(i)%x(3)=-0.48*box_size
          else
            s(i)%x(3)=0.48*box_size
          end if
          s(i)%u=0. ; s(i)%u1=0. ; s(i)%u2=0. !0 the velocity fields
          s(i)%a=0. ; s(i)%a1=0. ; s(i)%a2=0. !0 the acc fields
        end do
      case('corner')
        !particles in corners
        do i=1, SPH_count
          call random_number(s(i)%x(1))
          call random_number(s(i)%x(2)) ; call random_number(s(i)%x(3))
          s(i)%x(1)=0.1*(box_size*s(i)%x(1)-box_size/2.)
          s(i)%x(2)=0.1*(box_size*s(i)%x(2)-box_size/2.)
          s(i)%x(3)=0.1*(box_size*s(i)%x(3)-box_size/2.)
          if (i<floor(real(SPH_count)/2.)) then
            s(i)%x=s(i)%x+0.4*box_size
          else
            s(i)%x=s(i)%x-0.4*box_size
          end if
          s(i)%u=0. ; s(i)%u1=0. ; s(i)%u2=0. !0 the velocity fields
          s(i)%a=0. ; s(i)%a1=0. ; s(i)%a2=0. !0 the acc fields
        end do
      case default
        call fatal_error('setup_SPH','SPH_init not set to available value')
    end select
    s(:)%m=SPH_mass !set the masses of the particles
    !make an initial guess for the smoothing length
    s(:)%h=box_size/10.
    !now improve this guess
    call sph_get_h !sph.mod
  end subroutine
  !************************************************************
  !>evolve the SPH particles
  subroutine SPH_evolution
    implicit none
    integer :: i,j !for looping
    !firstly evolve h  due to motion of particles
    call sph_get_h !sph.mod
    !note this also sets rho and P at each particle
    !now get the acceleration at each particle
    call sph_get_acc
    do j=1,3
      !---------timestep the velocity--------------
      if (maxval(abs(s(i)%a1))==0) then
        s(:)%u(j)=s(:)%u(j)+dt*s(:)%a(j) !euler
      else if (maxval(abs(s(i)%a2))==0) then
        !first order adams-bashforth
        s(:)%u(j)=s(:)%u(j)+three_twos*dt*s(:)%a(j)-one_half*dt*s(:)%a1(j)
      else
        !second order adams-bashforth
        s(:)%u(j)=s(:)%u(j)+twenty_three_twelve*dt*s(:)%a(j)-four_thirds*dt*s(:)%a1(j)+five_twelths*dt*s(:)%a2(j)
      end if
      !---------timestep the position--------------
      if (maxval(abs(s(i)%u1))==0) then
        s(:)%x(j)=s(:)%x(j)+dt*s(:)%u(j) !euler
      else if (maxval(abs(s(i)%u2))==0) then
        !first order adams-bashforth
        s(:)%x(j)=s(:)%x(j)+three_twos*dt*s(:)%u(j)-one_half*dt*s(:)%u1(j)
      else
        !second order adams-bashforth
        s(:)%x(j)=s(:)%x(j)+twenty_three_twelve*dt*s(:)%u(j)-four_thirds*dt*s(:)%u1(j)+five_twelths*dt*s(:)%u2(j)
      end if
    end do
    s(i)%a2(:)=s(i)%a1(:) !store old acceleration 
    s(i)%a1(:)=s(i)%a(:)
    s(i)%u2(:)=s(i)%u1(:) !store old velocities 
    s(i)%u1(:)=s(i)%u(:)
    !----------------enforce periodicity----------------------
    if (periodic_bc) then
      do i=1,SPH_count
        !-------------x------------------     
        if (s(i)%x(1)>(box_size/2.)) then
          s(i)%x(1)=s(i)%x(1)-box_size
        else if (s(i)%x(1)<(-box_size/2.)) then
          s(i)%x(1)=s(i)%x(1)+box_size
        end if
        !-------------y------------------    write(*,'(a,a,a,e14.7)') ' initial setup ', trim(SPH_init), ' with mass ', SPH_mass
        if (s(i)%x(2)>(box_size/2.)) then
          s(i)%x(2)=s(i)%x(2)-box_size
        else if (s(i)%x(2)<(-box_size/2.)) then
          s(i)%x(2)=s(i)%x(2)+box_size
        end if
        !-------------z------------------
        if (s(i)%x(3)>(box_size/2.)) then
          s(i)%x(3)=s(i)%x(3)-box_size
        else if (s(i)%x(3)<(-box_size/2.)) then
          s(i)%x(3)=s(i)%x(3)+box_size
        end if
      end do
      !--------------------------------
    end if
    !finally check to see if we print
    if (mod(itime,shots)==0) then
      call diagnostics_SPH !sph.mod
      call print_SPH(itime/shots) !output.mod
    end if
  end subroutine
  !************************************************************
  !>get the acceleration of a particle
  subroutine sph_get_acc
    implicit none
    integer :: i,j
    integer :: peri, perj, perk !used to loop in periodic cases
    real, dimension(3) :: shift, force
    do i=1, SPH_count
      s(i)%a=0.
      do j=1, SPH_count
         if (i==j) cycle
         s(i)%a=s(i)%a-s(j)%m*(s(i)%P/(s(i)%rho**2))*sph_grad_W(s(i)%x,s(j)%x,s(i)%h)
         s(i)%a=s(i)%a-s(j)%m*(s(j)%P/(s(j)%rho**2))*sph_grad_W(s(i)%x,s(j)%x,s(j)%h)
         !periodic boundaries
         if (periodic_bc) then 
           do peri=-1,1 ; do perj=-1,1 ; do perk=-1,1
             if (peri==0.and.perj==0.and.perk==0) cycle
             shift=(/peri*box_size,perj*box_size,perk*box_size/)
             s(i)%a=s(i)%a-s(j)%m*(s(i)%P/(s(i)%rho**2))*sph_grad_W(s(i)%x,s(j)%x+shift,s(i)%h)
             s(i)%a=s(i)%a-s(j)%m*(s(j)%P/(s(j)%rho**2))*sph_grad_W(s(i)%x,s(j)%x+shift,s(j)%h)
           end do ; end do ;end do
        end if
      end do
      call get_forcing_gen(s(i)%x,force)
      s(i)%a=s(i)%a+force
    end do
  end subroutine
  !************************************************************
  !>get the correct smoothing length for each particle
  subroutine sph_get_h
    implicit none
    integer :: i
    do i=1, SPH_count
      s(i)%rho=0.
      call sph_get_rho(i)
      s(i)%h=(s(i)%rho/s(i)%m)**(-1./3.)
      call sph_get_rho(i)
      s(i)%h=(s(i)%rho/s(i)%m)**(-1./3.)
      call sph_get_rho(i)
      call sph_get_P(i)
    end do
  end subroutine
  !************************************************************
  !>get the density of the particle i based on the current
  !>smoothing length
  subroutine sph_get_rho(i)
    implicit none
    integer, intent(IN) :: i
    integer :: peri, perj, perk !used to loop in periodic cases
    real, dimension(3) :: shift
    real :: dist !distance between particles 
    integer :: j !for looping
    s(i)%rho=0.
    do j=1, sph_count
      dist=dist_gen(s(i)%x,s(j)%x)
      s(i)%rho=s(i)%rho+s(j)%m*sph_W(dist,s(i)%h)
      !periodic boundaries
      if (periodic_bc) then 
        do peri=-1,1 ; do perj=-1,1 ; do perk=-1,1
          if (peri==0.and.perj==0.and.perk==0) cycle
          shift=(/peri*box_size,perj*box_size,perk*box_size/)
          dist=dist_gen(s(i)%x,s(j)%x+shift)
          s(i)%rho=s(i)%rho+s(j)%m*sph_W(dist,s(i)%h)
        end do ; end do ;end do
      end if
    end do
  end subroutine 
  !************************************************************
  !>get the pressure at i using current density uses SPH_gamma from
  !>run.in as the adiabatic index (the default is 5/3)
  subroutine sph_get_P(i)
    implicit none
    integer, intent(IN) :: i
    s(i)%P=s(i)%rho**SPH_gamma
  end subroutine 
  !************************************************************
  !>get the derivative of the density of the particle i wrt h
  subroutine sph_get_drhodh(i)
    implicit none
    integer, intent(IN) :: i
    real :: dist !distance between particles 
    integer :: j !for looping
    s(i)%drhodh=0.
    do j=1, sph_count
      dist=dist_gen(s(i)%x,s(j)%x)
      s(i)%drhodh=s(j)%m*sph_dWdh(dist,s(i)%h)
    end do
  end subroutine  
  !************************************************************
  !>the smoothing kernal a cubic spline is used Springel 2010
  real function sph_W(r,h)
    implicit none
    real, intent(IN) :: r, h
    real :: q !ration of distance between particles and smoothing length
    q=r/(2*h)
    if (q<0) then
      sph_W=(8/pi)*(1-6*(q**2-q**3))
    else if (q<1) then
      sph_W=(8/pi)*2*(1-q)**3
    else 
      sph_W=0.
    end if
  end function
  !************************************************************
  !>derivative of the smoothing kernal wrt h
  real function sph_dWdh(r,h)
    implicit none
    real, intent(IN) :: r, h
    real :: q !ratio of distance between particles and smoothing length
    q=r/(2*h)
    if (q<0) then
      sph_dWdh=(r/(2*h**2))*(8/pi)*(12*q-18*q**2)
    else if (q<1) then
      sph_dWdh=(r/(2*h**2))*(8/pi)*6*(1-q)**2
    else 
      sph_dWdh=0.
    end if
  end function
  !************************************************************
  !>gradient of the smoothing kernal wrt h
  function sph_grad_W(xi,xj,h)
    implicit none
    real, dimension(3) :: sph_grad_W
    real, intent(IN) :: xi(3), xj(3) !particle positions
    real, intent(IN) :: h
    real :: r, q !ratio of distance between particles and smoothing length
    r=dist_gen(xi,xj)
    q=r/(2*h)
    if (q<0) then
      sph_grad_W=(xi-xj)*(8/pi)*(-12*q+18*q**2)/(2*r*h)
    else if (q<1) then
      sph_grad_W=-(xi-xj)*(8/pi)*6*((1-q)**2)/(2*r*h)
    else
      sph_grad_W=0.
    end if
  end function
  !************************************************************
  !>diagnostics for SPH code
  subroutine diagnostics_SPH
    implicit none
    real :: uinfo(SPH_count,2)
    real :: SPH_maxu, SPH_maxdu, SPH_urms
    integer :: i
    uinfo(:,1)=sqrt(s(:)%u(1)**2+s(:)%u(2)**2+s(:)%u(3)**2)
    uinfo(:,2)=sqrt(s(:)%a(1)**2+s(:)%a(2)**2+s(:)%a(3)**2)
    !in here we need to determine
    SPH_maxu=maxval(uinfo(:,1)) ; SPH_maxdu=maxval(uinfo(:,2))
    SPH_urms=sqrt(sum(uinfo(:,1)**2)/SPH_count)
    open(unit=78,file='data/SPH_ts.log',position='append')
    if (itime==shots) then
      write(78,*) '%-var--t----------maxu---------maxdu---------urms-----'
      if (particles_only) write(*,*) '%-var--t--------maxu------maxdu-------urms-----'
    end if
    write(78,'(i5.3,f6.2,f15.8,f15.8,f15.8)') &
    itime/shots,t,SPH_maxu,SPH_maxdu,SPH_urms
    if (particles_only) then
      write(*,'(i5.3,f6.2,f12.8,f13.5,f12.8)') &
      itime/shots,t,SPH_maxu,SPH_maxdu,SPH_urms
    end if
    close(78)
  end subroutine
end module
!>\page SPH 
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
