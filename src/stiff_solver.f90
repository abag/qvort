module stiff_solver
  !CONTAINS ALL THE ROUTINES FOR SOLVING THE STIFF ODE SOLVER REQIURED BY QUASI PARTICLES
  use cdata
  use hamiltonian
  real, allocatable, private :: eta(:) !used for adjusting timestep
  real, dimension(6,7), private :: BDF_coeff
  contains
  subroutine set_BDF_coeff()
    implicit none
    !SET THE BDF COEFFICIENTS ACCORDING TO:http://www.scholarpedia.org/article/Backward_differentiation_formulas
    BDF_coeff(1,:)=(/1., -1., 0., 0., 0., 0.,0./)
    BDF_coeff(2,:)=(/3./2., -2., 1./2., 0., 0., 0.,0./)
    BDF_coeff(3,:)=(/11./6., -3., 3./2., -1./3., 0., 0.,0./)
    BDF_coeff(4,:)=(/25./12., -4., 3., -4./3., 1./4., 0.,0./)
    BDF_coeff(5,:)=(/137./60., -5., 5., -10./3., 5./4., -1./5.,0./)
    BDF_coeff(6,:)=(/49./20., -6., 15./2., -20./3., 15./4., -6./5.,1./6./)
    write(*,*) 'set up backwards difference coefficients for stiff ode solver'
    allocate(eta(quasi_pcount))
    write(*,*) 'allocated adaptive timestep array'
  end subroutine
  !******************************************************************************************
  subroutine BDF(i,order)
    implicit none
    integer, intent(IN) :: order !the requested order of the BDF
    integer, intent(IN) :: i !the particle index
    real :: pdot(3), rdot(3) !momentum and position time derivatives
    integer :: useorder
    integer :: j !used for looping
    !first check if we can achieve the order asked for
    useorder=min(order,itime) !take the minimum of the order and the itime loop
    !get the momentum/position derivatives
    call velocity_quasip(i,rdot,pdot) !hamiltonian.mod
    do j=1,3
      g(i)%x(j)=(dt*rdot(j)-dot_product(BDF_coeff(useorder,2:7),g(i)%xold(1:6,j)))/BDF_coeff(useorder,1)
      g(i)%p(j)=(dt*pdot(j)-dot_product(BDF_coeff(useorder,2:7),g(i)%pold(1:6,j)))/BDF_coeff(useorder,1)
    end do
    !now store all the old values
    do j=6,2,-1
      g(i)%xold(j,:)=g(i)%xold(j-1,:)
      g(i)%pold(j,:)=g(i)%pold(j-1,:)
    end do
    !finally the most recent
    g(i)%xold(1,:)=g(i)%x ; g(i)%pold(1,:)=g(i)%p
    !now we want to alter the timestep
    if (useorder>2) then !can only do once we have a few values stored
      eta(i)=vector_angle(g(i)%pold(1,:),g(i)%pold(2,:)) !general.mod
    else
      eta(i)=0.
    end if
  end subroutine
  !*********************************************************************
  subroutine BDF_dt_adjust
    implicit none
    real :: dt_min=1E-17
    real :: dt_max=1E-9
    open(unit=78,file='data/qp_eta.log',position='append')
      write(78,*) t, maxval(eta), dt
    close(78)
    if (maxval(eta)>.00001) dt=dt/4.
    if (maxval(eta)<.000005) dt=1.2*dt
    if (dt>dt_max) dt=dt_max
    if (dt<dt_min) dt=dt_min
  end subroutine
endmodule
