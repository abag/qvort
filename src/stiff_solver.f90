!>all routines to solve stiff ode problems
!>in particular this is used for reflection of quasi particles
module stiff_solver
  use cdata
  use hamiltonian
  real, allocatable, private :: eta(:) !used for adjusting timestep
  real, dimension(6,7), private :: BDF_coeff !array of backwards diffference coefficients
  real, private :: t_int, dt_int ! internal time,timestep 
  contains
  !>set the coefficients for the backwards difference scheme, used only
  !>if solving quasi-particles.
  subroutine set_BDF_coeff()
    implicit none
    BDF_coeff(1,:)=(/1., -1., 0., 0., 0., 0.,0./)
    BDF_coeff(2,:)=(/3./2., -2., 1./2., 0., 0., 0.,0./)
    BDF_coeff(3,:)=(/11./6., -3., 3./2., -1./3., 0., 0.,0./)
    BDF_coeff(4,:)=(/25./12., -4., 3., -4./3., 1./4., 0.,0./)
    BDF_coeff(5,:)=(/137./60., -5., 5., -10./3., 5./4., -1./5.,0./)
    BDF_coeff(6,:)=(/49./20., -6., 15./2., -20./3., 15./4., -6./5.,1./6./)
    write(*,*) 'set up backwards difference coefficients for stiff ode solver'
    allocate(eta(quasi_pcount)) ; dt_int=dt !set internal timestep
    write(*,*) 'allocated adaptive timestep array'
  end subroutine
  !******************************************************************************************
  !>backwards difference method, usefule for stiff ODE's for a 
  !>detailed description of the method see
  !>http://www.scholarpedia.org/article/Backward_differentiation_formulas
  !>order is set in quasi_timestep.
  subroutine BDF(i,order)
    implicit none
    integer, intent(IN) :: order !the requested order of the BDF
    integer, intent(IN) :: i !the particle index
    real :: pdot(3), rdot(3) !momentum and position time derivatives
    real :: pdot_fut(3), rdot_fut(3) !future momentum and position time derivatives
    real :: predict_x(3), predict_p(3) !future prediction of p used to adjust timestep
    integer :: useorder !the order of the backwards difference scheme
    integer :: j !used for looping
    !t_int=t !intially set the internal time to be the global time
    !do while (t_int<t+dt) !loop until internal time is 'ready'
      !first check if we can achieve the order asked for
      !THIS NEEDS TO BE IMPROVED TO TAKE ACCOUNT FOR INTERNAL TIME LOOP
      useorder=min(order,itime) !take the minimum of the order and the itime loop
      !get the momentum/position derivatives
      call velocity_quasip(i,rdot,pdot) !hamiltonian.mod
      g(i)%rdot=rdot ; g(i)%pdot=pdot !store for printing to ts
      !predict the future value of p
      do j=1,3
        predict_p(j)=(dt_int*pdot(j)-dot_product(BDF_coeff(useorder,2:7),g(i)%pold(1:6,j)))/BDF_coeff(useorder,1)
        predict_x(j)=(dt_int*rdot(j)-dot_product(BDF_coeff(useorder,2:7),g(i)%xold(1:6,j)))/BDF_coeff(useorder,1)
      end do
      call velocity_quasip_gen(predict_x,predict_p,rdot_fut,pdot_fut) !hamiltonian.mod
      rdot=0.5*(rdot+rdot_fut) ; pdot=0.5*(pdot+pdot_fut)
      !determine the angle between p and predict_p
      if (useorder>2) then !can only do once we have a few values stored
        eta(i)=vector_angle(g(i)%pold(1,:),predict_p) !general.mod
      else
        eta(i)=0.
      end if
      !adjust the timestep based on this angle
      call BDF_dt_adjust
      !timestep r and p
      do j=1,3
        g(i)%x(j)=(dt_int*rdot(j)-dot_product(BDF_coeff(useorder,2:7),g(i)%xold(1:6,j)))/BDF_coeff(useorder,1)
        g(i)%p(j)=(dt_int*pdot(j)-dot_product(BDF_coeff(useorder,2:7),g(i)%pold(1:6,j)))/BDF_coeff(useorder,1)
      end do
      !now store all the old values
      do j=6,2,-1
        g(i)%xold(j,:)=g(i)%xold(j-1,:)
        g(i)%pold(j,:)=g(i)%pold(j-1,:)
      end do
      !finally the most recent
      g(i)%xold(1,:)=g(i)%x ; g(i)%pold(1,:)=g(i)%p
      t_int=t_int+dt_int !increment internal timestep
    !end do
  end subroutine
  !*********************************************************************
  !>adjust the timestep based on the angle between future momentum vector
  !>of quasi particle and current momentum vector
  subroutine BDF_dt_adjust
    implicit none
    real :: dt_min=1E-20
    real :: dt_max=1E-9
    open(unit=78,file='data/qp_eta.log',position='append')
      write(78,*) t_int, maxval(eta), dt_int
    close(78)
    if (maxval(eta)>.0000001) dt_int=dt_int/2.
    if (maxval(eta)<.00000005) dt_int=2*dt_int
    if (dt_int>dt_max) dt_int=dt_max
    if (dt_int<dt_min) dt_int=dt_min
  end subroutine
endmodule
