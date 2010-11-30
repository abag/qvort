module hamiltonian
  use cdata
  use general
  use timestep
  use output
  use tree
  !*******CONSTANTS ASSOCIATED WITH QUASI-PARTICLES****************
  real, private, parameter :: mbare=5.017e-24
  real, private, parameter :: mfermi=1.51e-23 
  real, parameter :: pfermi=8.281e-20 !not private as needed for initial conditions
  real, private, parameter :: efermi=2.27e-16
  real, private, parameter :: vfermi=5.48e3
  real, private, parameter :: delta_gap=2.43e-19
  real, private, parameter :: unitdt=2.9e-6 !unit time for matlab  ode solver 
  contains
  !***************************************************************************************
  subroutine velocity_quasip(i,rdot,pdot)
    !for each particle we 
    ! 1. define the hamiltonian
       !in order to do this we require the velocity (vortex) derivatives in all 3 dimension (7 evaluations)
    ! 2. intially update the particles using an euler step
    ! 3. build up to solving the stiff problem 
    implicit none
    integer, intent(IN) :: i
    real, intent(OUT) :: rdot(3), pdot(3)
    real :: u_sup(3), du_sup_x(3), du_sup_y(3), du_sup_z(3) !superfluid velocity and derivatives
    real :: u_plusx(3), u_minusx(3) !superfluid velocity at plus minus x
    real :: u_plusy(3), u_minusy(3) !superfluid velocity at plus minus y
    real :: u_plusz(3), u_minusz(3) !superfluid velocity at plus minus z
    real :: epsilonp !used in the hamiltonian
    real :: dist, min_dist !need the minimum distance between the particle and the vortex 
    integer :: peri, perj, perk !used to loop in periodic cases
    integer :: j !used for loops
    !Get the superfluid velocity at the position of the particle
    !account for all possible velocity fields
    u_sup=0. !must be zeroed for intially
    u_plusx=0. ; u_minusx=0. ; u_plusy=0. ; u_minusy=0.
    u_plusz=0. ; u_minusz=0.
    !first we need to find the minimum distance between the particle and the vortice's
    select case(velocity)
      case('LIA','BS')
        min_dist=100. !arbitrarily large
        do j=1,pcount !loop over all vortex points
          dist=dist_gen(g(i)%x,f(j)%x)
          if (dist<min_dist) min_dist=dist
        end do
      case('Tree')
        !we need a tree minimum distance algorithm
        call fatal_error('velocity_quasip','tree adaptation not ready yet')
    end select
    !now we check wether the particle is within the vortex core or outside
    if (min_dist<corea) then
      !empty core so 0 velocity
      u_sup=0. ; u_plusx=0. ; u_minusx=0.
      u_plusy=0. ; u_minusy=0.
      u_plusz=0. ; u_minusz=0.
    else
      select case(velocity)
        case('LIA','BS')
          call biot_savart_general(g(i)%x,u_sup) !timestep.mod
          call biot_savart_general(g(i)%x+(/delta,0.,0./),u_plusx) !timestep.mod
          call biot_savart_general(g(i)%x+(/-delta,0.,0./),u_minusx) !timestep.mod
          call biot_savart_general(g(i)%x+(/0.,delta,0./),u_plusy) !timestep.mod
          call biot_savart_general(g(i)%x+(/0.,-delta,0./),u_minusy) !timestep.mod
          call biot_savart_general(g(i)%x+(/0.,0.,delta/),u_plusz) !timestep.mod 
          call biot_savart_general(g(i)%x+(/0.,0.,-delta/),u_minusz) !timestep.mod
          if (periodic_bc) then
            !we must shift the mesh in all 3 directions, all 26 permutations needed!
            !needs adpating for periodic boundary conditons - see above
            do peri=-1,1 ; do perj=-1,1 ; do perk=-1,1
              if (peri==0.and.perj==0.and.perk==0) cycle
              call biot_savart_general_shift(g(i)%x,u_sup, &
                   (/peri*box_size,perj*box_size,perk*box_size/))!timestep.mod
            end do ; end do ;end do
          end if
        case('Tree')
          call fatal_error('velocity_quasip','tree adaptation not ready yet')
      end select
    end if
    !now we must find both rdot and pdot using the hamiltonian for the quasi-particle
    epsilonp=(dot_product(g(i)%p,g(i)%p)/(2.*mfermi))-efermi
    rdot(:)=(epsilonp/sqrt(epsilonp**2+delta_gap**2))*g(i)%p(:)/mfermi+u_sup
    !let us now calculate the derivative of the superfluid velocity field
    du_sup_x(:)=((u_plusx-u_minusx)/(2.*delta))
    du_sup_y(:)=((u_plusy-u_minusy)/(2.*delta))
    du_sup_z(:)=((u_plusz-u_minusz)/(2.*delta))
    !calculate pdot
    pdot(1)=-dot_product(g(i)%p,du_sup_x)
    pdot(2)=-dot_product(g(i)%p,du_sup_y)
    pdot(3)=-dot_product(g(i)%p,du_sup_z)
    !finally store the energy for output
    g(i)%energy=sqrt(epsilonp**2+delta_gap**2)+dot_product(g(i)%p,u_sup)
  end subroutine
  !******************************************************************************************************
  subroutine velocity_quasip_gen(pos,mom,rdot,pdot)
    !a generalised version of the above routine for use with a specified position 
    implicit none
    real, intent(IN) :: pos(3), mom(3)
    real, intent(OUT) :: rdot(3), pdot(3)
    real :: u_sup(3), du_sup_x(3), du_sup_y(3), du_sup_z(3) !superfluid velocity and derivatives
    real :: u_plusx(3), u_minusx(3) !superfluid velocity at plus minus x
    real :: u_plusy(3), u_minusy(3) !superfluid velocity at plus minus y
    real :: u_plusz(3), u_minusz(3) !superfluid velocity at plus minus z
    real :: epsilonp !used in the hamiltonian
    real :: dist, min_dist !need the minimum distance between the particle and the vortex 
    integer :: peri, perj, perk !used to loop in periodic cases
    integer :: j !used for loops
    !Get the superfluid velocity at the position of the particle
    !account for all possible velocity fields
    u_sup=0. !must be zeroed for intially
    u_plusx=0. ; u_minusx=0. ; u_plusy=0. ; u_minusy=0.
    u_plusz=0. ; u_minusz=0.
    !first we need to find the minimum distance between the particle and the vortice's
    select case(velocity)
      case('LIA','BS')
        min_dist=100. !arbitrarily large
        do j=1,pcount !loop over all vortex points
          dist=dist_gen(pos,f(j)%x)
          if (dist<min_dist) min_dist=dist
        end do
      case('Tree')
        !we need a tree minimum distance algorithm
        call fatal_error('velocity_quasip','tree adaptation not ready yet')
    end select
    !now we check wether the particle is within the vortex core or outside
    if (min_dist<corea) then
      !empty core so 0 velocity
      u_sup=0. ; u_plusx=0. ; u_minusx=0.
      u_plusy=0. ; u_minusy=0.
      u_plusz=0. ; u_minusz=0.
    else
      select case(velocity)
        case('LIA','BS')
          call biot_savart_general(pos,u_sup) !timestep.mod
          call biot_savart_general(pos+(/delta,0.,0./),u_plusx) !timestep.mod
          call biot_savart_general(pos+(/-delta,0.,0./),u_minusx) !timestep.mod
          call biot_savart_general(pos+(/0.,delta,0./),u_plusy) !timestep.mod
          call biot_savart_general(pos+(/0.,-delta,0./),u_minusy) !timestep.mod
          call biot_savart_general(pos+(/0.,0.,delta/),u_plusz) !timestep.mod 
          call biot_savart_general(pos+(/0.,0.,-delta/),u_minusz) !timestep.mod
          if (periodic_bc) then
            !we must shift the mesh in all 3 directions, all 26 permutations needed!
            !needs adpating for derivatives - see above
            do peri=-1,1 ; do perj=-1,1 ; do perk=-1,1
              if (peri==0.and.perj==0.and.perk==0) cycle
              call biot_savart_general_shift(pos,u_sup, &
                   (/peri*box_size,perj*box_size,perk*box_size/))!timestep.mod
            end do ; end do ;end do
          end if
        case('Tree')
          call fatal_error('velocity_quasip','tree adaptation not ready yet')
      end select
    end if
    !now we must find both rdot and pdot using the hamiltonian for the quasi-particle
    epsilonp=(dot_product(mom,mom)/(2.*mfermi))-efermi
    rdot(:)=(epsilonp/sqrt(epsilonp**2+delta_gap**2))*mom(:)/mfermi+u_sup
    !let us now calculate the derivative of the superfluid velocity field
    du_sup_x(:)=((u_plusx-u_minusx)/(2.*delta))
    du_sup_y(:)=((u_plusy-u_minusy)/(2.*delta))
    du_sup_z(:)=((u_plusz-u_minusz)/(2.*delta))
    !calculate pdot
    pdot(1)=-dot_product(mom,du_sup_x)
    pdot(2)=-dot_product(mom,du_sup_y)
    pdot(3)=-dot_product(mom,du_sup_z)
  end subroutine
  !*******************************************************************************************************
  subroutine velocity_quasip_jacobian(i,rdot,pdot,jac)
    !this routine gives not only rdot and pdot but also the full jacobian
    !ANDREW THIS IS NOT FINISHED!!!!!!!!!!!!!!!!!
    implicit none
    integer, intent(IN) :: i
    real, intent(OUT) :: rdot(3), pdot(3)
    real, intent(OUT) :: jac(6,6)
    real :: u_sup(3), du_sup_x(3), du_sup_y(3), du_sup_z(3) !superfluid velocity and derivatives
    real :: du_sup_x2(3), du_sup_y2(3), du_sup_z2(3) !superfluid velocity 2nd derivatives (1 var)
    real :: du_sup_xy(3), du_sup_yz(3), du_sup_zx(3) !superfluid velocity 2nd derivatives (2 var)
    !---------------------for partial derivatives------------------------
    real :: u_plusx(3), u_minusx(3) !superfluid velocity at plus minus x
    real :: u_plusy(3), u_minusy(3) !superfluid velocity at plus minus y
    real :: u_plusz(3), u_minusz(3) !superfluid velocity at plus minus z
    !----------------for partial derivatives of two variables------------
    real :: u_px_py(3), u_px_my(3), u_mx_py(3) , u_mx_my(3) !box permutations for pderiv
    real :: u_py_pz(3), u_py_mz(3), u_my_pz(3) , u_my_mz(3) !of two variables
    real :: u_pz_px(3), u_pz_mx(3), u_mz_px(3) , u_mz_mx(3)  
    !-------------------------------------------------------------------
    real :: epsilonp !used in the hamiltonian
    real :: dist, min_dist !need the minimum distance between the particle and the vortex 
    integer :: peri, perj, perk !used to loop in periodic cases
    integer :: j !used for loops
    !Get the superfluid velocity at the position of the particle
    !account for all possible velocity fields
    u_sup=0. !must be zeroed for intially
    u_plusx=0. ; u_minusx=0. ; u_plusy=0. ; u_minusy=0.
    u_plusz=0. ; u_minusz=0.
    u_px_py=0. ; u_px_my=0. ; u_mx_py=0. ; u_mx_my=0.  
    u_py_pz=0. ; u_py_mz=0. ; u_my_pz=0. ; u_my_mz=0.  
    u_pz_px=0. ; u_pz_mx=0. ; u_mz_px=0. ; u_mz_mx=0.  
    !first we need to find the minimum distance between the particle and the vortice's
    select case(velocity)
      case('LIA','BS')
        min_dist=100. !arbitrarily large
        do j=1,pcount !loop over all vortex points
          dist=dist_gen(g(i)%x,f(j)%x)
          if (dist<min_dist) min_dist=dist
        end do
      case('Tree')
        !we need a tree minimum distance algorithm
        call fatal_error('velocity_quasip','tree adaptation not ready yet')
    end select
    !now we check whether the particle is within the vortex core or outside
    if (min_dist<corea) then
      !empty core so 0 velocity - probably don't need these calls!
      u_sup=0. ; u_plusx=0. ; u_minusx=0.
      u_plusy=0. ; u_minusy=0.
      u_plusz=0. ; u_minusz=0.
    else
      select case(velocity)
        case('LIA','BS')
          !get the superfuld velocity at g(i)%x
          call biot_savart_general(g(i)%x,u_sup) !timestep.mod
          !these calls are for calculating partial derivatives (1st/2nd order) of 1 var.
          call biot_savart_general(g(i)%x+(/delta,0.,0./),u_plusx)
          call biot_savart_general(g(i)%x+(/-delta,0.,0./),u_minusx) 
          call biot_savart_general(g(i)%x+(/0.,delta,0./),u_plusy) 
          call biot_savart_general(g(i)%x+(/0.,-delta,0./),u_minusy) 
          call biot_savart_general(g(i)%x+(/0.,0.,delta/),u_plusz) 
          call biot_savart_general(g(i)%x+(/0.,0.,-delta/),u_minusz) 
          !these calls are for calculating 2nd order partial derivatives of 2 var.
          call biot_savart_general(g(i)%x+(/delta,delta,0./),u_px_py) !for d^2/dxdy
          call biot_savart_general(g(i)%x+(/delta,-delta,0./),u_px_my)
          call biot_savart_general(g(i)%x+(/-delta,delta,0./),u_mx_py)
          call biot_savart_general(g(i)%x+(/-delta,-delta,0./),u_mx_my)
          call biot_savart_general(g(i)%x+(/0.,delta,delta/),u_py_pz) ! for d^2/dydz
          call biot_savart_general(g(i)%x+(/0.,delta,-delta/),u_py_mz)
          call biot_savart_general(g(i)%x+(/0.,-delta,delta/),u_my_pz)
          call biot_savart_general(g(i)%x+(/0.,-delta,-delta/),u_my_mz)
          call biot_savart_general(g(i)%x+(/delta,0.,delta/),u_pz_px)! for d^2/dzdx
          call biot_savart_general(g(i)%x+(/-delta,0.,delta/),u_pz_mx)
          call biot_savart_general(g(i)%x+(/delta,0.,-delta/),u_mz_px)
          call biot_savart_general(g(i)%x+(/-delta,0.,-delta/),u_mz_mx)
          !if (periodic_bc) then
            !we must shift the mesh in all 3 directions, all 26 permutations needed!
            !THIS IS NOT READY - ALL THE PARTIAL DERIVATIVES NEED TO BE PUT IN HERE
          !  do peri=-1,1 ; do perj=-1,1 ; do perk=-1,1
          !    if (peri==0.and.perj==0.and.perk==0) cycle
          !    call biot_savart_general_shift(g(i)%x,u_sup, &
          !         (/peri*box_size,perj*box_size,perk*box_size/))!timestep.mod
          !  end do ; end do ;end do
          !end if
        case('Tree')
          call fatal_error('velocity_quasip','tree adaptation not ready yet')
      end select
    end if
    !now we must find both rdot and pdot using the hamiltonian for the quasi-particle
    epsilonp=(dot_product(g(i)%p,g(i)%p)/(2.*mfermi))-efermi
    rdot(:)=(epsilonp/sqrt(epsilonp**2+delta_gap**2))*g(i)%p(:)/mfermi+u_sup
    !let us now calculate the derivative of the superfluid velocity field
    du_sup_x(:)=((u_plusx-u_minusx)/(2.*delta))
    du_sup_y(:)=((u_plusy-u_minusy)/(2.*delta))
    du_sup_z(:)=((u_plusz-u_minusz)/(2.*delta))
    !and second derivatives of the superfluid velocity field
    du_sup_x2(:)=((u_plusx-2.*u_sup+u_minusx)/(delta**2))
    du_sup_y2(:)=((u_plusy-2.*u_sup+u_minusy)/(delta**2))
    du_sup_z2(:)=((u_plusz-2.*u_sup+u_minusz)/(delta**2))
    !finally second derivatives of two variables
    du_sup_xy(:)=((u_px_py+u_mx_my-u_px_my-u_mx_py)/(4.*delta**2))
    du_sup_yz(:)=((u_py_pz+u_my_mz-u_py_mz-u_my_pz)/(4.*delta**2))
    du_sup_zx(:)=((u_pz_px+u_mz_mx-u_pz_mx-u_mz_px)/(4.*delta**2))
    !calculate pdot
    pdot(1)=-dot_product(g(i)%p,du_sup_x)
    pdot(2)=-dot_product(g(i)%p,du_sup_y)
    pdot(3)=-dot_product(g(i)%p,du_sup_z)
    !store the energy for output
    g(i)%energy=sqrt(epsilonp**2+delta_gap**2)+dot_product(g(i)%p,u_sup)
    !finally determine the jacobian matrix, this is a 6x6 matrix x=(p,r)
    !note we follow the fortran convention first index is column, second row
    !this is 4 blocks of 3x3 matrices
    !d(pdot)/dp
    jac(1,1:3)=-du_sup_x(1:3)
    jac(2,1:3)=-du_sup_y(1:3)
    jac(3,1:3)=-du_sup_z(1:3)
    !d(rdot)/dr
    jac(4:6,4)=du_sup_x(1:3)
    jac(4:6,5)=du_sup_y(1:3)
    jac(4:6,6)=du_sup_z(1:3)
    !d(pdot)/dr
    jac(1,4)=-dot_product(g(i)%p,du_sup_x2) ; jac(1,5)=-dot_product(g(i)%p,du_sup_xy) ; jac(1,6)=-dot_product(g(i)%p,du_sup_zx)
    jac(2,4)=-dot_product(g(i)%p,du_sup_xy) ; jac(2,5)=-dot_product(g(i)%p,du_sup_y2) ; jac(2,6)=-dot_product(g(i)%p,du_sup_yz)
    jac(3,4)=-dot_product(g(i)%p,du_sup_zx) ; jac(3,5)=-dot_product(g(i)%p,du_sup_yz) ; jac(3,6)=-dot_product(g(i)%p,du_sup_z2)
    !finaly section that needs doing occupies jac(4:6,1:3)
    !d(rdot)/dp
  end subroutine
end module
