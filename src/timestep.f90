module timestep
  !ANYTHING RELATED TO VELOCITY SHOULD ENTER HERE
  use cdata
  use general
  use normal_fluid
  use forcing
  use tree
  contains
  !*************************************************
  subroutine pmotion()
    !implement adams bashforth time-stepping scheme to move particles
    implicit none
    real :: u(3) !dummy variable used to store velocities
    integer :: i
    do i=1, pcount
      if (f(i)%infront==0) cycle !check for 'empty' particles
      call calc_velocity(u,i)
      f(i)%u(:)=u(:) !store the velocity for time-step
    end do
    do i=1, pcount
      if (f(i)%infront==0) cycle !check for 'empty' particles
      if (maxval(abs(f(i)%u1))==0) then
        f(i)%x(:)=f(i)%x(:)+dt*f(i)%u(:) !euler
      else if (maxval(abs(f(i)%u2))==0) then
        !first order adams-bashforth
        f(i)%x(:)=f(i)%x(:)+three_twos*dt*f(i)%u(:)-one_half*dt*f(i)%u1(:)
      else
        !second order adams-bashforth
        f(i)%x(:)=f(i)%x(:)+twenty_three_twelve*dt*f(i)%u(:)-four_thirds*dt*f(i)%u1(:)+five_twelths*dt*f(i)%u2(:)
      end if
      f(i)%u2(:)=f(i)%u1(:) !store our old velocities 
      f(i)%u1(:)=f(i)%u(:)
    end do
    t=t+dt !time needs incrementing 
  end subroutine
  !*************************************************
  subroutine calc_velocity(u,i)
    !get the velocity of each particle subject to the superfluid velocity
    !plus any normal fluid/forcing
    implicit none
    integer, intent(IN) :: i
    real :: u(3), u_norm(3), u_force(3), u_bs(3) !velocities
    real :: curv, beta !LUA
    real :: f_dot(3), f_ddot(3) !first and second derivs
    integer :: peri, perj, perk !used to loop in periodic cases
    !what scheme are we using? (LIA/BS)
    select case(velocity)
      case('LIA')
        !use the local induction approximation
        call get_deriv_1(i,f_dot) !general.mod
        call get_deriv_2(i,f_ddot) !general.mod
        !IF WE WANT BETA TO VARY DO THIS********************
        !calculate the curvature
        !curv=sqrt(dot_product(f_ddot,f_ddot))
        !if (curv<epsilon(0.)) then
        !  curv=epsilon(0.) !we must check for zero curvature
        !end if
        !caluculate beta based on the curvature
        !beta=(quant_circ/(4.*pi))*log(1.213E8/curv)
        !***************************************************
        beta=1.3E-3
        u=beta*cross_product(f_dot,f_ddot) !general.mod
      case('BS')
        !full (nonlocal) biot savart
        !first get the local part (similar to LIA)
        call get_deriv_1(i,f_dot) !general.mod
        call get_deriv_2(i,f_ddot) !general.mod
        beta=(quant_circ/(4.*pi))*log(1.213E8*sqrt(distf(i,f(i)%infront)*distf(i,f(i)%behind)))
        u=beta*cross_product(f_dot,f_ddot) !general.mod
        !now we do the non-local part
        call biot_savart(i,u_bs)
        u=u+u_bs
      case('Tree')
        !tree approximation to biot savart
        !first get the local part (similar to LIA)
        call get_deriv_1(i,f_dot) !general.mod
        call get_deriv_2(i,f_ddot) !general.mod
        beta=(quant_circ/(4.*pi))*log(1.213E8*sqrt(distf(i,f(i)%infront)*distf(i,f(i)%behind)))
        u=beta*cross_product(f_dot,f_ddot) !general.mod
        !now walk the tree to get the non-local contribution
        u_bs=0. !zero u_bs
        call tree_walk(i,vtree,(/0.,0.,0./),u_bs) !tree.mod
        u=u+u_bs
        if (periodic_bc) then
          !we must shift the mesh in all 3 directions, all 26 permutations needed!
          u_bs=0. !zero u_bs
          do peri=-1,1 ; do perj=-1,1 ; do perk=-1,1
            if (peri==0.and.perj==0.and.perk==0) cycle
            call tree_walk(i,vtree,(/peri*box_size,perj*box_size,perk*box_size/),u_bs) !tree.mod
          end do ; end do ;end do
          u=u+u_bs
        end if
    end select
    !now account for mutual friction - test if alpha's are 0
    if ((sum(alpha)>epsilon(0.)).and.(t<normal_fluid_cutoff)) then
      call get_normal_velocity(f(i)%x,u_norm) !normal_fluid.mod
      u=u+alpha(1)*cross_product(f_dot,(u_norm-u))- &
          alpha(2)*cross_product(f_dot,cross_product(f_dot,(u_norm-u)))
    end if
    !forcing?
    call get_forcing(i,u_force)
    u=u+u_force
  end subroutine
  !**************************************************************************
  subroutine mesh_velocity
    !get the velocity at each point on the mesh for spectra etc.
    implicit none
    integer :: i,j,k
    integer :: peri, perj, perk !used to loop in periodic cases
    real :: normurms
    do k=1, mesh_size
      do j=1, mesh_size
        do i=1, mesh_size
          !superfluid velocity
          select case(velocity)
            case('LIA')
              if (itime==1) then
                write(*,*) 'WARNING, LIA is being used, no superfluid velocity will be calculated on mesh'
              end if
            case('BS')
              call biot_savart_general(mesh(k,j,i)%x,mesh(k,j,i)%u_sup)
            case('Tree')
              mesh(k,j,i)%u_sup=0. !must be zeroed for tree algorithms
              call tree_walk_general(mesh(k,j,i)%x,vtree,(/0.,0.,0./),mesh(k,j,i)%u_sup)
              if (periodic_bc) then
                !we must shift the mesh in all 3 directions, all 26 permutations needed!
                do peri=-1,1 ; do perj=-1,1 ; do perk=-1,1
                  if (peri==0.and.perj==0.and.perk==0) cycle
                  call tree_walk_general(mesh(k,j,i)%x,vtree, &
                       (/peri*box_size,perj*box_size,perk*box_size/),mesh(k,j,i)%u_sup) !tree.mod
                end do ; end do ;end do
              end if
          end select
          !normal fluid
          call get_normal_velocity(mesh(k,j,i)%x,mesh(k,j,i)%u_norm)
        end do
      end do
    end do
  end subroutine
  !**************************************************************************
  subroutine biot_savart(i,u)
    implicit none
    integer, intent(IN) :: i !the particle we want the velocity at
    real :: u(3) !the velocity at i (non-local)#
    real :: u_bs(3) !helper array
    real :: a_bs, b_bs, c_bs !helper variables 
    integer :: j !needed to loop over all particles
    u=0.
    do j=1, pcount
      !check that the particle is not empty/i/f(i)%behind
      if ((f(j)%infront==0).or.(i==j).or.(f(i)%behind==j)) cycle
      a_bs=distfsq(j,i) !distance squared between i and j
      b_bs=2.*dot_product((f(j)%x-f(i)%x),(f(j)%ghosti-f(j)%x))
      c_bs=dist_gen_sq(f(j)%ghosti,f(j)%x) !distance sqd between j, j+1
      !add non local contribution to velocity vector
      if (4*a_bs*c_bs-b_bs**2==0) cycle !avoid 1/0
      u_bs=cross_product((f(j)%x-f(i)%x),(f(j)%ghosti-f(j)%x))
      u_bs=u_bs*quant_circ/((2*pi)*(4*a_bs*c_bs-b_bs**2))
      u_bs=u_bs*((2*c_bs+b_bs)/sqrt(a_bs+b_bs+c_bs)-(b_bs/sqrt(a_bs)))
      u=u+u_bs !add on the non-local contribution of j
    end do
  end subroutine
  !**************************************************************************
  subroutine biot_savart_general(x,u)
    !calculate the velocity field at a point x induced by the vortices
    implicit none
    real, intent(IN) :: x(3)
    real, intent(OUT) :: u(3)
    real :: u_bs(3) !helper vector
    real :: a_bs, b_bs, c_bs !helper variables
    integer :: j !needed to loop over all particles 
    !what scheme are we using? (LIA/BS)
    integer :: i
    u=0.
    do i=1, pcount
      !check that the particle is not empty
      if (f(i)%infront==0) cycle
      a_bs=dist_gen_sq(x,f(i)%x) !distance squared between x and particle
      b_bs=2.*dot_product((f(i)%x-x),(f(i)%ghosti-f(i)%x))
      c_bs=dist_gen_sq(f(i)%ghosti,f(i)%x) !distance sqd between i, i+1
      !add non local contribution to velocity vector
      if ((4*a_bs*c_bs-b_bs**2)==0) cycle !avoid 1/0!
      u_bs=cross_product((f(i)%x-x),(f(i)%ghosti-f(i)%x))
      u_bs=u_bs*quant_circ/((2*pi)*(4*a_bs*c_bs-b_bs**2))
      u_bs=u_bs*((2*c_bs+b_bs)/sqrt(a_bs+b_bs+c_bs)-(b_bs/sqrt(a_bs)))
      u=u+u_bs !add on the non-local contribution of i
    end do
  end subroutine
end module
