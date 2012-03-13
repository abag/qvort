!>timestepping and velocity routines are contained in this module
module timestep
  use cdata
  use general
  use normal_fluid
  use hyperviscous
  use forcing
  use tree
  use mirror
  contains
  !*************************************************
  ! find the position and velocity of the
  ! centre of vorticity of all vortex points
  subroutine get_cov_data
    implicit none
    integer :: i
    do i=1,3
      cov%x(i)=sum(f(:)%x(i))/count(mask=f(:)%infront>0)
    end do
    if(itime.gt.1) then
      cov%u=(cov%x-cov%x_old)/dt
    else
      cov%u=0.0
    end if
    cov%x_old=cov%x
  end subroutine
  ! now find macro ring radii (R and a) - 2 routines
  subroutine get_macro_ring_radii_1 ! first routine to find macro ring radii (R and a)
    implicit none
    integer :: i
    do i=1,3   
      mrr1%x_max(i)=maxval(f(:)%x(i),mask=f(:)%infront>0)
      mrr1%x_min(i)=minval(f(:)%x(i),mask=f(:)%infront>0)
    end do    
    mrr1%x_max(4)=(maxval(f(:)%x(2),mask=f(:)%infront>0)+maxval(f(:)%x(2),mask=f(:)%infront>0))/2.;
    mrr1%x_min(4)=(minval(f(:)%x(2),mask=f(:)%infront>0)+minval(f(:)%x(2),mask=f(:)%infront>0))/2.;
    do i=1,4
      mrr1%x_spread(i)=mrr1%x_max(i)-mrr1%x_min(i)
    end do
    mrr1%r=(mrr1%x_spread(4)-mrr1%x_spread(1))/2.
    mrr1%a=mrr1%x_spread(1)/2.
    mrr1%ra=mrr1%r/mrr1%a
    if(itime.gt.1) then
      mrr1%r_u=(mrr1%r-mrr1%r_old)/dt
      mrr1%a_u=(mrr1%a-mrr1%a_old)/dt
    else
      mrr1%r_u=0.0
      mrr1%a_u=0.0
    end if
    mrr1%r_old=mrr1%r
    mrr1%a_old=mrr1%a   
  end subroutine
  subroutine get_macro_ring_radii_2 ! second routine to find macro ring radii (R and a)
    implicit none
    integer :: i
    allocate(mrr2(pcount))
    do i=1,3
      mrr2(:)%r(i)=sqrt((f(:)%x(i)-cov%x(i))**2)
    end do
    mrr2(:)%r(4)=sqrt((f(:)%x(2)-cov%x(2))**2 + (f(:)%x(3)-cov%x(3))**2) 
    mrr2(:)%r(5)=sqrt((f(:)%x(1)-cov%x(1))**2 + (f(:)%x(2)-cov%x(2))**2 + (f(:)%x(3)-cov%x(3))**2)
    do i=1,5
      avg_r(i)=sum(mrr2(:)%r(i))/count(mask=f(:)%infront>0)
      mrr2(:)%a(i)=abs(mrr2(:)%r(i)-avg_r(i))
      avg_a(i)=sum(mrr2(:)%a(i))/count(mask=f(:)%infront>0)
      avg_ra(i)=avg_r(i)/avg_a(i)
    end do
    if(itime.gt.1) then
      avg_r_u=(avg_r-avg_r_old)/dt
      avg_a_u=(avg_a-avg_a_old)/dt
    else
      avg_r_u=0.0
      avg_a_u=0.0
    end if
    avg_r_old=avg_r
    avg_a_old=avg_a
    deallocate(mrr2)
  end subroutine
  !*************************************************
  !>implement adams bashforth time-stepping scheme to move particles
  !>\f[
  !! \mathbf{s}_{i}^{n+1}=\mathbf{s}_{i}^{n}+\frac{\Delta t}{12}(23\mathbf{u}_{i}^{n}
  !! -16\mathbf{u}_{i}^{n-1}+5\mathbf{u}_{i}^{n-2})+\mathcal{O}(\Delta t^4)
  !>\f]
  !>adaptive timestep can be set in run.in in which case a comparison between 2nd and 3rd
  !>order scheme is used to estimate error and hence adjust the timestep
  subroutine pmotion()
    implicit none
    real :: u(3) !dummy variable used to store velocities
    real :: adap_x(3) !dummy variable to check timestep
    real :: dummy_max_error, max_error !maximum error between 2nd and 3rd order method
    real :: rot_r, rot_theta !for differential rotation
    integer :: i
    !intialise forcing every timestep incase there is a random element
    call randomise_forcing
    !begin by testing if we have special velocity field Rotate
    select case(velocity)
      case('Rotate')
        call differential_rotation
        return !exit the routine once applied
    end select
    !any quantaties that need zeroing before velocity calculation here    
    max_error=0. 
    hyp_power_dissipate=0.
    !now loop over all points and get the velocity field
    !$omp parallel do private(i,u)
    do i=1, pcount
      if (f(i)%infront==0) cycle !check for 'empty' particles
      call calc_velocity(u,i)
      f(i)%u(:)=u(:) !store the velocity for time-step
    end do
    !$omp end parallel do
    !$omp parallel do private(i)
    do i=1, pcount
      if (f(i)%infront==0) cycle !check for 'empty' particles
      if (maxval(abs(f(i)%u1))==0) then
        f(i)%x(:)=f(i)%x(:)+dt*f(i)%u(:) !euler
      else if (maxval(abs(f(i)%u2))==0) then
        !first order adams-bashforth
        f(i)%x(:)=f(i)%x(:)+three_twos*dt*f(i)%u(:)-one_half*dt*f(i)%u1(:)
      else if (maxval(abs(f(i)%u3))==0) then
        !second order adams-bashforth
        f(i)%x(:)=f(i)%x(:)+twenty_three_twelve*dt*f(i)%u(:)-four_thirds*dt*f(i)%u1(:)+five_twelths*dt*f(i)%u2(:)
       else
        !3rd order adams-bashforth
        f(i)%x(:)=f(i)%x(:)+dt*(55.*f(i)%u(:)-59.*f(i)%u1(:)+37.*f(i)%u2(:)-9.*f(i)%u3(:))/24.
        if (dt_adapt) then !adaptive timestep (read in from cdata)
          adap_x(:)=f(i)%x(:)+three_twos*dt*f(i)%u(:)-one_half*dt*f(i)%u1(:)
          dummy_max_error=dist_gen(f(i)%x,adap_x)
          if (dummy_max_error>max_error) max_error=dummy_max_error
        end if
      end if
      f(i)%u3(:)=f(i)%u2(:) !store our old velocities 
      f(i)%u2(:)=f(i)%u1(:)  
      f(i)%u1(:)=f(i)%u(:)
    end do
    !$omp end parallel do
    !adjust timestep
    if (dt_adapt) then
      !at present done every 5 timesteps this needs to be experimented with
      if (mod(itime,5)==0) then
        dt=dt*(1E-4/(2*max_error))**(1./3.)
        open(unit=34,file='data/adaptive_error.log',position='append')
          write(34,*) t, max_error, (1E-6/(2*max_error))**(1./3.), dt
        close(34)
        !min and max timestep need to be set in run.in
        if (dt<1E-8) dt=1E-8 !min timestep
        if (dt>1E-5) dt=1E-5 !max timestep
      end if
    end if 
  end subroutine
  !*************************************************
  !>get the velocity of each particle subject to the superfluid velocity
  !>plus any normal fluid/forcing
  subroutine calc_velocity(u,i)
    implicit none
    integer, intent(IN) :: i
    real :: u(3), u_norm(3), u_force(3), u_bs(3), u_mir(3), u_B(3)!velocities
    real :: cov(3), cov_vel(3) !centre of vorticity
    real :: curv, beta !LIA
    real :: hyperviscous_alpha !for hyperviscosity
    real :: f_dot(3), f_ddot(3) !first and second derivs
    integer :: peri, perj, perk !used to loop in periodic cases
   ! integer :: miri, mirj, mirk !used to loop in mirror cases
    !what scheme are we using? (LIA/BS)
    select case(velocity)
      case('Off')
        !no superfluid velocity
        u=0.
      case('LIA')
        !use the local induction approximation
        call get_deriv_1(i,f_dot) !general.mod
        call get_deriv_2(i,f_ddot) !general.mod
        !****************do we vary beta?*******************
        !calculate the curvature
        if (fixed_LIA_beta) then
          beta=1.3E-3
        else
          curv=sqrt(dot_product(f_ddot,f_ddot))
          if (curv<epsilon(0.)) then
            curv=epsilon(0.) !we must check for zero curvature
          else
            curv=1./curv
          end if
          !caluculate beta based on the curvature
          beta=(quant_circ/(4.*pi))*log(4.6*curv/corea)
        end if
        !***************************************************
        u=beta*cross_product(f_dot,f_ddot) !general.mod
        f(i)%u_sup=u !store just the superfluid velcotity
      case('BS')
        !full (nonlocal) biot savart
        !first get the local part (similar to LIA)
        call get_deriv_1(i,f_dot) !general.mod
        call get_deriv_2(i,f_ddot) !general.mod
        beta=(quant_circ/(4.*pi))*log((1.1/corea)*sqrt(dist_gen(f(i)%x,f(i)%ghosti)*dist_gen(f(i)%x,f(i)%ghostb)))
        u=beta*cross_product(f_dot,f_ddot) !general.mod
        !now we do the non-local part
        u_bs=0. !always 0 before calling the routine
        call biot_savart(i,u_bs)
        u=u+u_bs
        !we now need to insert periodic and mirror boundary conditions here
        if (periodic_bc) then
          !we must shift the mesh in all 3 directions
          u_bs=0. !zero u_bs
          do peri=-1,1 ; do perj=-1,1 ; do perk=-1,1
            if (peri==0.and.perj==0.and.perk==0) cycle
            call biot_savart_shift(i,u_bs,(/peri*box_size,perj*box_size,perk*box_size/))
          end do ; end do ;end do
          u=u+u_bs
        else if (periodic_bc_notx) then
          !we must shift the mesh in 2 directions (yz)
          u_bs=0. !zero u_bs
          do perj=-1,1 ; do perk=-1,1
            if (perj==0.and.perk==0) cycle
            call biot_savart_shift(i,u_bs,(/0.,perj*box_size,perk*box_size/))
          end do ;end do
          u=u+u_bs
        else if (periodic_bc_notxy) then
          !we must shift the mesh in z direction
          u_bs=0. !zero u_bs
          do perk=-1,1
            if (perk==0) cycle
            call biot_savart_shift(i,u_bs,(/0.,0.,perk*box_size/))
          end do
          u=u+u_bs
        end if
        if (mirror_bc) then
          !get the contribution of the image vortices
          u_mir=0.
          call biot_savart_mirror(i,u_mir) !mirror.mod
          u=u+u_mir
        end if
        f(i)%u_sup=u !store just the superfluid velcotity
      case('Tree')
        !tree approximation to biot savart
        !first get the local part (similar to LIA)
        call get_deriv_1(i,f_dot) !general.mod
        call get_deriv_2(i,f_ddot) !general.mod
        beta=(quant_circ/(4.*pi))*log((1.1/corea)*sqrt(dist_gen(f(i)%x,f(i)%ghosti)*dist_gen(f(i)%x,f(i)%ghostb)))
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
        else if (periodic_bc_notx) then
          !we must shift the mesh in 2 directions, 9 permutations needed!
          u_bs=0. !zero u_bs
          do perj=-1,1 ; do perk=-1,1
            if (perj==0.and.perk==0) cycle
            call tree_walk(i,vtree,(/0.,perj*box_size,perk*box_size/),u_bs) !tree.mod
          end do ;end do
          u=u+u_bs
        else if (periodic_bc_notxy) then
          !we must shift the mesh in z direction, 3 permutations needed, not so bad
          u_bs=0. !zero u_bs
          do perk=-1,1
            if (perk==0) cycle
            call tree_walk(i,vtree,(/0.,0.,perk*box_size/),u_bs) !tree.mod
          end do
          u=u+u_bs
        end if
        f(i)%u_sup=u !store just the superfluid velcotity
    end select
    !hyperviscosity
    if (hyperviscosity) then
      call get_hyp_alpha(sqrt(dot_product(f_ddot,f_ddot)),hyperviscous_alpha) !hyperviscous.mod
      if (mod(itime,shots)==0) then!only calculate P_{mf} every shots timesteps
        !P_{mf}=-\alpha*\rho_s\Gamma \int |s' \times v_s| d\xi
        hyp_power_dissipate=hyp_power_dissipate-hyperviscous_alpha*quant_circ*&
                          vector_norm(cross_product(f_dot,u))*dist_gen(f(i)%x,f(i)%ghosti)
      end if
      if (hyperviscous_alpha>epsilon(0.)) then
        u=u+hyperviscous_alpha*cross_product(f_dot,(-u))
      end if
    end if
    !now account for mutual friction - test if alpha's are 0
    select case(velocity)
      case('Off')
        call get_normal_velocity(f(i)%x,u_norm) !normal_fluid.mod
        u=u_norm
      case default
        !check that either of the mutual friction coefficients are >0
        if ((abs(alpha(1))>epsilon(0.)).or.(abs(alpha(2))>epsilon(0.))) then
          if (t<normal_fluid_cutoff) then !cutoff time set in run.in
            !omp critical
            call get_normal_velocity(f(i)%x,u_norm) !normal_fluid.mod
            !omp end critical
            ! \todo this could be improved calculating same thing twice 
            f(i)%u_mf=alpha(1)*cross_product(f_dot,(u_norm-u))- &
                      alpha(2)*cross_product(f_dot,cross_product(f_dot,(u_norm-u)))
            u=u+f(i)%u_mf !this way we store the mutual friction velocity
          end if
        end if
    end select
    !forcing?
    if (t<force_cutoff) then
      call get_forcing(i,u_force)
      u=u+u_force
    end if
    if (mirror_bc) then
      !check the flux through the boundaries is 0
      call mirror_flux_check(i,u) !mirror.mod
    end if
    if (sticky_z_boundary) then
      !particles at top/bottom of box are fixed
      if ((abs(f(i)%x(3))-box_size/2.)>-1.5*delta) then
        !particle is sufficiently close to top boundary stick
        u=0.
      end if
    end if
  end subroutine
  !**************************************************************************
  !>if the mesh size is set to be larger than 0 in run.in then we calculate both
  !>the normal and superfluid velocities at the mesh points to print to file
  subroutine mesh_velocity
    !get the velocity at each point on the mesh for spectra etc.
    implicit none
    integer :: i,j,k
    integer :: peri, perj, perk !used to loop in periodic cases
    real :: normurms
   !$omp parallel do private(i,j,k,peri,perj,perk)
    do k=1, mesh_size
      if (mesh_size>=64) then
        write(*,'(a,i4.1,a,i4.1)') 'drawing mesh section ', k, ' of ', mesh_size
      end if
      do j=1, mesh_size
        do i=1, mesh_size
          !superfluid velocity
          select case(velocity)
            case('Off, LIA')
              if (itime==1) then
                write(*,*) 'WARNING, LIA is being used, no superfluid velocity will be calculated on mesh'
              end if
            case('BS')
              mesh(k,j,i)%u_sup=0. !always 0 before making initial BS call
              call biot_savart_general(mesh(k,j,i)%x,mesh(k,j,i)%u_sup)
              if (periodic_bc) then
                !we must shift the mesh in all 3 directions, all 26 permutations needed!
                do peri=-1,1 ; do perj=-1,1 ; do perk=-1,1
                  if (peri==0.and.perj==0.and.perk==0) cycle
                  call biot_savart_general_shift(mesh(k,j,i)%x,mesh(k,j,i)%u_sup, &
                  (/peri*box_size,perj*box_size,perk*box_size/)) !timestep.mod
                end do ; end do ; end do 
              else if (periodic_bc_notx) then
                !as above but we do not need the x permutations
                do perj=-1,1 ; do perk=-1,1
                  if (perj==0.and.perk==0) cycle
                  call biot_savart_general_shift(mesh(k,j,i)%x,mesh(k,j,i)%u_sup, &
                  (/0.,perj*box_size,perk*box_size/)) !timestep.mod
                end do ; end do 
              else if (periodic_bc_notxy) then
                !as above but we do not need the x permutations
                do perk=-1,1
                  if (perk==0) cycle
                  call biot_savart_general_shift(mesh(k,j,i)%x,mesh(k,j,i)%u_sup, &
                  (/0.,0.,perk*box_size/)) !timestep.mod
                end do 
              end if
            case('Tree')
              mesh(k,j,i)%u_sup=0. !must be zeroed for all algorithms
              call tree_walk_general(mesh(k,j,i)%x,vtree,(/0.,0.,0./),mesh(k,j,i)%u_sup)
              if (periodic_bc) then
                !we must shift the mesh in all 3 directions, all 26 permutations needed!
                do peri=-1,1 ; do perj=-1,1 ; do perk=-1,1
                  if (peri==0.and.perj==0.and.perk==0) cycle
                  call tree_walk_general(mesh(k,j,i)%x,vtree, &
                       (/peri*box_size,perj*box_size,perk*box_size/),mesh(k,j,i)%u_sup) !tree.mod
                end do ; end do ; end do
              else if (periodic_bc_notx) then
                !as above but we do not need the x permutations
                do perj=-1,1 ; do perk=-1,1
                  if (perj==0.and.perk==0) cycle
                  call tree_walk_general(mesh(k,j,i)%x,vtree, &
                       (/0.,perj*box_size,perk*box_size/),mesh(k,j,i)%u_sup) !tree.mod
                end do ; end do
              else if (periodic_bc_notxy) then
                !as above but we do not need the x/y permutations
                do perk=-1,1
                  if (perk==0) cycle
                  call tree_walk_general(mesh(k,j,i)%x,vtree, &
                       (/0.,0.,perk*box_size/),mesh(k,j,i)%u_sup) !tree.mod
                end do
              end if
          end select
          !normal fluid
          !omp critical
          call get_normal_velocity(mesh(k,j,i)%x,mesh(k,j,i)%u_norm)
          !omp end critical
        end do
      end do
    end do
    !$omp end parallel do
  end subroutine
  !**************************************************************************
  !>the desingularised biot savart integral 
  !>\f[
  !>\frac{d\mathbf{s}_i}{dt}=\frac{\Gamma}{4\pi} \ln \left(\frac{\sqrt{\ell_i
  !>\ell_{i+1}}}{a}\right)\mathbf{s}_i' \times \mathbf{s}_i'' 
  !>+\frac{\Gamma}{4 \pi} \oint_{\cal L'} \frac{(\mathbf{s}_i-\mathbf{r}) }
  !>{\vert \mathbf{s}_i - \mathbf{r} \vert^3}
  !>\times {\bf d}\mathbf{r}
  !>\f] note the LIA part is calculated in calc_velocity
  subroutine biot_savart(i,u)
    implicit none
    integer, intent(IN) :: i !the particle we want the velocity at
    real :: u(3) !the velocity at i (non-local)#
    real :: u_bs(3) !helper array
    real :: a_bs, b_bs, c_bs !helper variables
    integer :: j !needed to loop over all particles
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
  !>as above but shifts the particles (by a vector shift) for periodicity
  subroutine biot_savart_shift(i,u,shift)
    implicit none
    integer, intent(IN) :: i !the particle we want the velocity at
    real :: u(3) !the velocity at i (non-local)#
    real :: u_bs(3) !helper array
    real :: a_bs, b_bs, c_bs !helper variables
    real :: shift(3) !moves all the particles for periodic_bc
    integer :: j !needed to loop over all particles
    do j=1, pcount
      !check that the particle is not empty or the particle behind
      !if ((f(j)%infront==0).or.(i==j).or.(f(i)%behind==j)) cycle
      if ((f(j)%infront==0).or.(f(i)%behind==j)) cycle
      a_bs=dist_gen_sq(f(i)%x,f(j)%x+shift) !distance squared between i and j+shift
      b_bs=2.*dot_product((f(j)%x+shift-f(i)%x),(f(j)%ghosti-f(j)%x))
      c_bs=dist_gen_sq(f(j)%ghosti,f(j)%x) !distance sqd between j, j+1
      !add non local contribution to velocity vector
      if (4*a_bs*c_bs-b_bs**2==0) cycle !avoid 1/0
      u_bs=cross_product((f(j)%x+shift-f(i)%x),(f(j)%ghosti-f(j)%x))
      u_bs=u_bs*quant_circ/((2*pi)*(4*a_bs*c_bs-b_bs**2))
      u_bs=u_bs*((2*c_bs+b_bs)/sqrt(a_bs+b_bs+c_bs)-(b_bs/sqrt(a_bs)))
      u=u+u_bs !add on the non-local contribution of j
    end do
  end subroutine
  !**************************************************************************
  !>calculate the velocity field at a point x induced by the vortices now using
  !>more general biot savart integral
  !>\f[\frac{\Gamma}{4 \pi} \oint_{\cal L} \frac{(\mathbf{s}_i-\mathbf{r}) }
  !>{\vert \mathbf{s}_i - \mathbf{r} + \epsilon \vert^3}
  !>\times {\bf d}\mathbf{r}.
  !>\f]
  !>This is only used to calculate the velocity at points on a cartesian mesh.
  subroutine biot_savart_general(x,u)
    implicit none
    real, intent(IN) :: x(3)
    real, intent(OUT) :: u(3)
    real :: u_bs(3) !helper vector
    real :: a_bs, b_bs, c_bs !helper variables
    integer :: j !needed to loop over all particles 
    integer :: i
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
  !**************************************************************************
  !>as above but allows a shift for periodicity
  !>\todo makes above redundent so look into removal
  subroutine biot_savart_general_shift(x,u,shift)
    implicit none
    real, intent(IN) :: x(3)
    real, intent(OUT) :: u(3)
    real :: u_bs(3) !helper vector
    real :: shift(3) !moves all the particles for periodic_bc
    real :: a_bs, b_bs, c_bs !helper variables
    integer :: j !needed to loop over all particles 
    integer :: i
    do i=1, pcount
      !check that the particle is not empty
      if (f(i)%infront==0) cycle
      a_bs=dist_gen_sq(x,f(i)%x+shift) !distance squared between x and particle
      b_bs=2.*dot_product((f(i)%x+shift-x),(f(i)%ghosti-f(i)%x))
      c_bs=dist_gen_sq(f(i)%ghosti,f(i)%x) !distance sqd between i, i+1
      !add non local contribution to velocity vector
      if ((4*a_bs*c_bs-b_bs**2)==0) cycle !avoid 1/0!
      u_bs=cross_product((f(i)%x+shift-x),(f(i)%ghosti-f(i)%x))
      u_bs=u_bs*quant_circ/((2*pi)*(4*a_bs*c_bs-b_bs**2))
      u_bs=u_bs*((2*c_bs+b_bs)/sqrt(a_bs+b_bs+c_bs)-(b_bs/sqrt(a_bs)))
      u=u+u_bs !add on the non-local contribution of i
    end do
  end subroutine
  !**************************************************************************
  !>apply differential rotation to all the particles
  subroutine differential_rotation()
    implicit none
    real :: rot_r, rot_theta !for cylindrical coords
    integer :: i !for looping
    do i=1, pcount
      if (f(i)%infront==0) cycle !check for 'empty' particles
      rot_r=sqrt(f(i)%x(1)**2+f(i)%x(3)**2)
      rot_theta=atan2(f(i)%x(3),f(i)%x(1))
      rot_theta=rot_theta+dt*(0.728*exp(-1.6666*(6.*rot_r/box_size)**2))
      if (rot_theta>pi) rot_theta=rot_theta-2*pi
      f(i)%x(1)=rot_r*cos(rot_theta) 
      f(i)%x(3)=rot_r*sin(rot_theta) 
    end do
  end subroutine
end module
