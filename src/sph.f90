!>SPH smooth paricle hydrodynamics - see \ref SPH
module sph
  use cdata
  use general
  use output
  use sph_tree
  use sph_kernel
  implicit none
  integer,parameter,private :: neighnumb=25 !number of neighbours to use
  real,parameter :: h_eta=0.5*(0.75*neighnumb/pi)**(0.333333333333) !works with above
  !**********SPH MESH STRUCTURE***********
  !>mesh structure - at present only density stored
  !>runs from -box_size /2 to box_size/2 to with mesh_size^3
  !>points, note by setting mesh_size>0 mesh is switched on
  !>@param x position of mesh point
  !>@param rho density
  type SPH_grid
    real :: x(3)
    real :: rho 
  end type
  !>3D allocatable mesh
  type(SPH_grid), allocatable :: SPH_mesh(:,:,:)
  contains
  !************************************************************
  !>evolve the SPH particles
  subroutine SPH_evolution
    implicit none
    integer :: i !for looping
    !----------------create tree-------------------
    if (SPH_theta>0) call create_SPH_tree !sph_tree.mod
    !----------------------------------------------
    !do operations involving neighbouring particles
    call sph_neighbour_loop 
    !now get the acc. due to non local (gravity) forces
    call sph_get_non_loc
    do i=1,SPH_count
      !---------timestep the velocity--------------
      if (maxval(abs(s(i)%a1))==0) then
        s(i)%u=s(i)%u+dt*s(i)%a !euler
      else if (maxval(abs(s(i)%a2))==0) then
        !first order adams-bashforth
        s(i)%u=s(i)%u+three_twos*dt*s(i)%a-one_half*dt*s(i)%a1
      else
        !second order adams-bashforth
        s(i)%u=s(i)%u+twenty_three_twelve*dt*s(i)%a-four_thirds*dt*s(i)%a1+five_twelths*dt*s(i)%a2
      end if
      !---------timestep the position--------------
      if (maxval(abs(s(i)%u1))==0) then
        s(i)%x=s(i)%x+dt*s(i)%u !euler
      else if (maxval(abs(s(i)%u2))==0) then
        !first order adams-bashforth
        s(i)%x=s(i)%x+three_twos*dt*s(i)%u-one_half*dt*s(i)%u1
      else
        !second order adams-bashforth
        s(i)%x=s(i)%x+twenty_three_twelve*dt*s(i)%u-four_thirds*dt*s(i)%u1+five_twelths*dt*s(i)%u2
      end if
      s(i)%a2(:)=s(i)%a1(:) ; s(i)%a1(:)=s(i)%a(:) !store old acceleration 
      s(i)%u2(:)=s(i)%u1(:) ; s(i)%u1(:)=s(i)%u(:)!store old velocities 
    end do
    !timestep
    call time_check_SPH
    !check to see if we print
    if (mod(itime,shots)==0) then
      call diagnostics_SPH !sph.mod
      call print_SPH(itime/shots) !output.mod
    end if
    !do we have a mesh?
    if (mod(itime,mesh_shots)==0) then
      if (SPH_mesh_size>0) then
        call mesh_print_SPH(itime/mesh_shots)  !sph.mod
      end if
    end if
    if (SPH_theta>0) then
      call SPH_empty_tree(stree) !empty the tree to avoid a memory leak
      deallocate(stree%parray) ; deallocate(stree) ; nullify(stree)
    end if
  end subroutine
  !************************************************************
  !>do all operations which only involve nearest neighbours
  subroutine sph_neighbour_loop
    implicit none
    integer :: i, j
    real, dimension(3) :: rij_hat
    real :: dist, rho_bar, nu_sig
    !get the nearest neighbours for all particles
    if (SPH_theta>0) then
      call SPH_neighbour_list !sph.mod
    else
      call SPH_neighbour_list_notree !sph_tree.mod
    end if
    do i=1, SPH_count
      s(i)%rho=0. ; s(i)%drhodh=0. !0 these quantities
      !%%%%%%%%%%%%%%%%% DENSITY LOOP %%%%%%%%%%%%%%%%%
      s_NN(i)%current => s_NN(i)%list
      do while (associated( s_NN(i)%current ))
        j=s_NN(i)%current%i
        dist=dist_gen(s(i)%x,s(j)%x)
        !---------------density----------------
        s(i)%rho=s(i)%rho+s(j)%m*sph_W(dist,s(i)%h)
        !---------------density deriv wrt h----------------
        s(i)%drhodh=s(i)%drhodh+s(j)%m*sph_dWdh(dist,s(i)%h)
        !---------move down neighbour list-----
        s_NN(i)%previous => s_NN(i)%current
        s_NN(i)%current => s_NN(i)%current%next
      end do
      if (s(i)%rho<epsilon(0.)) then
        s(i)%u=0. ; s(i)%a=0. 
        cycle !exit the loop when particles have left the box
      end if
      s(i)%h=h_eta*(s(i)%m/s(i)%rho)**(0.333333333333)   
      !---------set pressure-------
      s(i)%P=s(i)%rho**SPH_gamma
      !set the correction term according to springel 2010 (review)
      s(i)%f=1./(1+(h_eta*s(i)%h/(3.*s(i)%rho))*s(i)%drhodh)
      s(i)%a=0. ; s(i)%divu=0.
      !%%%%%%%%%%%%%%%%% VELOCITY LOOP %%%%%%%%%%%%%%%%%
      s_NN(i)%current => s_NN(i)%list
      do while (associated( s_NN(i)%current ))
        j=s_NN(i)%current%i
        if (i/=j) then
          dist=dist_gen(s(i)%x,s(j)%x)
          !---------------hydrodynaics----------------
          s(i)%a=s(i)%a-s(i)%f*s(j)%m*(s(i)%P/(s(i)%rho**2))*sph_grad_W(s(i)%x,s(j)%x,s(i)%h)
          s(i)%a=s(i)%a-s(j)%f*s(j)%m*(s(j)%P/(s(j)%rho**2))*sph_grad_W(s(i)%x,s(j)%x,s(j)%h)
          !---------------artificial viscosity-----------------
          rij_hat=(s(i)%x-s(j)%x) !unit vector between particles
          rij_hat=rij_hat/sqrt(dot_product(rij_hat,rij_hat)) 
          rho_bar=0.5*(s(i)%rho+s(j)%rho) !mean of density
          !get \nu_{sig}=c_i+c_j-v_{ij} \cdot \hat{r}_{ij}
          nu_sig=sqrt(SPH_gamma*s(i)%P/s(i)%rho)+sqrt(SPH_gamma*s(j)%P/s(j)%rho)-dot_product(s(i)%u-s(j)%u,rij_hat)
          s(i)%a=s(i)%a+nu_sig*(s(j)%m/rho_bar)*dot_product(s(i)%u-s(j)%u,rij_hat)*sph_grad_W(s(i)%x,s(j)%x,s(i)%h)/2.
          s(i)%a=s(i)%a+nu_sig*(s(j)%m/rho_bar)*dot_product(s(i)%u-s(j)%u,rij_hat)*sph_grad_W(s(i)%x,s(j)%x,s(j)%h)/2.
          !----------------compressibility------------
          s(i)%divu=s(i)%divu+(s(j)%m/s(i)%rho)*dot_product((s(i)%u-s(j)%u),sph_grad_W(s(i)%x,s(j)%x,s(i)%h))      
        end if
        !---------deallocate neighbour list-----
        s_NN(i)%previous => s_NN(i)%current
        s_NN(i)%current => s_NN(i)%current%next
        deallocate( s_NN(i)%previous )
      end do
    end do
  end subroutine
  !************************************************************
  !>get the acceleration of a particle due to nonlocal forces
  subroutine sph_get_non_loc
    implicit none
    integer :: i,j
    real :: a(3)
    real :: dist
    if (SPH_G>0.) then
      do i=1, SPH_count
        if (SPH_theta>0) then
          !Ge tree approximation
          a=0. !0 this before we start the walk
          call SPH_gravity_tree_walk(i,stree,a) !sph_tree.mod
          s(i)%a=s(i)%a+a
        else
          !Brute force N^2 
          do j=1, SPH_count
            if (i==j) cycle
            dist=dist_gen(s(i)%x,s(j)%x)
            s(i)%a=s(i)%a-(s(i)%x-s(j)%x)*SPH_G*s(j)%m*sph_phi(dist,s(i)%h)/(2.*dist)
            s(i)%a=s(i)%a-(s(i)%x-s(j)%x)*SPH_G*s(j)%m*sph_phi(dist,s(j)%h)/(2.*dist)
          end do
        end if
      end do
    end if
  end subroutine
  !********************************************************
  !generate linked list of nearest neighbours
  subroutine SPH_neighbour_list_notree
    integer :: i, j
    integer :: counter
    real :: dist
    do i=1, SPH_count
      counter=0
      nullify(s_NN(i)%list)
      allocate(s_NN(i)%list)
      nullify(s_NN(i)%list%next)
      s_NN(i)%current => s_NN(i)%list
      do j=1, SPH_count
        dist=dist_gen(s(i)%x,s(j)%x)
        if (dist<2*s(i)%h) then
          counter=counter+1
          if (counter==1) then
            !the first neighbour is added to the top of the list
            s_NN(i)%current%i=j
          else
            !subsequent neighbours linked to this 
            allocate( s_NN(i)%current%next )
            nullify( s_NN(i)%current%next%next )
            s_NN(i)%current%next%i=j
            s_NN(i)%current => s_NN(i)%current%next
          end if
        end if
      end do
      s(i)%ncount=counter
    end do
  end subroutine
  !************************************************************
  !>timestep-check SPH
  subroutine time_check_SPH
    implicit none
    integer :: i
    real :: sphdt,sphdt1,sphdt2, minsphdt
    minsphdt=dt
    do i=1, SPH_count
      if (s(i)%rho<epsilon(0.)) cycle
      sphdt2=0.2*s(i)%h/sqrt(SPH_gamma*s(i)%P/s(i)%rho)
      sphdt1=sqrt(s(i)%h/(sqrt(dot_product(s(i)%a,s(i)%a))+0.0001))
      sphdt=min(sphdt1,sphdt2)
      if (sphdt<minsphdt) then
        minsphdt=sphdt
      end if
    end do
    dt=minsphdt
  end subroutine
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
      write(78,*) '%-var-----t--------dt-------------maxu----------maxdu &
                  --------urms------minh----max_neigh-----fbar--maxdivu'
      if (particles_only) write(*,*) '%-var-----t------dt-----------maxu &
                          --------maxdu------urms------minh---max_neigh----fbar--maxdivu'
    end if
    write(78,'(i5.3,f9.2,e12.4,f15.8,f15.8,f15.8,f15.8,i12.4,f12.8,f12.8)') &
    itime/shots,t,dt,SPH_maxu,SPH_maxdu,SPH_urms,minval(s%h),maxval(s%ncount),sum(s%f)/SPH_count,maxval(s%divu)
    if (particles_only) then
      write(*,'(i5.3,f9.2,e12.4,f12.4,f10.2,f10.5,f10.4,i10.4,f10.4,f10.4)') &
      itime/shots,t,dt,SPH_maxu,SPH_maxdu,SPH_urms,minval(s%h),maxval(s%ncount),sum(s%f)/SPH_count,maxval(s%divu)
    end if
    close(78)
  end subroutine
  !************************************************************
  !>setup the particles in as set by SPH_init in run.in
  subroutine setup_SPH
    implicit none
    integer :: i, j
    real :: rs, rthet, rphi !for some intial conditions
    !check that the mass of the particles has been set in run.in
    if (SPH_mass<epsilon(0.)) call fatal_error('setup_SPH',&
                              'mass of SPH particles not set')
    allocate(s(SPH_count)) !allocate main SPH particle array
    allocate(s_NN(SPH_count)) !allocate main SPH particle array
    write(*,'(a)') ' --------------------SPH--------------------' 
    write(*,'(a,i7.1,a)') ' initialising ', SPH_count, ', particles '
    write(*,'(a,a,a,e14.7)') ' initial setup ', trim(SPH_init), ', with mass ', SPH_mass
    write(*,'(a,f10.5)') ' adiabatic index, gamma= ', SPH_gamma
    write(*,'(a,f10.5)') ' gravitational constant= ', SPH_G
    if (SPH_theta>epsilon(0.)) write(*,'(a,f10.5)') ' using tree, opening angle= ', SPH_theta
    select case(SPH_init)
      case('sphere')
        !particles in a uniform sphere - http://en.wikipedia.org/wiki/N-sphere
        write(*,'(a,f10.5)') ' radius of sphere= ', SPH_init_r*box_size
        do i=1, SPH_count
          s(i)%x(1)=rnorm(0.,1.) ; s(i)%x(2)=rnorm(0.,1.) ; s(i)%x(3)=rnorm(0.,1.)
          call random_number(rs)
          s(i)%x=box_size*SPH_init_r*0.5*(rs**0.333333333333)*s(i)%x/sqrt(dot_product(s(i)%x,s(i)%x))
          s(i)%u=0. ; s(i)%u1=0. ; s(i)%u2=0. !0 the velocity fields
          s(i)%a=0. ; s(i)%a1=0. ; s(i)%a2=0. !0 the acc fields
          s(i)%i=i !needed for neighbour list
        end do
      case('sphere+noise')
        !particles in uniform sphere with noise in velocity field
        write(*,'(a,f10.5)') ' radius of sphere= ', SPH_init_r*box_size
        do i=1, SPH_count
          s(i)%x(1)=rnorm(0.,1.) ; s(i)%x(2)=rnorm(0.,1.) ; s(i)%x(3)=rnorm(0.,1.)
          call random_number(rs)
          s(i)%x=box_size*SPH_init_r*0.5*(rs**0.333333333333)*s(i)%x/sqrt(dot_product(s(i)%x,s(i)%x))
          call random_number(s(i)%u) ; s(i)%u=2*(2*s(i)%u-1.)
          s(i)%u1=s(i)%u ; s(i)%u2=s(i)%u !set the velocity fields
          s(i)%a=0. ; s(i)%a1=0. ; s(i)%a2=0. !0 the acc fields
          s(i)%i=i !needed for neighbour list
        end do
      case('random_sphere')
        !particles in a random sphere
        write(*,'(a,f10.5)') ' radius of sphere= ', SPH_init_r*box_size
        do i=1, SPH_count
          call random_number(rs)
          call random_number(rthet) ; call random_number(rphi)
          rs=rs*box_size*SPH_init_r ; rthet=rs*(2*rthet-1.) ; rphi=rphi*2*pi
          s(i)%x(1)=sqrt(rs**2-rthet**2)*cos(rphi)
          s(i)%x(2)=sqrt(rs**2-rthet**2)*sin(rphi)
          s(i)%x(3)=rthet
          s(i)%u=0. ; s(i)%u1=0. ; s(i)%u2=0. !0 the velocity fields
          s(i)%a=0. ; s(i)%a1=0. ; s(i)%a2=0. !0 the acc fields
          s(i)%i=i !needed for neighbour list
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
          s(i)%i=i !needed for neighbour list
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
          s(i)%i=i !needed for neighbour list
        end do
      case('centre-plane')
        !particles in y-z plane at centre of box
        do i=1, SPH_count
          call random_number(s(i)%x(2)) ; call random_number(s(i)%x(3))
          s(i)%x(1)=0.
          s(i)%x(2)=box_size*s(i)%x(2)-box_size/2.
          s(i)%x(3)=box_size*s(i)%x(3)-box_size/2.
          s(i)%u=0. ; s(i)%u1=0. ; s(i)%u2=0. !0 the velocity fields
          s(i)%a=0. ; s(i)%a1=0. ; s(i)%a2=0. !0 the acc fields
          s(i)%i=i !needed for neighbour list
        end do
      case('figure8')
        if (SPH_count/=3) call fatal_error('SPH','SPH_count must be 3')
        s(1)%x(1)=0.9700000436
        s(2)%x(1)=-0.9700000436
        s(1)%x(2)=-0.2430875
        s(2)%x(2)=0.2430875
        s(1)%x(3)=0.
        s(2)%x(3)=0.
        s(3)%x=0.
        do i=1,SPH_count
          s(i)%u=0. ; s(i)%u1=0. ; s(i)%u2=0. !0 the velocity fields
          s(i)%a=0. ; s(i)%a1=0. ; s(i)%a2=0. !0 the acc fields
          s(i)%i=i !needed for neighbour list
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
          s(i)%i=i !needed for neighbour list
        end do
      case('loop')
        !this initial condition is set to interface with the filamemt
        !rs is based on filament separation
        rs=(0.75*SPH_count*delta)/(2*pi) !75% of potential size    
        do i=1, SPH_count      
          s(i)%x(1)=rs*sin(pi*real(2*i-1)/SPH_count)
          s(i)%x(2)=rs*cos(pi*real(2*i-1)/SPH_count)
          s(i)%x(3)=0.
          s(i)%u=0. ; s(i)%u1=0. ; s(i)%u2=0. !0 the velocity fields
          s(i)%a=0. ; s(i)%a1=0. ; s(i)%a2=0. !0 the acc fields
          s(i)%i=i !needed for neighbour list
        end do
      case default
        call fatal_error('setup_SPH','SPH_init not set to available value')
    end select
    s(:)%m=SPH_mass !set the masses of the particles
    !make an initial guess for the smoothing length
    s(:)%h=box_size/10.
    !improve the guess for h
    !first build the tree if needed
    if (SPH_theta>0) call create_SPH_tree !sph_tree.mod
    !now iteratively improve smoothing length
    do j=1, 5 !do this twenty times
      call sph_neighbour_loop
    end do
    if (SPH_theta>0) then
      call SPH_empty_tree(stree) !empty the tree to avoid a memory leak
      deallocate(stree%parray) ; deallocate(stree) ; nullify(stree)
    end if
    !finally set up SPH_mesh
    if (SPH_mesh_size>0) then
      call setup_SPH_mesh !init.mod
    else
      write(*,*) 'no SPH mesh allocated'
    end if
    nf_compressible=.true.
    !open(unit=23,file='data/sph_correction.log',status='replace')
    !do j=1, SPH_count
    !  write(23,*) s(j)%
    !end do 
    !close(23)
  end subroutine
  !****************************************************************
  subroutine setup_SPH_mesh
    implicit none
    integer :: i,j,k
    real :: res
    if (mod(SPH_mesh_size,2)/=0) then
      call fatal_error('sph.mod:setup_mesh', &
      'mesh size must be a multiple of 2')
    end if
    if (SPH_mesh_size<16) call warning_message('sph.mod:setup_mesh','warning mesh size is small')
    res=real(box_size)/SPH_mesh_size
    allocate(SPH_mesh(SPH_mesh_size,SPH_mesh_size,SPH_mesh_size))
    do k=1, SPH_mesh_size
      do j=1, SPH_mesh_size
        do i=1, SPH_mesh_size
          SPH_mesh(k,j,i)%x(1)=res*real(2*i-1)/2.-(box_size/2.)
          SPH_mesh(k,j,i)%x(2)=res*real(2*j-1)/2.-(box_size/2.)
          SPH_mesh(k,j,i)%x(3)=res*real(2*k-1)/2.-(box_size/2.)
          SPH_mesh(k,j,i)%rho=0. !clear the density slot
        end do
      end do
    end do
    write(*,*) 'allocated SPH mesh'
  end subroutine
  !****************************************************************
  subroutine mesh_print_SPH(filenumber)
    implicit none
    character (len=40) :: print_file
    integer, intent(IN) :: filenumber
    integer :: i,j,k,si
    real :: dist
    !first set the elements of the mesh
    SPH_mesh%rho=0.
    do k=1, SPH_mesh_size
      do j=1, SPH_mesh_size
        do i=1, SPH_mesh_size
          do si=1, SPH_count
            dist=dist_gen(s(si)%x,SPH_mesh(k,j,i)%x)
            SPH_mesh(k,j,i)%rho=SPH_mesh(k,j,i)%rho+s(si)%m*sph_W(dist,s(si)%h)
          end do
        end do
      end do
    end do
    write(unit=print_file,fmt="(a,i3.3,a)")"./data/SPH_mesh",filenumber,".dat"
    open(unit=92,file=print_file,form='unformatted',status='replace',access='stream')
      write(92) t
      write(92) SPH_mesh(SPH_mesh_size/2,SPH_mesh_size/2,1:SPH_mesh_size)%x(1)
      write(92) SPH_mesh%rho
    close(92)
    if (vapor_print) then
      write(unit=print_file,fmt="(a,i3.3,a)")"./data/SPH_vapor",filenumber,".dat"
      open(unit=92,file=print_file,form='unformatted',status='replace',access='stream')
        write(92) SPH_mesh%rho
      close(92)
    end if
  end subroutine
end module
