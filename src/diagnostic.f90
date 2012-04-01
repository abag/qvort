!> diagnostic routines purely for vortex filament, particle diagnostics contained within quasip.mod
module diagnostic
  use cdata
  use general
  use topology
  contains
  !>dummy routine to call all diagnostic routines
  subroutine calculate_diagnostics()
    implicit none
    if (mod(itime, 5)==0) then
      if (closest_distance) call get_min_distance !diagnostics.mod    
    end if
    if (mod(itime, shots)==0) then
      call velocity_info !diagnostics.mod
      call curv_info !diagnostics.mod
      if (energy_inf) call energy_info !diagnostics.mod
      if (topo_inf) call get_topo_info !diagnostics.mod
      if (torsion_hist) call get_torsion_hist !diagnostics.mod      
      if (particle_plane_inf) call get_particle_plane_info !diagnostics.mod
      if (anisotropy_params) call get_anisotropy_info !diagnostics.mod
      if (full_loop_counter) call get_full_loop_count!diagnostics.mod
      if (mod(itime, mesh_shots)==0) then
        if (boxed_vorticity) call get_boxed_vorticity !diganostics.mod
        if (sep_inf) call get_sep_inf !diganostics.mod
        if (line_sep_inf) call get_line_sep_inf !diganostics.mod
        if (recon_time_info) call print_recon_time_info!diagnostics.mod
      end if 
    end if
  end subroutine
  !*************************************************
  !> Get the minimum spearation between vortices
  subroutine get_min_distance()
    implicit none
    open(unit=72,file='data/min_dist.log',position='append')
      write(72,*) t, minval(f(:)%closestd,mask=f(:)%infront>0),&
      -maxval(f(:)%x(3), mask=((f(:)%x(3)<0).and.(f(:)%x(1)>0.).and.(f(:)%x(1)<0.1).and.f(:)%infront>0))& 
      +minval(f(:)%x(3), mask=((f(:)%x(3)>0).and.(f(:)%x(1)>0.).and.(f(:)%x(1)<0.1).and.f(:)%infront>0))
    close(72)
  end subroutine
  !*************************************************
  !> The routine to calculate number of loops with their sizes and dumps them 
  !>to file to plot histograms.
  subroutine get_full_loop_count()
    implicit none
    real,allocatable :: line(:,:,:)
    integer :: line_count
    integer ::  next, next_old
    integer,allocatable :: counter(:)
    real,allocatable :: counter_rad(:)
    integer :: i, j, l, m
    logical :: unique
    character (len=30) :: line_file
    allocate(counter(ceiling(pcount/5.)))
    allocate(counter_rad(ceiling(pcount/5.)))
    allocate(line(ceiling(pcount/5.),pcount,4))
    counter=0 ; counter_rad=0. ; line_count=0
    !create file to print to 
    if (mod(itime,mesh_shots)==0) then
      write(unit=line_file,fmt="(a,i3.3,a)")'./data/loop_size',itime/mesh_shots,".log"
      open(unit=97,file=line_file,action="write",position="append")
    end if
    do i=1, pcount !determine starting position
      if (f(i)%infront/=0) then
        next=i
        exit
      end if
    end do
    next_old=next
    do l=1, 1000 !limited to 500 loops at present - should be do while
      do i=1, pcount
        line(l,i,1)=f(next)%x(1)
        line(l,i,2)=f(next)%x(2)
        line(l,i,3)=f(next)%x(3)
        line(l,i,4)=next
        counter_rad(l)=counter_rad(l)+dist_gen(f(next)%x,f(f(next)%infront)%x)
        next=f(next)%infront
        counter(l)=counter(l)+1
        if (next==next_old) exit
        if (next==0) exit
      end do
      ! make line a loop
      line_count=line_count+1
      if (mod(itime,mesh_shots)==0) then
        write(97,*) line_count, counter(l), counter_rad(l)
      end if
      !check size of line
      if (sum(counter)<count(mask=f(:)%infront>0)) then
        !need new staring point
        do i=1, pcount
          unique=.true.
          do m=1, l
          do j=1, counter(m)
            if (i==int(line(m,j,4)).or.f(i)%infront==0) then
              unique=.false.
            end if
          end do
          end do
          if (unique) then
            next=i
            next_old=next
            exit
          end if
        end do
      else
        exit
      end if 
    end do
    if (mod(itime,mesh_shots)==0) then
      close(97)
    end if
    open(unit=72,file='data/loop_counter.log',position='append')
      write(72,*) t, line_count
    close(72)
    deallocate(counter,counter_rad,line)
  end subroutine
  !*************************************************
  !>get the maximum velocity and change of velocity
  !>on the filament
  subroutine velocity_info()
    implicit none
    real, allocatable :: uinfo(:,:)
    allocate(uinfo(pcount,5))
    uinfo(:,1)=sqrt(f(:)%u(1)**2+f(:)%u(2)**2+f(:)%u(3)**2)
    uinfo(:,2)=sqrt((f(:)%u1(1)-f(:)%u2(1))**2+&
                    (f(:)%u1(2)-f(:)%u2(2))**2+&
                    (f(:)%u1(3)-f(:)%u2(3))**2)
    maxu=maxval(uinfo(:,1)) ; maxdu=maxval(uinfo(:,2))
    uinfo(:,3)=f(:)%u(1) ; uinfo(:,4)=f(:)%u(2) ; uinfo(:,5)=f(:)%u(3)
    !rather than use pcount more correct to use count(mask=f(:)%infront>0)
    !---------------------------------------------------------------------------
    open(unit=34,file='./data/basic_velocity_info.log',position='append')
      write(34,'(7e13.4)') t, maxval(uinfo(:,1)),&
                              sum(uinfo(:,1))/count(mask=f(:)%infront>0),&
                              minval(uinfo(:,1),mask=uinfo(:,1)>0),&
                              sum(uinfo(:,3))/count(mask=f(:)%infront>0),&
                              sum(uinfo(:,4))/count(mask=f(:)%infront>0),&
                              sum(uinfo(:,5))/count(mask=f(:)%infront>0)
    close(34)
    if (initf=='macro_ring') then
      !additional diagnostics for macro_ring initial_condition
      !---------------------------------------------------------------------------
      open(unit=34,file='./data/centre_of_vorticity_info.log',position='append')
        write(34,'(7e13.4)') t, cov%x, cov%u
      close(34)
      !---------------------------------------------------------------------------
      open(unit=34,file='./data/R_spread_info_1.log',position='append')
        write(34,'(18e13.4)') t, mrr1%x_max, mrr1%x_min, mrr1%x_spread, mrr1%r, mrr1%a, mrr1%r_u,   mrr1%a_u, mrr1%ra
      close(34)
      open(unit=34,file='./data/R_spread_info_2.log',position='append')
        write(34,'(26e13.4)') t, avg_r, avg_a, avg_r_u, avg_a_u, avg_ra
      close(34)        
    end if   
    deallocate(uinfo)
  end subroutine
  !*************************************************
  !>get the boxed vorticity by summing circulation vectors in
  !> a mesh
  subroutine get_boxed_vorticity()
    implicit none
    real, allocatable :: vort_mesh(:,:,:,:)
    real, allocatable :: vmeshx(:)
    real :: organised_length
    integer :: i,j,k !for looping over mesh
    integer :: vi !for looping over filament
    allocate(vort_mesh(boxed_vorticity_size,boxed_vorticity_size,boxed_vorticity_size,3))
    allocate(vmeshx(boxed_vorticity_size+1))
    !--------------------set the dimensions of the box-----------------
    do i=1,boxed_vorticity_size+1
      vmeshx(i)=-box_size/2.+box_size*(i-1)/boxed_vorticity_size;
    end do
    !------------------now loop over mesh and set the vorticity---------------
    vort_mesh=0. !0 this to begin with 
    do i=1,boxed_vorticity_size ; do j=1, boxed_vorticity_size ; do k=1, boxed_vorticity_size
      do vi=1,pcount
        if (f(vi)%infront==0) cycle !check for 'empty' particles
        !does particle lie within meshpoint?
        if ((f(vi)%x(1)>vmeshx(i)).and.(f(vi)%x(1)<=vmeshx(i+1))) then
          if ((f(vi)%x(2)>vmeshx(j)).and.(f(vi)%x(2)<=vmeshx(j+1))) then
            if ((f(vi)%x(3)>vmeshx(k)).and.(f(vi)%x(3)<=vmeshx(k+1))) then
              vort_mesh(k,j,i,:)=vort_mesh(k,j,i,:)+(f(vi)%ghosti-f(vi)%x);
            end if
          end if
        end if
      end do 
    end do ; end do ; end do
    !now sum the vorticity mesh - must be square
    organised_length=sum(sqrt(vort_mesh(:,:,:,1)**2+vort_mesh(:,:,:,3)**2+vort_mesh(:,:,:,3)**2))
    !print this length to file
    open(unit=72,file='data/organised_line_length.log',position='append')
      write(72,*) t, organised_length
    close(72)
    deallocate(vort_mesh) ; deallocate(vmeshx)
  end subroutine
  !*************************************************
  !>get information on (an)isotropy of filaments
  subroutine get_anisotropy_info()
    implicit none
    real :: l_para, l_perp !parallel and perpendicular densitys
    real :: l_l !3rd measure
    real,dimension(3) :: r_para, r_perp !unit vectors for r parallel and perpendicular
    real, dimension(3) :: sdot, sddot !1st and 2nd derivative
    integer :: i,j !for looping
    !set r parallel and perpendicular
    r_para=(/1.,0.,0./) ; r_perp=(/0.,1./sqrt(2.),1./sqrt(2.)/)
    l_para=0. ; l_perp=0. ; l_l=0. !0 these
    do i=1, pcount
      if (f(i)%infront==0) cycle !check for 'empty' particles
      call get_deriv_1(i,sdot)
      call get_deriv_2(i,sddot)
      l_para=l_para+(1.-(dot_product(sdot,r_para)**2))*dist_gen(f(i)%x,f(i)%ghosti)
      l_perp=l_perp+(1.-(dot_product(sdot,r_perp)**2))*dist_gen(f(i)%x,f(i)%ghosti)
      l_l=l_l+dot_product(cross_product(sdot,sddot),r_para)
    end do
    !normalise by total line_length
    l_para=l_para/total_length ; l_perp=l_perp/total_length
    !normalise 3rd quantity
    l_l=l_l*(sqrt(box_size/total_length))**3
    open(unit=72,file='data/anisotropy.log',position='append')
      write(72,*) t, l_para, l_perp, l_l
    close(72)
  end subroutine
  !*************************************************
  !>calculate the energy of the vortex filament using
  !>the trapezium rule, not valid with periodic b.c.'s
  !>\f[ E=\frac{1}{2} \int_V \mathbf{u}^2 dV=
  !>=\frac{\Gamma}{2} \oint_{\cal L} \mathbf{u} \cdot \mathbf{s} \times \mathbf{s}' d\xi \f]
  subroutine energy_info()
    implicit none
    real :: sdot(3), sdoti(3)
    integer :: infront
    integer :: i
    !compute the integral
    energy=0.
    do i=1, pcount
      if (f(i)%infront==0) cycle !check for 'empty' particles
      infront=f(i)%infront
      call get_deriv_1(i,sdot)
      call get_deriv_1(infront,sdoti)
      energy=energy+0.5*dist_gen(f(i)%x,f(i)%ghosti)*&
      (dot_product(f(i)%u,cross_product(f(i)%x,sdot))+&
       dot_product(f(infront)%u,cross_product(f(infront)%x,sdoti)))
    end do
    !multiply by terms outside the integral
    energy=energy*(quant_circ/2.)
  end subroutine
  !*************************************************
  !>print to file difference in time between reconnections
  !>that specific points have been involved in
  subroutine print_recon_time_info()
    implicit none
    integer :: i
    open(unit=49,file='./data/recon_times.log',position='append')
    do i=1, pcount
      if (f(i)%infront==0) cycle !check for 'empty' particles
      if ((f(i)%t_recon(1)>epsilon(0.)).and.&
          (f(i)%t_recon(2)>epsilon(0.)).and.&
          (f(i)%t_recon(1)-f(i)%t_recon(2)>epsilon(0.))) then
        write(49,*) f(i)%t_recon(1)-f(i)%t_recon(2)
      end if
    end do
    close(49)
  end subroutine
  !*************************************************
  !> Call all topological routines from topology module
  subroutine get_topo_info()
    implicit none
    call get_linking_number !topology.mod
    call get_writhing_number !topology.mod
  end subroutine
  !*************************************************
  !> Get information on particles in the plane
  subroutine get_particle_plane_info()
    implicit none
    integer :: totalN
    integer :: i !for looping
    totalN=0
    do i=1, pcount
      if (f(i)%infront==0) cycle !empty particle
      !particles in the xy plane at z=0
      if (f(i)%x(3)*f(i)%ghosti(3)<0) totalN=totalN+1
      !particles in the yz plane at x=0
      if (f(i)%x(1)*f(i)%ghosti(1)<0) totalN=totalN+1
      !particles in the xz plane at y=0
      if (f(i)%x(2)*f(i)%ghosti(2)<0) totalN=totalN+1
    end do 
    open(unit=72,file='data/particle_plane.log',position='append')
      write(72,*) t, totalN
    close(72)
  end subroutine
  !*************************************************
  !> get normal/superfluid velocity for a 1D strip
  subroutine one_dim_vel(filenumber)
    use tree
    use timestep
    implicit none
    integer, intent(IN) :: filenumber
    real :: u_norm(3), u_sup(3)=0., x(3)=0.
    integer :: i, peri, perj, perk
    character (len=40) :: print_file
    if (one_dim<1) return
    write(unit=print_file,fmt="(a,i4.4,a)")"./data/vel_slice_1D",filenumber,".log"
    open(unit=32,file=print_file)
    do i=1, one_dim
      x=0.
      select case (one_dim_direction)
        case('x')
          x(1)=((2.*i-1)/(2.*one_dim))*box_size-box_size/2.
        case('y')
          x(2)=((2.*i-1)/(2.*one_dim))*box_size-box_size/2.
        case('z')
          x(1)=-0.005
          x(3)=((2.*i-1)/(2.*one_dim))*box_size-box_size/2.
      end select
      !superfluid velocity
      select case(velocity)
        case('BS')
          u_sup=0. !always 0 before making initial BS call
          call biot_savart_general(x,u_sup)
          if (periodic_bc) then
            !we must shift the mesh in all 3 directions, all 26 permutations needed!
            do peri=-1,1 ; do perj=-1,1 ; do perk=-1,1
              if (peri==0.and.perj==0.and.perk==0) cycle
              call biot_savart_general_shift(x,u_sup, &
              (/peri*box_size,perj*box_size,perk*box_size/)) !timestep.mod
            end do ; end do ;end do
          end if
        case('Tree')
          u_sup=0. !must be zeroed for all algorithms
          call tree_walk_general(x,vtree,(/0.,0.,0./),u_sup)
          if (periodic_bc) then
            !we must shift the mesh in all 3 directions, all 26 permutations needed!
            do peri=-1,1 ; do perj=-1,1 ; do perk=-1,1
              if (peri==0.and.perj==0.and.perk==0) cycle
              call tree_walk_general(x,vtree, &
                   (/peri*box_size,perj*box_size,perk*box_size/),u_sup) !tree.mod
            end do ; end do ;end do
          end if
      end select
      call get_normal_velocity(x,u_norm)
      select case (one_dim_direction)
        case('x')
          write(32,*) x(1), u_sup, u_norm 
        case('y')
          write(32,*) x(2), u_sup, u_norm 
        case('z')
          write(32,*) x(3), u_sup, u_norm 
      end select
    end do
    close(32)
  end subroutine
  !*************************************************
  !> get normal/superfluid velocity one the 1D lattice structure
  subroutine one_dim_lattice_vel(filenumber)
    use tree
    use timestep
    implicit none
    integer, intent(IN) :: filenumber
    integer :: i, j, peri, perj, perk
    character (len=40) :: print_file
    if (one_dim_lattice<1) return
    write(unit=print_file,fmt="(a,i4.4,a)")"./data/1D_lattice",filenumber,".dat"
    open(unit=32,file=print_file,status='replace',form='unformatted',access='stream')
    do j=1, one_dim_lattice_count
      do i=1, one_dim_lattice
        !superfluid velocity
        select case(velocity)
          case('BS')
            lat_mesh_1D(j,i)%u_sup=0.!always 0 before making initial BS call
            call biot_savart_general(lat_mesh_1D(j,i)%x,lat_mesh_1D(j,i)%u_sup)
            if (periodic_bc) then
              !we must shift the mesh in all 3 directions, all 26 permutations needed!
              do peri=-1,1 ; do perj=-1,1 ; do perk=-1,1
                if (peri==0.and.perj==0.and.perk==0) cycle
                call biot_savart_general_shift(lat_mesh_1D(j,i)%x,lat_mesh_1D(j,i)%u_sup, &
                (/peri*box_size,perj*box_size,perk*box_size/)) !timestep.mod
              end do ; end do ;end do
            end if
          case('Tree')
            lat_mesh_1D(j,i)%u_sup=0. !must be zeroed for all algorithms
            call tree_walk_general(lat_mesh_1D(j,i)%x,vtree,(/0.,0.,0./),lat_mesh_1D(j,i)%u_sup)
            if (periodic_bc) then
              !we must shift the mesh in all 3 directions, all 26 permutations needed!
              do peri=-1,1 ; do perj=-1,1 ; do perk=-1,1
                if (peri==0.and.perj==0.and.perk==0) cycle
                call tree_walk_general(lat_mesh_1D(j,i)%x,vtree, &
                     (/peri*box_size,perj*box_size,perk*box_size/),lat_mesh_1D(j,i)%u_sup) !tree.mod
              end do ; end do ;end do
            end if
        end select
        call get_normal_velocity(lat_mesh_1D(j,i)%x,lat_mesh_1D(j,i)%u_norm)
        write(32) lat_mesh_1D(j,i)%x,lat_mesh_1D(j,i)%u_sup,lat_mesh_1D(j,i)%u_norm
      end do
    end do
    close(32)
  end subroutine
  !*************************************************
  !> get normal/superfluid velocity for a 2D slice
  subroutine two_dim_vel(filenumber)
    use tree
    use timestep
    implicit none
    integer, intent(IN) :: filenumber
    integer :: i, j, peri, perj, perk
    character (len=40) :: print_file
    real :: u(3)
    if (two_dim<1) return  
    !$omp parallel do private(i,j,peri,perj,perk,u)
    do j=1, two_dim
      do i=1, two_dim
        !superfluid velocity
        select case(velocity)
          case('BS')
            mesh2D(j,i)%u_sup=0. !always 0 before making initial BS call
            call biot_savart_general(mesh2D(j,i)%x,mesh2D(j,i)%u_sup)
            if (periodic_bc) then
              !we must shift the tree-mesh in all 3 directions, all 26 permutations needed!
              do peri=-1,1 ; do perj=-1,1 ; do perk=-1,1
                if (peri==0.and.perj==0.and.perk==0) cycle
                call biot_savart_general_shift(mesh2D(j,i)%x,mesh2D(j,i)%u_sup, &
                (/peri*box_size,perj*box_size,perk*box_size/)) !timestep.mod
              end do ; end do ;end do
            end if
          case('Tree')
            mesh2D(j,i)%u_sup=0. !always 0 before making initial BS call
            call tree_walk_general(mesh2D(j,i)%x,vtree,(/0.,0.,0./),mesh2D(j,i)%u_sup)
            if (periodic_bc) then
              !we must shift the tree-mesh in all 3 directions, all 26 permutations needed!
              do peri=-1,1 ; do perj=-1,1 ; do perk=-1,1
                if (peri==0.and.perj==0.and.perk==0) cycle
                call tree_walk_general(mesh2D(j,i)%x,vtree, &
                     (/peri*box_size,perj*box_size,perk*box_size/),mesh2D(j,i)%u_sup) !tree.mod
              end do ; end do ;end do
            end if
        end select
        !normal fluid velocity
        !$OMP critical
        call get_normal_velocity(mesh2D(j,i)%x,u)
        !$OMP end critical
        mesh2D(j,i)%u_norm=u
      end do
    end do
    !$omp end parallel do
    write(unit=print_file,fmt="(a,i4.4,a)")"./data/vel_slice_2D",filenumber,".dat"
    open(unit=32,file=print_file,status='replace',form='unformatted',access='stream')
      write(32) mesh2D(1,1:two_dim)%x(1)
      write(32) mesh2D(1:two_dim,1:two_dim)%u_norm(1)
      write(32) mesh2D(1:two_dim,1:two_dim)%u_norm(2)
      write(32) mesh2D(1:two_dim,1:two_dim)%u_norm(3)
      write(32) mesh2D(1:two_dim,1:two_dim)%u_sup(1)
      write(32) mesh2D(1:two_dim,1:two_dim)%u_sup(2)
      write(32) mesh2D(1:two_dim,1:two_dim)%u_sup(3)
    close(32)
  end subroutine
  !*************************************************
  !>caculate the mean, min, max curvature of the filament
  !>if set in run.in will also bin the curvatures to plot a 
  !>histogram
  subroutine curv_info()
    use kernel_density
    implicit none
    real, allocatable :: curvi(:)
    integer :: i, j
    !------------histogram parameters below---------------------
    !warning - at present the matlab routine to plot the histogram 
    !is not adpative therefore if the number of bins is changed 
    !the script must also be changed - warning
    integer, parameter :: bin_num=10 !number of bins used in KDE
    real :: kdensity(bin_num), kdmesh(bin_num), bwidth !all for KDE
    allocate(curvi(pcount)) !allocate this array pcount size
    !$omp parallel do private(i)
    do i=1, pcount
      if (f(i)%infront==0) then
        curvi(i)=0. !check for 'empty' particles
      else
        curvi(i)=curvature(i) !general.mod
      end if
    end do
    !$omp end parallel do
    !compute the average/min/max of this array
    kappa_bar=sum(curvi)/count(mask=f(:)%infront>0)
    kappa_max=maxval(curvi)
    kappa_min=minval(curvi,mask=curvi>0)
    !experimental creation of a histogram using kernel density estimation
    if (curv_hist) then
      !set the mesh over which we calculate densities
      do j=1, bin_num
        kdmesh(j)=kappa_min+(j-1)*(kappa_max-kappa_min)/(bin_num-1)
      end do
      call qsort(curvi) !sort the data in ascending order
      call get_kernel_density(curvi,pcount,kdmesh,bin_num,kdensity,bwidth)
      open(unit=79,file='data/curv_pdf.log',position='append')
      write(79,*) '%----------------t=',t,'---------------'
      do j=1, bin_num
        write(79,*) kdmesh(j), kdensity(j)
      end do 
      close(79)
    end if
    deallocate(curvi) !deallocate helper array
  end subroutine
  !*************************************************
  !>caculate the mean, min, max torsion of the filament
  !>and will also bin the torsions to plot a histogram
  subroutine get_torsion_hist()
    use kernel_density
    implicit none
    real, allocatable :: torsioni(:)
    integer :: i, j
    real :: tors_min, tors_max, tors_bar
    real :: ddot_i(3), ddot_infront(3), ddot_behind(3)
    real :: sdot(3), sddot(3), sdddot(3)
    real :: disti, distb
    !------------histogram parameters below---------------------
    !warning - at present the matlab routine to plot the histogram 
    !is not adpative therefore if the number of bins is changed 
    !the script must also be changed - warning
    integer, parameter :: bin_num=10 !number of bins used in KDE
    real :: kdensity(bin_num), kdmesh(bin_num), bwidth !all for KDE
    allocate(torsioni(pcount)) !allocate this array pcount size
    do i=1, pcount
      if (f(i)%infront==0) then
        torsioni(i)=0. !check for 'empty' particles
      else
        call get_deriv_1(i,sdot) !derivatives.mod
        call get_deriv_2(i,ddot_i) ; sddot=ddot_i
        call get_deriv_2(f(i)%infront,ddot_infront) !derivatives.mod
        call get_deriv_2(f(i)%behind,ddot_behind) !derivatives.mod
        disti=dist_gen(f(i)%x,f(i)%ghosti) ; distb=dist_gen(f(i)%x,f(i)%ghostb)
        sdddot=distb*ddot_infront+ &
              (disti-distb)*ddot_i- &
               disti*ddot_behind
        sdddot=sdddot/(2.*distb*disti)
        torsioni(i)=dot_product(cross_product(sdot,sddot),sdddot)
        torsioni(i)=torsioni(i)/dot_product(cross_product(sdot,sddot),&
                                           cross_product(sdot,sddot))
      end if
    end do
    !compute the average/min/max of this array
    tors_bar=sum(torsioni)/count(mask=f(:)%infront>0)
    tors_max=maxval(torsioni)
    tors_min=minval(torsioni,mask=torsioni>0)
    open(unit=76,file='data/torsion.log',position='append')
      write(76,*) tors_bar, tors_min, tors_max
    close(76)
    !set the mesh over which we calculate densities
    do j=1, bin_num
      kdmesh(j)=tors_min+(j-1)*(tors_max-tors_min)/(bin_num-1)
    end do
    call qsort(torsioni) !sort the data in ascending order
    call get_kernel_density(torsioni,pcount,kdmesh,bin_num,kdensity,bwidth)
    open(unit=79,file='data/torsion_pdf.log',position='append')
    write(79,*) '%----------------t=',t,'---------------'
    do j=1, bin_num
      write(79,*) kdmesh(j), kdensity(j)
    end do 
    close(79)
    deallocate(torsioni) !deallocate helper array
  end subroutine
  !*************************************************
  !>caculate the mean, min, max of separation of the vortex points
  !>will also bin the separations to plot a histogram
  subroutine get_sep_inf()
    implicit none
    real, allocatable :: sep_array(:) !point seperations
    real :: sep_min, sep_max, sep_bar
    integer :: i, j, counter
    !------------histogram parameters below---------------------
    !warning - at present the matlab routine to plot the histogram 
    !is not adpative therefore if the number of bins is changed 
    !the script must also be changed - warning
    integer, parameter :: bin_num=20 !number of bins
    real :: bin_size, bins(bin_num) !the bins of the histograms
    real :: probability(bin_num) !the probability of that bin
    allocate(sep_array(pcount*(pcount-1))) !allocate this array pcount*(pcount-1) size
    sep_array=0. !0 initially
    counter=1
    do i=1, pcount
      if (f(i)%infront==0) cycle
      do j=i, pcount
        if (f(j)%infront==0) cycle
        sep_array(counter)=distf(i,j)
        counter=counter+1
      end do
    end do
    !compute the average/min/max of this array
    sep_bar=sum(sep_array)/count(mask=sep_array>0.)
    sep_max=maxval(sep_array)
    sep_min=minval(sep_array,mask=sep_array>0)
    !output this to file
    open(unit=72,file='data/sep_info.log',position='append')
      write(72,*) t, sep_bar,sep_max,sep_min
    close(72)
    !now create of a histogram

    !first we set the bins
    bin_size=(sep_max-sep_min)/bin_num
    bins=0.
    !now loop over sep_array and bins and create the histogram
    do i=1, pcount*(pcount-1)
      do j=1, bin_num
        if (sep_array(i)>sep_min+(j-1)*bin_size) then
          if (sep_array(i)<sep_min+j*bin_size) then
            bins(j)=bins(j)+1
          end if
        end if
      end do
    end do
    bins=bins/(count(mask=sep_array>0.)*bin_size) !normalise 
    open(unit=79,file='data/sep_pdf.log',position='append')
    write(79,*) '%----------------t=',t,'---------------'
    do j=1, bin_num
      write(79,*) sep_min+(j-1)*bin_size, bins(j) 
    end do 
    close(79)

    deallocate(sep_array) !deallocate helper array
  end subroutine
  !*************************************************
  !>caculate the mean, min, max of minimum separation of the vortex filaments
  !>will also bin the separations to plot a histogram
  subroutine get_line_sep_inf()
    use tree !needs tree module
    implicit none
    real :: line_sep_min, line_sep_max, line_sep_bar
    integer :: i, j
    !------------histogram parameters below---------------------
    !warning - at present the matlab routine to plot the histogram 
    !is not adpative therefore if the number of bins is changed 
    !the script must also be changed - warning
    integer, parameter :: bin_num=10 !number of bins
    real :: bin_size, bins(bin_num) !the bins of the histograms
    real :: probability(bin_num) !the probability of that bin
    !find the nearest point on another loop
    call pclose_tree_loop !tree.mod
    if (maxval(f(:)%closestd_loop,mask=(f(:)%infront>0))>2*box_size) then
      call warning_message('get_line_sep_info','loops are too diffuse to make this routine feasible')
    end if
    !compute the average/min/max of this array
    line_sep_min=minval(f(:)%closestd_loop,mask=f(:)%infront>0)
    line_sep_max=maxval(f(:)%closestd_loop,mask=((f(:)%infront>0).and.(f(:)%closestd_loop<2.*box_size)))
    line_sep_bar=sum(f(:)%closestd_loop,mask=((f(:)%infront>0).and.(f(:)%closestd_loop<2.*box_size)))&
                 /count(mask=(f(:)%infront>0).and.(f(:)%closestd_loop<2.*box_size))
    !output this to file
    open(unit=72,file='data/line_sep_info.log',position='append')
      write(72,*) t, line_sep_bar,line_sep_max,line_sep_min
    close(72)
    !now create of a histogram

    !first we set the bins
    bin_size=(line_sep_max-line_sep_min)/bin_num
    bins=0.
    !now loop over sep_array and bins and create the histogram
    do i=1, pcount
      do j=1, bin_num
        if (f(i)%closestd_loop>line_sep_min+(j-1)*bin_size) then
          if (f(i)%closestd_loop<line_sep_min+j*bin_size) then
            bins(j)=bins(j)+1
          end if
        end if
      end do
    end do
    !normalise the histogram
    bins=bins/(count(mask=(f(:)%infront>0).and.(f(:)%closestd_loop<2.*box_size))*bin_size)
    open(unit=79,file='data/line_sep_pdf.log',position='append')
    write(79,*) '%----------------t=',t,'---------------'
    do j=1, bin_num
      write(79,*) line_sep_min+(j-1)*bin_size, bins(j) 
    end do 
    close(79)
  end subroutine
  !*************************************************
  !>get the structure functions of the mesh velocities
  !>@warning not completed yet
  subroutine structure_function() !THIS NEEDS TESTING!!!!
    implicit none
    real, allocatable :: sfunction(:)
    integer :: i, j, k
    integer :: lag
    integer :: peri, perj, perk
    allocate(sfunction(mesh_size-1))
    sfunction=0 !empty this
    do k=1, mesh_size ; do j=1, mesh_size ; do i=1, mesh_size
      do lag=1, mesh_size-1
        if (i-lag<1) then        
          peri=mesh_size+i-lag
        else
          peri=i-lag
        end if
        if (j-lag<1) then        
          perj=mesh_size+j-lag
        else
          perj=j-lag
        end if
        if (k-lag<1) then        
          perk=mesh_size+k-lag
        else
          perk=k-lag
        end if
        sfunction(lag)=sfunction(lag)+(mesh(k,j,i)%u_sup(1)-mesh(perk,j,i)%u_sup(1))
        sfunction(lag)=sfunction(lag)+(mesh(k,j,i)%u_sup(1)-mesh(k,perj,i)%u_sup(1))
        sfunction(lag)=sfunction(lag)+(mesh(k,j,i)%u_sup(1)-mesh(k,j,peri)%u_sup(1))
      end do
    end do ; end do ; end do
    deallocate(sfunction)
  end subroutine
end module
