!> diagnostic routines purely for vortex filament, particle diagnostics contained within quasip.mod
module diagnostic
  use cdata
  use general
  use topology
  contains
  !*************************************************
  !>get the maximum velocity and change of velocity
  !>on the filament
  subroutine velocity_info()
    implicit none
    real, allocatable :: uinfo(:,:)
    allocate(uinfo(pcount,2))
    uinfo(:,1)=sqrt(f(:)%u(1)**2+f(:)%u(2)**2+f(:)%u(3)**2)
    uinfo(:,2)=sqrt((f(:)%u1(1)-f(:)%u2(1))**2+&
                    (f(:)%u1(2)-f(:)%u2(2))**2+&
                    (f(:)%u1(3)-f(:)%u2(3))**2)
    maxu=maxval(uinfo(:,1)) ; maxdu=maxval(uinfo(:,2))
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
  !>=\Gamma \oint_{\cal L} \mathbf{u} \cdot \mathbf{s} \times \mathbf{s}' d\xi \f]
  subroutine energy_info()
    implicit none
    real :: sdot(3), sdoti(3)
    integer :: infront
    integer :: i
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
      x(1)=((2.*i-1)/(2.*one_dim))*box_size-box_size/2.
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
      write(32,*) x(1), u_sup, u_norm 
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
    real :: u_norm(3), u_sup(3)=0., x(3)=0.
    integer :: i, j, peri, perj, perk
    character (len=40) :: print_file
    if (two_dim<1) return
    write(unit=print_file,fmt="(a,i4.4,a)")"./data/vel_slice_2D",filenumber,".dat"
    open(unit=32,file=print_file,status='replace',form='unformatted',access='stream')
    do i=1, two_dim
      do j=1, two_dim
        x(1)=((2.*i-1)/(2.*two_dim))*box_size-box_size/2.
        x(2)=((2.*j-1)/(2.*two_dim))*box_size-box_size/2.
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
        write(32) x(1), x(2), u_sup, u_norm 
      end do
    end do
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
    do i=1, pcount
      if (f(i)%infront==0) then
        curvi(i)=0. !check for 'empty' particles
      else
        curvi(i)=curvature(i) !general.mod
      end if
    end do
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
