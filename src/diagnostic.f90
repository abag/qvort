!> diagnostic routines purely for vortex filament, if acting as a magnetic field diagnostics done within the module
!> mag, particle diagnostics contained within quasip.mod
module diagnostic
  use cdata
  use general
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
    write(unit=print_file,fmt="(a,i3.3,a)")"./data/vel_slice_1D",filenumber,".log"
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
      write(32,*) x(1), u_sup, u_norm 
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
    write(unit=print_file,fmt="(a,i3.3,a)")"./data/vel_slice_2D",filenumber,".dat"
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
    implicit none
    real, allocatable :: curvi(:)
    integer :: i, j
    !------------histogram parameters below---------------------
    !warning - at present the matlab routine to plot the histogram 
    !is not adpative therefore if the number of bins is changed 
    !the script must also be changed - warning
    integer, parameter :: bin_num=10 !number of bins
    real :: bin_size, bins(bin_num) !the bins of the histograms
    real :: probability(bin_num) !the probability of that bin
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
    !experimental creation of a histogram
    if (curv_hist) then
      !first we set the bins
      bin_size=(kappa_max-kappa_min)/bin_num
      bins=0.
      !now loop over particles and bins and create the histogram
      do i=1, pcount
        if (f(i)%infront==0) cycle
        do j=1, bin_num
          if (curvi(i)>kappa_min+(j-1)*bin_size) then
            if (curvi(i)<kappa_min+j*bin_size) then
              bins(j)=bins(j)+1
            end if
          end if
        end do
      end do
      bins=bins/(count(mask=f(:)%infront>0)*bin_size) !normalise 
      open(unit=79,file='data/curv_pdf.log',position='append')
      write(79,*) '%----------------t=',t,'---------------'
      do j=1, bin_num
        write(79,*) kappa_min+(j-1)*bin_size, bins(j) 
      end do 
      close(79)
    end if
    deallocate(curvi) !deallocate helper array
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
