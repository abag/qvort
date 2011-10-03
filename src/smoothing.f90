!>All routines used to smooth the field onto the mesh, using tree methods.
!>At present smoothing the filaments onto a mesh can only be done with a
!>gaussian kernal:-
!>\f[ \hat{G}(\mathbf{x}) = \int (2\pi \sigma^2)^{-3/2} e^{-\frac{|\mathbf{x}-\mathbf{x}'|}{2\sigma^2}} G(\mathbf{x}') d\mathbf{x}'\f]
module smoothing
  use cdata 
  use general
  use tree
  implicit none
  !> mesh points on the smoothed mesh
  !>@param x the position of the mesh point
  !>@param w the smoothed field either vorticity 
  type smoothing_grid
     private
     real :: x(3)
     real :: w(3)
  end type
  !use the abbreviation sm (smoothed mesh)
  !>resolution of the smoothed mesh
  real, private :: sm_res 
  !>the smoothing length, set in terms of \f$\delta\f$ in run.in
  real, private :: sm_sigma 
  !>the array to calculate smoothed fields on  
  type(smoothing_grid), allocatable, private :: sm(:,:,:)
  contains
  !************************************************************
  !>set up the smoothing mesh, if set in run.in
  subroutine setup_smoothing_mesh
    implicit none
    integer :: i, j, k
    write(*,'(a)') ' -------------------------SMOOTHING----------------------------' 
    sm_sigma=smoothing_length*delta
    if (tree_theta>0) then
      write(*,'(a,i4.2,a)') ' creating a ', sm_size,' ^3 mesh to smooth \omega/B field onto'
      if (smoothing_interspace) then
        write(*,'(a)') ' using adaptive smoothing length based on inter-vortex spacing'
      else
        write(*,'(a,f8.5,a)') ' smoothing length is ', sm_sigma
      end if
    else
      call fatal_error('smoothing.mod','tree theta must be +ve to use smoothing')
    end if
    if (sm_size<16) call warning_message('smoothing.mod','sm_size<16 results will be poor')
    !print dimensions to file for matlab
    open(unit=77,file='./data/sm_dims.log',status='replace')
      write(77,*) sm_size
    close(77)
    allocate(sm(sm_size,sm_size,sm_size))
    sm_res=(real(box_size)/sm_size)
    do k=1, sm_size  ; do j=1, sm_size ; do i=1, sm_size
      sm(k,j,i)%x(1)=sm_res*real(2*i-1)/2.-(box_size/2.)
      sm(k,j,i)%x(2)=sm_res*real(2*j-1)/2.-(box_size/2.)
      sm(k,j,i)%x(3)=sm_res*real(2*k-1)/2.-(box_size/2.)
      sm(k,j,i)%w(:)=0. !0 the vorticity array
    end do ; end do ; end do
  end subroutine
  !**********************************************************************
  !>print the smoothed mesh to a binary file for plotting with matlab/paraview
  subroutine print_smooth_mesh(filenumber)
    implicit none
    integer, intent(IN) :: filenumber
    real, allocatable :: vapor_array(:,:,:)
    character (len=50) :: print_file
    write(unit=print_file,fmt="(a,i3.3,a)")"./data/smoothed_field",filenumber,".dat"
    open(unit=92,file=print_file,form='unformatted',status='replace',access='stream')
      write(92) t
      write(92) sm(sm_size/2,sm_size/2,1:sm_size)%x(1)
      write(92) sm(:,:,:)%w(1)
      write(92) sm(:,:,:)%w(2)
      write(92) sm(:,:,:)%w(3)
    close(92)
    if (vapor_print) then
      allocate(vapor_array(sm_size, sm_size, sm_size))
      vapor_array(:,:,:)=sqrt(sm(:,:,:)%w(1)**2+&
                              sm(:,:,:)%w(2)**2+&
                              sm(:,:,:)%w(3)**2)
      write(unit=print_file,fmt="(a,i3.3,a)")"./data/vap_smooth",filenumber,".dat"
      open(unit=93,file=print_file,form='unformatted',status='replace',access='stream')
        write(93) vapor_array
      close(93)
      deallocate(vapor_array) 
    end if
  end subroutine
  !************************************************************
  !>get the smoothed field on the smoothed mesh 
  subroutine get_smoothed_field
    implicit none
    integer :: i, j, k
    integer :: peri, perj, perk !used to loop in periodic cases
    real :: w(3)
    !define smoothing length 
    if (smoothing_interspace) then
      sm_sigma=sqrt((box_size**3)/total_length)
      !smoothing lenght can be checked against ts.m
      open(unit=78,file='data/smoothing_length.log',position='append')
        write(78,*) t, sm_sigma
      close(78)
    end if
    do k=1, sm_size  ; do j=1, sm_size ; do i=1, sm_size
      w=0. !0 this before each call
      call tree_smooth(sm(k,j,i)%x,vtree,(/0.,0.,0./),w)
      if (periodic_bc) then
        !we must shift the mesh in all 3 directions
        do peri=-1,1 ; do perj=-1,1 ; do perk=-1,1
          if (peri==0.and.perj==0.and.perk==0) cycle
          call tree_smooth(sm(k,j,i)%x,vtree,(/peri*box_size,perj*box_size,perk*box_size/),w)
        end do ; end do ;end do
      end if
      sm(k,j,i)%w=w
    end do ; end do ; end do
  end subroutine
  !**********************************************************
  !>given a specific position x, get the smoothed field wsmooth
  recursive subroutine tree_smooth(x,vtree,shift,wsmooth)
    implicit none
    real, intent(IN) :: x(3) !the position we want the velocity at
    real, intent(IN) :: shift(3) !shift tree (periodicity)
    type (node), pointer :: vtree !the tree
    real :: wsmooth(3) !smoothed vorticity field at point x
    real :: smoothing_factor
    real :: dist, dist2, theta !helper variables
    if (vtree%pcount==0) return !empty box no use
    !work out distances opening angles etc. here
    dist2=(x(1)-(vtree%centx+shift(1)))**2+&
         (x(2)-(vtree%centy+shift(2)))**2+&
         (x(3)-(vtree%centz+shift(3)))**2
    dist=sqrt(dist2)
    theta=vtree%width/dist !the most simple way we can improve this
    if (vtree%pcount==1.or.theta<tree_theta) then
      !use the contribution of this cell
      smoothing_factor=quant_circ*(1./((2.*pi*(sm_sigma**2))**(3/2)))*exp(-dist2/(2.*(sm_sigma**2)))
      wsmooth=wsmooth+vtree%circ*smoothing_factor
    else
      !open the box up and use the child cells
      call tree_smooth(x,vtree%fbl,shift,wsmooth) ; call tree_smooth(x,vtree%fbr,shift,wsmooth)
      call tree_smooth(x,vtree%ftl,shift,wsmooth) ; call tree_smooth(x,vtree%ftr,shift,wsmooth)
      call tree_smooth(x,vtree%bbl,shift,wsmooth) ; call tree_smooth(x,vtree%bbr,shift,wsmooth)
      call tree_smooth(x,vtree%btl,shift,wsmooth) ; call tree_smooth(x,vtree%btr,shift,wsmooth)
    end if 
  end subroutine
end module
