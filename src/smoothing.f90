module smoothing
  ! CONTAINS ALL THE ROUTINES USED TO SMOOTH THE FIELD ONTO A MESH
  ! USING TREE CODE METHODS - THIS IS EITHER VORTICITY OR MAGNETIC FIELD
  ! VORTICITY USED TO LOOK FOR BUNDLING OF FILAMENTS, MAGNETIC FIELD JxB
  use cdata 
  use general
  use tree
  implicit none
  type smoothing_grid
     private
     real :: x(3) !position
     real :: w(3) !vorticity/magnetic field 
  end type
  !use the abbreviation sm (smoothed mesh)
  integer, private, parameter :: sm_size=64 !to be used with below
  real, private :: sm_res !resolution of mesh
  real, private :: sm_sigma !smoothing length
  type(smoothing_grid), allocatable, private :: sm(:,:,:)
  contains
  !************************************************************
  subroutine setup_smoothing_mesh
    implicit none
    integer :: i, j, k
    sm_sigma=smoothing_length*delta
    if (tree_theta>0) then
      write(*,'(a,i4.2,a)') ' creating a ', sm_size,' ^3 mesh to smooth \omega/B field onto'
      write(*,'(a,f8.5,a)') ' smooting length is ', sm_sigma
    else
      call fatal_error('smoothing.mod','tree theta must be +ve to use smoothing')
    end if
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
      sm(k,j,i)%w(:)=0.
    end do ; end do ; end do
    !0 the vorticity/magnetic field array
  end subroutine
  !**********************************************************************
  subroutine print_smooth_mesh(filenumber)
    !print the mesh to a binary file
    implicit none
    integer, intent(IN) :: filenumber
    character (len=50) :: print_file
    write(unit=print_file,fmt="(a,i3.3,a)")"./data/smoothed_field",filenumber,".dat"
    open(unit=92,file=print_file,form='unformatted',status='replace',access='stream')
      write(92) t
      write(92) sm(sm_size/2,sm_size/2,1:sm_size)%x(1)
      write(92) sm(:,:,:)%w(1)
      write(92) sm(:,:,:)%w(2)
      write(92) sm(:,:,:)%w(3)
    close(92)
  end subroutine
  !************************************************************
  subroutine get_smoothed_field
    implicit none
    integer :: i, j, k
    real :: w(3)
    do k=1, sm_size  ; do j=1, sm_size ; do i=1, sm_size
      w=0. !0 this before each call
      !this is not truly periodic yet, need to loop over 27 directions
      call tree_smooth(sm(k,j,i)%x,vtree,(/0.,0.,0./),w)
      sm(k,j,i)%w=w
    end do ; end do ; end do
    !now need to print this to file
  end subroutine
  !**********************************************************
  recursive subroutine tree_smooth(x,vtree,shift,wsmooth)
    implicit none
    real, intent(IN) :: x(3) !the position we want the velocity at
    real, intent(IN) :: shift(3) !shift tree (periodicity)
    type (node), pointer :: vtree !the tree
    real :: wsmooth(3) !smoothed vorticity/magnetic field at point x
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
      smoothing_factor=(1./sqrt(2.*pi*(sm_sigma**2)))*exp(-dist2/(2.*(sm_sigma**2)))
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
