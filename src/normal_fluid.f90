module normal_fluid
  !NORMAL FLUID COMPONENT IN THE EQUATION OF MOTION
  use cdata
  use general
  use ksmodel
    !parameters used in get_normal_fluid
    real, parameter, private :: vel_xflow=1.
    real, parameter, private :: abc_A=1., abc_B=1., abc_C=1.
    real, private :: abc_k
    !normal fluid mesh - all grid based normal velocities use this (e.g. Navier Stokes future)
    type norm_fluid_grid
      private
      real :: x(3) !position
      real :: u(3) !velocity 
      real :: v !scalar - compressible flow uses this
      real :: div !divergence - check div of vel field
    end type
    !use the abbreviation nfm (normal fluid mesh)
    integer, private, parameter :: nfm_size=64 !to be used with below
    real, private :: nfm_res, nfm_inv_res !resolution/inv resolution of normal fluid mesh
    type(norm_fluid_grid), allocatable, private :: nfm(:,:,:)
    contains 
    !************************************************************
    subroutine setup_normal_fluid
      !setup everything needed to use a normal fluid
      implicit none
      abc_k=2.*pi/box_size
      write(*,*) 'normal fluid velocity field is: ', trim(normal_velocity)
      select case(normal_velocity)
        case('xflow')
          write(*,'(a,f6.3)') ' u(x)=', vel_xflow
        case('ABC')
          write(*,'(a,f6.3,a,f6.3,a,f6.3)') ' A=', abc_A, ' B=', abc_B, ' C=', abc_C
        case('KS')
          call setup_KS !ksmodel.mod
        case('compressible')
          if (periodic_bc) then 
            write(*,'(a,i4.2,a)') 'creating gradient of a random field  on a', nfm_size,'^3 mesh'
            call setup_compressible !normal_fluid.mod
          else
            call fatal_error('normal_fluid.mod','comrpressible v needs periodic b.c.')
          end if
      end select
      write(*,'(a,f6.4,a,f6.4)') ' mutual friction coefficients alpha=',alpha(1),' alpha`=',alpha(2) 
      if (normal_fluid_cutoff<100.) then
        write(*,'(a,f6.3)') ' normal fluid is turned off when t=',normal_fluid_cutoff 
      end if
    end subroutine
    !************************************************************
    subroutine get_normal_velocity(x,u)
      !get the normal fluid at a position x
      implicit none
      real, intent(IN) :: x(3) !position of the particle
      real, intent(OUT) :: u(3) !velocity at x
      select case(normal_velocity)
        case('zero')
          u=0. ! no flow
        case('xflow')
          u=0. ; u(1)=vel_xflow !flow in the x direction
        case('ABC')
          !commonly used toy model - turbulence/dynamo
          u(1)=abc_B*cos(abc_k*x(2))+abc_C*sin(abc_k*x(3))
          u(2)=abc_C*cos(abc_k*x(3))+abc_A*sin(abc_k*x(1))
          u(3)=abc_A*cos(abc_k*x(1))+abc_B*sin(abc_k*x(2))
        case('KS')
          !multi-scale model of turbulence
          call get_KS_flow(x,u)
        case('compressible')
          !compressible flow on a mesh, needs interpolation
          call nfm_interpolation(x,u)
        case default
          call fatal_error('normal_fluid.mod:get_normal_fluid', &
          'correct parameter for normal_veloctity not set')
      end select
    end subroutine
    !**********************************************************
    subroutine nfm_interpolation(x,u)
      implicit none
      real, intent(IN) :: x(3) !position
      real, intent(OUT) :: u(3) !velocity
      real :: dx, dy, dz 
      integer :: i !for looping
      integer :: ix=1, iy=1, iz=1
      integer :: midpt
      midpt=nfm_size/2
      !find the closest meshpoints

      do i=1, nfm_size-1
        if ((x(1)>nfm(midpt,midpt,i)%x(1)).and.(x(1)<nfm(midpt,midpt,i+1)%x(1))) ix=i
        if ((x(2)>nfm(midpt,i,midpt)%x(2)).and.(x(2)<nfm(midpt,i+1,midpt)%x(2))) iy=i
        if ((x(3)>nfm(i,midpt,midpt)%x(3)).and.(x(3)<nfm(i+1,midpt,midpt)%x(3))) iz=i
      end do

      !find normalised disatances in all three spacial dimensions
      dx=abs(x(1)-nfm(iz,iy,ix)%x(1))/nfm_res
      dy=abs(x(2)-nfm(iz,iy,ix)%x(2))/nfm_res
      dz=abs(x(3)-nfm(iz,iy,ix)%x(3))/nfm_res

      !trilinear interpolation
      u=(1.-dx)*(1.-dy)*(1.-dz)*nfm(iz,iy,ix)%u+&
           dx*(1.-dy)*(1.-dz)*nfm(iz,iy,ix+1)%u+&
           (1.-dx)*dy*(1.-dz)*nfm(iz,iy+1,ix)%u+&
           (1.-dx)*(1.-dy)*dz*nfm(iz+1,iy,ix)%u+&
              dx*dy*(1.-dz)*nfm(iz,iy+1,ix+1)%u+&
              dx*(1.-dy)*dz*nfm(iz+1,iy,ix+1)%u+&
              (1.-dx)*dy*dz*nfm(iz+1,iy+1,ix)%u+&
                 dx*dy*dz*nfm(iz+1,iy+1,ix+1)%u
    end subroutine
    !**********************************************************
    subroutine setup_compressible
      implicit none
      integer :: i, j, k
      integer :: im, jm, km, ip, jp, kp
      !print dimensions to file for matlab
      open(unit=77,file='./data/nfm_dims.log',status='replace')
        write(77,*) nfm_size
      close(77)
      allocate(nfm(nfm_size,nfm_size,nfm_size))
      nfm_res=(real(box_size)/nfm_size)
      nfm_inv_res=1./nfm_res
      do k=1, nfm_size  ; do j=1, nfm_size ; do i=1, nfm_size
        nfm(k,j,i)%x(1)=nfm_res*real(2*i-1)/2.-(box_size/2.)
        nfm(k,j,i)%x(2)=nfm_res*real(2*j-1)/2.-(box_size/2.)
        nfm(k,j,i)%x(3)=nfm_res*real(2*k-1)/2.-(box_size/2.)
        !first set the scalar field
        call random_number(nfm(k,j,i)%v)
        nfm(k,j,i)%v=nfm(k,j,i)%v*2.-1. ; nfm(k,j,i)%v=0.01*nfm(k,j,i)%v
      end do ; end do ; end do
      ! now we must take the gradient of this
      do k=1, nfm_size  ; do j=1, nfm_size ; do i=1, nfm_size
        if (k==1) then !PERIODICITY
          kp=k+1 ; km=nfm_size
        else if (k==nfm_size) then
          kp=1 ; km=k-1
        else
          kp=k+1 ; km=k-1
        end if 
        if (j==1) then !PERIODICITY
          jp=j+1 ; jm=nfm_size
        else if (j==nfm_size) then
          jp=1 ; jm=j-1
        else
          jp=j+1 ; jm=j-1
        end if 
        if (i==1) then !PERIODICITY
          ip=i+1 ; im=nfm_size
        else if (i==nfm_size) then
          ip=1 ; im=i-1
        else
          ip=i+1 ; im=i-1
        end if
        !get the velocity field - gradient of scalar field
        nfm(k,j,i)%u(1)=(nfm(k,j,ip)%v-nfm(k,j,im)%v)*nfm_inv_res
        nfm(k,j,i)%u(2)=(nfm(k,jp,i)%v-nfm(k,jm,i)%v)*nfm_inv_res
        nfm(k,j,i)%u(3)=(nfm(kp,j,i)%v-nfm(km,j,i)%v)*nfm_inv_res
      end do ; end do ; end do
      ! calculate the divergence of this field
      do k=1, nfm_size  ; do j=1, nfm_size ; do i=1, nfm_size
        if (k==1) then !PERIODICITY
          kp=k+1 ; km=nfm_size
        else if (k==nfm_size) then
          kp=1 ; km=k-1
        else
          kp=k+1 ; km=k-1
        end if 
        if (j==1) then !PERIODICITY
          jp=j+1 ; jm=nfm_size
        else if (j==nfm_size) then
          jp=1 ; jm=j-1
        else
          jp=j+1 ; jm=j-1
        end if 
        if (i==1) then !PERIODICITY
          ip=i+1 ; im=nfm_size
        else if (i==nfm_size) then
          ip=1 ; im=i-1
        else
          ip=i+1 ; im=i-1
        end if
        !take the gradient of a random scalar field
        nfm(k,j,i)%div=(nfm(k,j,ip)%u(1)-nfm(k,j,im)%u(1))*nfm_inv_res + &
                                (nfm(k,jp,i)%u(2)-nfm(k,jm,i)%u(2))*nfm_inv_res + &
                                (nfm(kp,j,i)%u(3)-nfm(km,j,i)%u(3))*nfm_inv_res
      end do ; end do ; end do
      write(*,'(a)') 'velocity field calculated, printing to ./data/compressible_mesh.dat'
      open(unit=92,file='./data/compressible_mesh.dat',form='unformatted',status='replace',access='stream')
        write(92) nfm(nfm_size/2,nfm_size/2,1:nfm_size)%x(1)
        write(92) nfm(:,:,:)%v
        write(92) nfm(:,:,:)%div
        write(92) nfm(:,:,:)%u(1)
        write(92) nfm(:,:,:)%u(2)
        write(92) nfm(:,:,:)%u(3)
      close(92)
    end subroutine
end module
