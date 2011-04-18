!>normal fluid component in the equation of motion, if the filament acts as a vortex
!>velocity enters as \f$\mathbf{u}=\alpha \mathbf{s}' \times (\mathbf{u}_\mathrm{n}-\mathbf{u}_\mathrm{s})+
!>\alpha' \mathbf{s}' \times [\mathbf{s}' \times (\mathbf{u}_\mathrm{n}-\mathbf{u}_\mathrm{s})]\f$
!>if the filament is a magnetic flux tube or a material line then the normaly velocity is used directly.
!>for further information on the normal fluid options see \ref NF
module normal_fluid
  use cdata
  use general
  use ksmodel
    !parameters used in get_normal_fluid
    real, parameter, private :: vel_xflow=1.
    real, parameter, private :: abc_A=1., abc_B=1., abc_C=1.
    real, private :: norm_k
    real, private :: urms_norm
    !>normal fluid mesh - all grid based normal velocities use this (e.g. Navier Stokes future)
    !>@param x the position on the grid
    !>@param u the velocity
    !>@param v a scalar used in the compressible flow case
    !>@param div the divergence of the velocity field at x
    type norm_fluid_grid
      private
      real :: x(3)
      real :: u(3)
      real :: v
      real :: div 
    end type
    !use the abbreviation nfm (normal fluid mesh)
    !>the size of the normal fluid mesh, put in run.in soon
    integer, private, parameter :: nfm_size=64 
    real, private :: nfm_res, nfm_inv_res !resolution/inv resolution of normal fluid mesh
    !>the normal fluid mesh
    type(norm_fluid_grid), allocatable, private :: nfm(:,:,:)
    contains 
    !************************************************************
    !>setup everything needed to use a normal fluid
    subroutine setup_normal_fluid
      !in here we print to file the timescale of the flow
      implicit none
      norm_k=2.*pi/box_size !wavenumber used in a number of normal fluid
      write(*,*) 'normal fluid velocity field is: ', trim(normal_velocity)
      select case(normal_velocity)
        case('xflow')
          urms_norm=vel_xflow
          write(*,'(a,f6.3)') ' u(x)=', vel_xflow
          open(unit=77,file='./data/normal_timescale.log',status='replace')
            write(77,*) box_size/vel_xflow !time-taken to cross box
          close(77)
        case('ABC')
          write(*,'(a,f6.3,a,f6.3,a,f6.3)') ' A=', abc_A, ' B=', abc_B, ' C=', abc_C
          call setup_gen_normalf !normal_fluid.mod
          open(unit=77,file='./data/normal_timescale.log',status='replace')
            write(77,*) box_size/urms_norm !size of box scaled by urms
          close(77)
        case('taylor-green')
          call setup_gen_normalf !normal_fluid.mod
          open(unit=77,file='./data/normal_timescale.log',status='replace')
            write(77,*) box_size/urms_norm !size of box scaled by urms
          close(77)
        case('shear')
          call setup_gen_normalf !normal_fluid.mod
          open(unit=77,file='./data/normal_timescale.log',status='replace')
            write(77,*) box_size/urms_norm !size of box scaled by urms
          close(77)
        case('sod')
          call setup_gen_normalf !normal_fluid.mod
          open(unit=77,file='./data/normal_timescale.log',status='replace')
            write(77,*) box_size/urms_norm !size of box scaled by urms
          close(77)
        case('galloway-proctor')
          call setup_gen_normalf !normal_fluid.mod
          open(unit=77,file='./data/normal_timescale.log',status='replace')
            write(77,*) box_size/urms_norm !size of box scaled by urms
          close(77)
        case('KS')
          call setup_KS !ksmodel.mod
          call print_KS_Mesh !normal_fluid.mod
          !Normal vel timescale written in setup_KS (line 109)
        case('compressible')
          if (periodic_bc) then 
            write(*,'(a,i4.2,a)') ' creating gradient of a random field  on a', nfm_size,'^3 mesh'
            call setup_compressible !normal_fluid.mod
            open(unit=77,file='./data/normal_timescale.log',status='replace')
              write(77,*) box_size/urms_norm !size of box scaled by urms
            close(77)
          else
            call fatal_error('normal_fluid.mod','comrpressible v needs periodic b.c.')
          end if
      end select
      write(*,'(a,f6.4,a,f6.4)') ' mutual friction coefficients alpha=',alpha(1),' alpha`=',alpha(2) 
      if (normal_fluid_cutoff<100.) then
        write(*,'(a,f6.3)') ' normal fluid is turned off when t=',normal_fluid_cutoff 
      end if
      write(*,'(a,f8.4)') ' normal fluid rms velocity is: ', urms_norm
      write(*,'(a)') ' normal fluid timescale printed to ./data/normal_timescale.log'
    end subroutine
    !************************************************************
    !>get the velocity (u) of the normal fluid at a position (x)
    !>see the normal fluid page for more information
    subroutine get_normal_velocity(x,u)
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
          u(1)=abc_B*cos(norm_k*x(2))+abc_C*sin(norm_k*x(3))
          u(2)=abc_C*cos(norm_k*x(3))+abc_A*sin(norm_k*x(1))
          u(3)=abc_A*cos(norm_k*x(1))+abc_B*sin(norm_k*x(2))
        case('taylor-green')
          !The Taylor-Green Vortex
          u(1)=sin(norm_k*x(1))*cos(norm_k*x(2))*cos(norm_k*x(3))
          u(2)=-cos(norm_k*x(1))*sin(norm_k*x(2))*cos(norm_k*x(3))
          u(3)=0.
        case('shear')
          u(1)=exp(-(6*x(3)/box_size)**2)
        case('sod')
          if (x(1)<box_size/3.) then
            u(1)=0.5*tanh(10.*(x(1)+box_size/5.)/box_size)+0.5
          else 
            u(1)=0.
          end if
        case('galloway-proctor')
          u(1)=-sin(norm_k*x(2)+cos(t))
          u(2)=-cos(norm_k*x(1)+sin(t))
          u(3)=sin(norm_k*x(1)+sin(t))+cos(norm_k*x(2)+cos(t))
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
    !>(linear) interpolation of the normal fluid velocity from a mesh
    !>needed for certain selections of velocity
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
    !>setup a compressible velocity field by taking the gradient 
    !>of a random scalar field
    subroutine setup_compressible
      implicit none
      real :: u_rms
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
      u_rms=0. !0 the root mean squared velocity
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
        u_rms=u_rms+(nfm(k,j,i)%u(1)**2+nfm(k,j,i)%u(2)**2+nfm(k,j,i)%u(3)**2)
      end do ; end do ; end do
      u_rms=sqrt(u_rms/(nfm_size**3))
      !now we normalise this field by u_rms
      nfm%u(1)=nfm%u(1)/u_rms ; nfm%u(2)=nfm%u(2)/u_rms ; nfm%u(3)=nfm%u(3)/u_rms
      urms_norm=1.!urms_norm must be 1 as we have scaled it to be!
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
      write(*,'(a)') ' velocity field calculated, printing to ./data/compressible_mesh.dat'
      open(unit=92,file='./data/compressible_mesh.dat',form='unformatted',status='replace',access='stream')
        write(92) nfm(nfm_size/2,nfm_size/2,1:nfm_size)%x(1)
        write(92) nfm(:,:,:)%v
        write(92) nfm(:,:,:)%div
        write(92) nfm(:,:,:)%u(1)
        write(92) nfm(:,:,:)%u(2)
        write(92) nfm(:,:,:)%u(3)
      close(92)
    end subroutine 
    !**********************************************************
    !>set up simple normal fluid models (ones with an analytic form)
    !>and print this field (discretized on a cubic mesh) to a binary file
    subroutine setup_gen_normalf
      implicit none
      integer :: i, j, k
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
      end do ; end do ; end do
      urms_norm=0. !0 the root mean squared velocity
      do k=1, nfm_size  ; do j=1, nfm_size ; do i=1, nfm_size
        !get the velocity field - shearing wave
        call get_normal_velocity(nfm(k,j,i)%x,nfm(k,j,i)%u)
        urms_norm=urms_norm+(nfm(k,j,i)%u(1)**2+nfm(k,j,i)%u(2)**2+nfm(k,j,i)%u(3)**2)
      end do ; end do ; end do
      urms_norm=sqrt(urms_norm/(nfm_size**3))
      write(*,'(a)') ' velocity field calculated, printing to ./data/norm_init_mesh.dat'
      open(unit=92,file='./data/norm_init_mesh.dat',form='unformatted',status='replace',access='stream')
        write(92) nfm(nfm_size/2,nfm_size/2,1:nfm_size)%x(1)
        write(92) nfm(:,:,:)%u(1)
        write(92) nfm(:,:,:)%u(2)
        write(92) nfm(:,:,:)%u(3)
      close(92)
    end subroutine
    !**********************************************************
    !>print the KS velocity field (on a cubic mesh) to a binary file
    subroutine print_KS_Mesh
      implicit none
      integer :: i, j, k
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
      end do ; end do ; end do
      urms_norm=0. !0 the root mean squared velocity
      do k=1, nfm_size  ; do j=1, nfm_size ; do i=1, nfm_size
        !get the velocity field - ABC flow
        call get_KS_flow(nfm(k,j,i)%x,nfm(k,j,i)%u)       
        urms_norm=urms_norm+(nfm(k,j,i)%u(1)**2+nfm(k,j,i)%u(2)**2+nfm(k,j,i)%u(3)**2)
      end do ; end do ; end do
      urms_norm=sqrt(urms_norm/(nfm_size**3))
      write(*,'(a)') ' velocity field calculated, printing to ./data/KS_mesh.dat'
      open(unit=92,file='./data/KS_mesh.dat',form='unformatted',status='replace',access='stream')
        write(92) nfm(nfm_size/2,nfm_size/2,1:nfm_size)%x(1)
        write(92) nfm(:,:,:)%u(1)
        write(92) nfm(:,:,:)%u(2)
        write(92) nfm(:,:,:)%u(3)
      close(92)
    end subroutine
end module
!>\page NF Normal fluid velocity field
!!Normal fluid velocity field is set in run.in throught the parameter normal_velocity\n
!!Options are:\n
!!- \p zero - no flow, friction only
!!- \p ABC - \f$\mathbf{u}=(\cos y+\sin z,\sin x+\cos z, \cos x+\sin y)\f$ \n
!!- \p xflow - \f$\mathbf{u}=(u_x,0,0)\f$ \n
!!- \p taylor-green \n
!!- \p KS - multi-scale model of turbulence see KSmodel
!!- \p shear - \f$\mathbf{u}=(0,0,u_0\exp(-z^2))\f$ \n
!!- \p galloway-proctor  \f$\mathbf{u}=(\sin(y+\cos t)),-\cos(x+\sin t),\sin(x+\sin(t))+\cos(y+\cos(t)))\f$ \n
!!- \p compressible - a compressible flow \f$\mathbf{u}=\nabla\phi\f$, where \f$\phi\f$ is a random scalar field
