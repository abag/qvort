module forcing
  !THIS MODULE CONTAINS ALL THE ROUTINES TO EMPLOY FORCING IN THE CODE
  use cdata
  type fgrid
    real :: x(3) !position
    real :: u(3) !velocity 
    real :: v !scalar
  end type
  real, private :: force_direction(3)=0.
  integer, private, parameter :: fmesh_size=128 !to be used with below
  type(fgrid), allocatable, private :: forcing_mesh(:,:,:) 
  contains
  subroutine setup_forcing()
    !essentially check all the necessary conditions to use forcing are set in run.in
    implicit none
    integer :: i, j, k
    integer :: im, jm, km, ip, jp, kp
    select case(force)
      case('off')
        write(*,*) 'no forcing employed'
      case('top_boundary')
        !check initial conditions
         select case(initf)
           case('single_line', 'line_motion')
             write(*,'(a,f6.3,a,f6.3)') 'forcing top boundary with amplitude', force_amp, &
             'frequency', force_freq
           case default
             call fatal_error('forcing.mod:setup_forcing', &
             'incorrect initf for forcing to be applied') !cdata.mod
         end select
      case('box_shake')
        write(*,'(a,f5.3,a,i4.3,a)') 'box shaking forcing with amplitude ', force_amp, &
        ' forcing direction changed every ', ceiling(force_freq), 'timesteps'
      case('delta_corr')
        write(*,'(a,f5.3)') 'delta correlated (time/space) forcing with amplitude ', force_amp
      case('compressible')
        write(*,'(a,f5.3)') 'creating gradient of a random field  in a box with amplitude ', force_amp
        allocate(forcing_mesh(fmesh_size,fmesh_size,fmesh_size))
        do k=1, fmesh_size  ; do j=1, fmesh_size ; do i=1, fmesh_size
          forcing_mesh(k,j,i)%x(1)=(real(box_size)/fmesh_size)*real(2*i-1)/2.-(box_size/2.)
          forcing_mesh(k,j,i)%x(2)=(real(box_size)/fmesh_size)*real(2*j-1)/2.-(box_size/2.)
          forcing_mesh(k,j,i)%x(3)=(real(box_size)/fmesh_size)*real(2*k-1)/2.-(box_size/2.)
          !set the velocity slots - for safety make u_sup=0.
          call random_number(forcing_mesh(k,j,i)%v)
          forcing_mesh(k,j,i)%v=forcing_mesh(k,j,i)%v*2.-1.
          !scale by amplitude
          forcing_mesh(k,j,i)%v=forcing_mesh(k,j,i)%v*force_amp
        end do ; end do ; end do
        ! now we must take the gradient of this
        do k=1, fmesh_size  ; do j=1, fmesh_size ; do i=1, fmesh_size
          if (k==1) then !PERIODICITY
            kp=k+1 ; km=fmesh_size
          else if (k==fmesh_size) then
            kp=1 ; km=k-1
          else
            kp=k+1 ; km=k-1
          end if 
          if (j==1) then !PERIODICITY
            jp=j+1 ; jm=fmesh_size
          else if (j==fmesh_size) then
            jp=1 ; jm=j-1
          else
            jp=j+1 ; jm=j-1
          end if 
          if (i==1) then !PERIODICITY
            ip=i+1 ; im=fmesh_size
          else if (i==fmesh_size) then
            ip=1 ; im=i-1
          else
            ip=i+1 ; im=i-1
          end if
          !take the gradient of a random scalar field
          forcing_mesh(k,j,i)%u(1)=(forcing_mesh(k,j,ip)%v-forcing_mesh(k,j,im)%v)*(fmesh_size/real(box_size))
          forcing_mesh(k,j,i)%u(2)=(forcing_mesh(k,jp,i)%v-forcing_mesh(k,jm,i)%v)*(fmesh_size/real(box_size))
          forcing_mesh(k,j,i)%u(3)=(forcing_mesh(kp,j,i)%v-forcing_mesh(km,j,i)%v)*(fmesh_size/real(box_size))
        end do ; end do ; end do
        write(*,'(a)') 'velocity field calculated, printing to ./data/compressible_mesh.dat'
        open(unit=92,file='./data/compressible_mesh.dat',form='unformatted',status='replace',access='stream')
          write(92) forcing_mesh(fmesh_size/2,fmesh_size/2,1:fmesh_size)%x(1)
          write(92) forcing_mesh(:,:,:)%v
          write(92) forcing_mesh(:,:,:)%u(1)
          write(92) forcing_mesh(:,:,:)%u(2)
          write(92) forcing_mesh(:,:,:)%u(3)
        close(92)
        stop
      case default
        call fatal_error('forcing.mod:setup_forcing', &
        'incorrect forcing parameter set in run.in') !cdata.mod
      end select
  end subroutine
  !***********************************************
  subroutine get_forcing(i,u)
    !force the particles - an additional velocity
    implicit none
    integer, intent(IN) :: i
    real, intent(OUT) :: u(3)
    select case(force)
      case('off')
        u=0.
      case('top_boundary')
        if ((f(i)%x(3)-box_size/2.)>-1.5*delta) then
          !particle is sufficiently close to top boundary to force
          u=0.
          !sinusoidal forcing in the x direction
          u(1)=delta*force_amp*sin(force_freq*t/(2*pi))
        end if
      case('box_shake')
        if (i==1) then
          !do we need to generate a new forcing direction?
          if (mod(itime,ceiling(force_freq))==0) then
            call random_number(force_direction)
            force_direction=force_direction*2.-1.
            !normalise
            force_direction=force_direction/sqrt(dot_product(force_direction,force_direction))
          end if
        end if
        u=force_direction*force_amp
      case('delta_corr')
        !at each time and position generate a new random vector
        call random_number(force_direction)
        force_direction=force_direction*2.-1.
        !normalise it
        force_direction=force_direction/sqrt(dot_product(force_direction,force_direction))
        !use this as the forcing
        u=force_direction*force_amp
      case('compressible')
        !at each time and position generate a new random vector
        call random_number(force_direction)
        force_direction=force_direction*2.-1.
        !normalise it
        force_direction=force_direction/sqrt(dot_product(force_direction,force_direction))
        !use this as the forcing
        u=force_direction*force_amp
    end select
  end subroutine
end module
