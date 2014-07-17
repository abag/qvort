module killing_sphere
  use cdata
  use general
  use diagnostic
  use timestep
  contains
  !******************************************************************
  !> remove loops which hit the boundaries
  subroutine setup_killing_sphere()
    implicit none
    !check boundary conditions
    select case(boundary)
        case('open')
          !all OK
        case default
          call fatal_error('killing_sphere.mod','wrong boundary - must be open')
    end select
    if (killing_radius>epsilon(0.)) then
      write(*,*) 'Employing a "killing sphere":'
      write(*,*) 'all loops at a distance greater than, ', killing_radius, ' will be removed'
      write(*,*) 'The energy lost due to these removed loops will be printed to file'
    else
      call fatal_error('killing_sphere.mod','killing_radius must be larger than 0.')
    end if
  end subroutine
  !******************************************************************
  !> remove loops which hit the boundaries
  subroutine enforce_killing_sphere()
    implicit none
    integer :: i, j
    real :: dist, u(3)
    real :: energy1, energy2
    do i=1, pcount
      if (f(i)%infront==0) cycle !empty particle
      dist=sqrt(f(i)%x(1)**2+f(i)%x(2)**2+f(i)%x(3)**2)
      if (i==1) then
        print*, dist
      end if
      if (dist>killing_radius) then
        !compute energy before hand
        !$omp parallel do private(i,u)
        do j=1, pcount
          if (f(j)%infront==0) cycle !check for 'empty' particles
          call calc_velocity(u,j)
          f(j)%u(:)=u(:) !store the velocity for time-step
        end do
        !$omp end parallel do
        call energy_info
        energy1=energy
        print*, energy1
        call ksphere_loop_remove(i)
        !now compute energy afterwards
        !$omp parallel do private(i,u)
        do j=1, pcount
          if (f(j)%infront==0) cycle !check for 'empty' particles
          call calc_velocity(u,j)
          f(j)%u(:)=u(:) !store the velocity for time-step
        end do
        !$omp end parallel do
        call energy_info
        energy2=energy
        open(unit=32,file='./data/killing_sphere_E_loss.log',position='append')
          write(32,*) energy1-energy2
        close(32)
      end if
    end do
  end subroutine
  !**************************************************
  !>remove loops that have left the box as they pass through the killing sphere
  !>this is very similar to loop_killer in line.f90
  !!a better way would be have a separate loop count
  !!and loop removal routine in general.mod
  subroutine ksphere_loop_remove(particle)
    implicit none
    integer :: particle, next
    integer :: store_next
    integer :: i, counter
    next=particle 
    boundary_loop_remove_count=boundary_loop_remove_count+1
    do i=1, pcount
      store_next=f(next)%infront
      if (store_next/=particle) then
        boundary_loop_remove_length=boundary_loop_remove_length&
                                   +distf(next,store_next)
      end if
      call clear_particle(next) !general.mod
      next=store_next
      if (next==particle) then
        exit  
      end if
    end do
  end subroutine
end module