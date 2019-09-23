module killing_sphere
  use cdata
  use general
  use diagnostic
  use timestep
  implicit none
  real,private :: centre_of_killing(3)=0.
  contains
  !******************************************************************
  !> initialise the killing sphere
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
      write(*,'(a,f6.4,a)') 'all loops at a distance greater than, ', killing_radius, ' will be removed'
      write(*,*) 'The energy lost due to these removed loops will be printed to file'
    else
      call fatal_error('killing_sphere.mod','killing_radius must be larger than 0.')
    end if
    if (adaptive_killing_sphere) then
      write(*,*) 'killing sphere will adapt to centre of mass and std of tangle'
    end if
    if (moving_box) then
      write(*,*) 'box will move with the centre of mass of the tangle'
    end if
  end subroutine
    !******************************************************************
  !> change to a frame moving with the tangle
  subroutine killing_moving_frame
    implicit none
    !check boundary conditions
    real:: sum_x(3), count_x
    !first compute the centre of "mass"
    sum_x(1)=sum(f(:)%x(1), mask=f(:)%infront>0)
    sum_x(2)=sum(f(:)%x(2), mask=f(:)%infront>0)
    sum_x(3)=sum(f(:)%x(3), mask=f(:)%infront>0)
    count_x=count(mask=f(:)%infront>0)
    !set this to be the centre of killing
    centre_of_killing=sum_x/count_x
    f(:)%x(1)=f(:)%x(1)-centre_of_killing(1)
    f(:)%x(2)=f(:)%x(2)-centre_of_killing(2)
    f(:)%x(3)=f(:)%x(3)-centre_of_killing(3)
  end subroutine
  !******************************************************************
  !> remove loops which hit the sphere
  subroutine enforce_killing_sphere()
    implicit none
    integer :: i, j
    real :: dist, u(3)
    real :: energy1, energy2
    !do we adapt the killing sphere?
    if (adaptive_killing_sphere) call adapt_killing_sphere !killing_sphere.mod
    do i=1, pcount
      if (f(i)%infront==0) cycle !empty particle
      dist=sqrt((f(i)%x(1)-centre_of_killing(1))**2+&
                (f(i)%x(2)-centre_of_killing(2))**2+&
                (f(i)%x(3)-centre_of_killing(3))**2)
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
  !******************************************************************
  !> adapt the killing sphere
  subroutine adapt_killing_sphere()
    implicit none
    real:: sum_x(3), count_x,x_std(3), sum_x2(3)
    !first compute the centre of "mass"
    sum_x(1)=sum(f(:)%x(1), mask=f(:)%infront>0)
    sum_x(2)=sum(f(:)%x(2), mask=f(:)%infront>0)
    sum_x(3)=sum(f(:)%x(3), mask=f(:)%infront>0)
    count_x=count(mask=f(:)%infront>0)
    !set this to be the centre of killing
    centre_of_killing=sum_x/count_x
    !now the standard deviation
    sum_x2(1)=sum((f(:)%x(1)-centre_of_killing(1))**2, mask=f(:)%infront>0)
    sum_x2(2)=sum((f(:)%x(2)-centre_of_killing(2))**2, mask=f(:)%infront>0)
    sum_x2(3)=sum((f(:)%x(3)-centre_of_killing(3))**2, mask=f(:)%infront>0)
    x_std(1)=sqrt(sum_x2(1)/(count_x-1))
    x_std(2)=sqrt(sum_x2(2)/(count_x-1))
    x_std(3)=sqrt(sum_x2(3)/(count_x-1))
    !now we set the killing radius to be 3 times the maximum standard deviation
    !as is standard in finding outliers
    killing_radius=3*maxval(x_std)
    if (mod(itime,shots)==0) then
      open(unit=37,file='./data/adaptive_killing_radius.log', position='append')
        write(37,*) centre_of_killing, killing_radius
      close(37)
    end if
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
