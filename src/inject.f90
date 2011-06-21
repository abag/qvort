!>the routines contained within this module will inject loops into the code every nth timestep
!>(n is selected in run.in via inject_skip), the size of the loop is set using inject_size and the 
!>topology of the injected loop is set via loop_type, again in run.in.
module inject
  use cdata
  use general
  use periodic
  contains
  !>check all the necessary conditions to inject vortices are set in run.in
  subroutine setup_vortex_injection()
    implicit none
    !check inject size is not too small
    select case(inject_type)
      case('off')
        return !leave the routine
      case default
        if (inject_size<5) call fatal_error('inject.mod','inject size is too small')
        write(*,'(a)') ' ------------------VORTEX INJECTION-------------------' 
        write(*,'(a,i4.1,a,i3.1,a)') ' loops will be injected every ', inject_skip, ' timesteps with ', inject_size, ' points'
        write(*,'(a,a)') ' inject type is set to: ', trim(inject_type)
        select case(inject_type)
          case('rand-yz-loop2')
            write(*,'(a,f7.4,a)') ' rotation applied to loops: ', rotation_factor, '*2\pi'
        end select
        !check if we have set an injection stop time in run.in
        if (inject_stop<1E6) then
          write(*,'(a,f10.4)') ' inject will stop after t= ', inject_stop
        end if
    end select
  end subroutine
  !*******************************************************************
  !>the routine which injects the vortices into the code
  subroutine vortex_inject
    implicit none    
    type(qvort), allocatable, dimension(:) :: tmp
    real :: radius !used for loops
    real :: rand1, rand2 !random numbers
    real :: anglex,angley,anglez
    real,dimension(3)::dummy_xp_1, dummy_xp_2, dummy_xp_3, dummy_xp_4
    integer :: old_pcount
    integer :: i, j
    !check that inject_size is not 0
    if (inject_size==0) return !leave routine if it is
    !check if inject is set to off
    select case(inject_type)
      case('off')
        return !leave routine
    end select
    !first store the old point count
    old_pcount=pcount
    !now we must resize the array 
    !increment by inject_size using a dummy array tmp
    allocate(tmp(size(f)+inject_size)) ; tmp(:)%infront=0 !0 the infront array
    !copy accross information
    tmp(1:size(f)) = f
    !now deallocate tmp and transfer particles back to f
    !move_alloc is new intrinsic function (fortran 2003)
    call move_alloc(from=tmp,to=f)
    pcount=pcount+inject_size !increase the particle count
    !---------------now check what the inject_type is--------------
    select case(inject_type)
      case('xy-loop')
        radius=(0.75*inject_size*delta)/(2*pi) !75% of potential size
        !loop over particles setting spatial and 'loop' position
        do i=old_pcount+1, pcount
          f(i)%x(1)=radius*sin(pi*real(2*i-1)/inject_size)
          f(i)%x(2)=radius*cos(pi*real(2*i-1)/inject_size)
          f(i)%x(3)=-box_size/2.
          if (i==old_pcount+1) then
            f(i)%behind=pcount ; f(i)%infront=i+1
          else if (i==pcount) then 
            f(i)%behind=i-1 ; f(i)%infront=old_pcount+1
          else
            f(i)%behind=i-1 ; f(i)%infront=i+1
          end if
          !zero the stored velocities
          f(i)%u1=0. ; f(i)%u2=0.
        end do
      case('yz-loop')
        radius=(0.75*inject_size*delta)/(2*pi) !75% of potential size
        !loop over particles setting spatial and 'loop' position
        do i=old_pcount+1, pcount
          f(i)%x(1)=-box_size/2.
          f(i)%x(2)=radius*cos(pi*real(2*i-1)/inject_size)
          f(i)%x(3)=radius*sin(pi*real(2*i-1)/inject_size)
          if (i==old_pcount+1) then
            f(i)%behind=pcount ; f(i)%infront=i+1
          else if (i==pcount) then 
            f(i)%behind=i-1 ; f(i)%infront=old_pcount+1
          else
            f(i)%behind=i-1 ; f(i)%infront=i+1
          end if
          !zero the stored velocities
          f(i)%u1=0. ; f(i)%u2=0.
        end do
      case('rand-yz-loop')
        radius=(0.75*inject_size*delta)/(2*pi) !75% of potential size
        !loop over particles setting spatial and 'loop' position
        call random_number(rand1)
        call random_number(rand2)
        rand1=box_size*(rand1*2.-1.)/20.
        rand2=box_size*(rand2*2.-1.)/20.
        do i=old_pcount+1, pcount
          f(i)%x(1)=-box_size/2.
          f(i)%x(2)=radius*cos(pi*real(2*i-1)/inject_size)+rand1
          f(i)%x(3)=radius*sin(pi*real(2*i-1)/inject_size)+rand2
          if (i==old_pcount+1) then
            f(i)%behind=pcount ; f(i)%infront=i+1
          else if (i==pcount) then 
            f(i)%behind=i-1 ; f(i)%infront=old_pcount+1
          else
            f(i)%behind=i-1 ; f(i)%infront=i+1
          end if
          !zero the stored velocities
          f(i)%u1=0. ; f(i)%u2=0.
        end do
      case('rand-yz-loop2')
        radius=(0.75*inject_size*delta)/(2*pi) !75% of potential size
        !loop over particles setting spatial and 'loop' position
        call random_number(rand1)
        call random_number(rand2)
        call random_number(anglex)
        call random_number(angley)
        call random_number(anglez)
        anglex=(2.*anglex-1.)*2*pi*rotation_factor
        angley=(2.*angley-1.)*2*pi*rotation_factor
        anglez=(2.*anglez-1.)*2*pi*rotation_factor
        rand1=box_size*(rand1*2.-1.)/5.
        rand2=box_size*(rand2*2.-1.)/5.
        do i=old_pcount+1, pcount
          dummy_xp_1(1)=0.
          dummy_xp_1(2)=radius*cos(pi*real(2*i-1)/inject_size)
          dummy_xp_1(3)=radius*sin(pi*real(2*i-1)/inject_size)

          dummy_xp_2(1)=dummy_xp_1(1)
          dummy_xp_2(2)=dummy_xp_1(2)*cos(anglex)+dummy_xp_1(3)*sin(anglex)
          dummy_xp_2(3)=dummy_xp_1(2)*(-sin(anglex))+dummy_xp_1(3)*cos(anglex)

          dummy_xp_3(1)=dummy_xp_2(1)*cos(angley)+dummy_xp_2(3)*(-sin(angley))
          dummy_xp_3(2)=dummy_xp_2(2)
          dummy_xp_3(3)=dummy_xp_2(1)*sin(angley)+dummy_xp_2(3)*cos(angley)

          dummy_xp_4(1)=dummy_xp_3(1)*cos(anglez)+dummy_xp_3(2)*sin(anglez)
          dummy_xp_4(2)=dummy_xp_3(1)*(-sin(anglez))+dummy_xp_3(2)*cos(anglez)
          dummy_xp_4(3)=dummy_xp_3(3)
    
          f(i)%x(1)=dummy_xp_4(1)-box_size/2.2
          f(i)%x(2)=dummy_xp_4(2)+rand1
          f(i)%x(3)=dummy_xp_4(3)+rand2
          if (i==old_pcount+1) then
            f(i)%behind=pcount ; f(i)%infront=i+1
          else if (i==pcount) then 
            f(i)%behind=i-1 ; f(i)%infront=old_pcount+1
          else
            f(i)%behind=i-1 ; f(i)%infront=i+1
          end if
          !zero the stored velocities
          f(i)%u1=0. ; f(i)%u2=0.
        end do
      case('focus-yz-loop')
        radius=(0.75*inject_size*delta)/(2*pi) !75% of potential size
        !loop over particles setting spatial and 'loop' position
        call random_number(rand1)
        call random_number(rand2)
        call random_number(anglex)
        call random_number(angley)
        call random_number(anglez)
        anglex=(2.*anglex-1.)*2*pi/10.
        angley=(2.*angley-1.)*2*pi/10.
        anglez=(2.*anglez-1.)*2*pi/10.
        rand1=box_size*(rand1*2.-1.)/20.
        rand2=box_size*(rand2*2.-1.)/20.
        do i=old_pcount+1, pcount
          dummy_xp_1(1)=-box_size/2.
          dummy_xp_1(2)=radius*cos(pi*real(2*i-1)/inject_size)+rand1
          dummy_xp_1(3)=radius*sin(pi*real(2*i-1)/inject_size)+rand2

          dummy_xp_2(1)=dummy_xp_1(1)
          dummy_xp_2(2)=dummy_xp_1(2)*cos(anglex)+dummy_xp_1(3)*sin(anglex)
          dummy_xp_2(3)=dummy_xp_1(2)*(-sin(anglex))+dummy_xp_1(3)*cos(anglex)

          dummy_xp_3(1)=dummy_xp_2(1)*cos(angley)+dummy_xp_2(3)*(-sin(angley))
          dummy_xp_3(2)=dummy_xp_2(2)
          dummy_xp_3(3)=dummy_xp_2(1)*sin(angley)+dummy_xp_2(3)*cos(angley)

          dummy_xp_4(1)=dummy_xp_3(1)*cos(anglez)+dummy_xp_3(2)*sin(anglez)
          dummy_xp_4(2)=dummy_xp_3(1)*(-sin(anglez))+dummy_xp_3(2)*cos(anglez)
          dummy_xp_4(3)=dummy_xp_3(3)
    
          f(i)%x(1)=dummy_xp_4(1)
          f(i)%x(2)=dummy_xp_4(2)
          f(i)%x(3)=dummy_xp_4(3)
          if (i==old_pcount+1) then
            f(i)%behind=pcount ; f(i)%infront=i+1
          else if (i==pcount) then 
            f(i)%behind=i-1 ; f(i)%infront=old_pcount+1
          else
            f(i)%behind=i-1 ; f(i)%infront=i+1
          end if
          !zero the stored velocities
          f(i)%u1=0. ; f(i)%u2=0.
        end do
      case default
        call fatal_error('inject.mod','inject_type not set to useable value')    
    end select
  end subroutine
end module
