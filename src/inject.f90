!>the routines contained within this module will inject loops into the code every nth timestep
!>(n is selected in run.in via inject_skip). This is useful to model specific
!>experiments.
!>The size of the loop is set using inject_size and the 
!>topology of the injected loop is set via loop_type, again in run.in.
module inject
  use cdata
  use general
  use periodic
  contains
  !>check all the necessary conditions to inject vortices are set in run.in
  subroutine setup_vortex_injection()
    implicit none
    !check requirements of specific conditions
    select case(inject_type)
      case('off')
        return !leave the routine
      case('edge_pulse')
        if (mod(inject_size,6)/=0) then
          call fatal_error('inject.mod','inject size must be a multiple of 6')
        end if
    end select
    !check inject size is not too small
    if (inject_size<5) call fatal_error('inject.mod','inject size is too small')
    !print to screen pertinent details
    write(*,'(a)') ' ------------------VORTEX INJECTION-------------------' 
    write(*,'(a,i4.1,a,i3.1,a)') ' loops will be injected every ', inject_skip, ' timesteps with ', inject_size, ' points'
    write(*,'(a,a)') ' inject type is set to: ', trim(inject_type)
    select case(inject_type)
      case('rand-yz-loop-rot','rand-xyz-loop')
        write(*,'(a,f7.4,a)') ' rotation applied to loops: ', rotation_factor, '*2\pi'
    end select
    !check if we have set an injection stop time in run.in
    if (inject_stop<1E6) then
      write(*,'(a,f10.4)') ' inject will stop after t= ', inject_stop
    end if
  end subroutine
  !*******************************************************************
  !>the routine which injects the vortices into the code
  subroutine vortex_inject
    implicit none    
    type(qvort), allocatable, dimension(:) :: tmp !used to copy array to
    real :: radius !used for loops
    real :: rand1, rand2, rand3 !random numbers
    real :: anglex,angley,anglez !for rotating
    real,dimension(3)::dummy_xp_1, dummy_xp_2, dummy_xp_3, dummy_xp_4
    integer :: old_pcount
    integer :: i, j !to loop
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
      case('edge_pulse')!six rings at each face of the computational box
        radius=(0.75*inject_size*delta)/(6*2*pi) !75% of potential size 
        !loop over particles (6 loops) setting spatial and 'loop' position
        !first loop - xy plane top of box
        do i=old_pcount+1, old_pcount+inject_size/6
          f(i)%x(1)=radius*sin(pi*real(2*i-1)/(inject_size/6))
          f(i)%x(2)=radius*cos(pi*real(2*i-1)/(inject_size/6))
          f(i)%x(3)=box_size/2.01
          if (i==old_pcount+1) then
            f(i)%behind=old_pcount+inject_size/6 ; f(i)%infront=i+1
          else if (i==old_pcount+inject_size/6) then 
            f(i)%behind=i-1 ; f(i)%infront=old_pcount+1
          else
            f(i)%behind=i-1 ; f(i)%infront=i+1
          end if
          !zero the stored velocities
          f(i)%u1=0. ; f(i)%u2=0. ; f(i)%u3=0.
        end do
        !second loop- xy plane bottom of box
        do i=old_pcount+inject_size/6+1, old_pcount+inject_size/3
          f(i)%x(1)=-radius*sin(pi*real(2*i-1)/(inject_size/6))
          f(i)%x(2)=radius*cos(pi*real(2*i-1)/(inject_size/6))
          f(i)%x(3)=-box_size/2.01
          if (i==(old_pcount+inject_size/6+1)) then
            f(i)%behind=old_pcount+inject_size/3 ; f(i)%infront=i+1
          else if (i==old_pcount+inject_size/3) then 
            f(i)%behind=i-1 ; f(i)%infront=old_pcount+inject_size/6+1
          else
            f(i)%behind=i-1 ; f(i)%infront=i+1
          end if
          !zero the stored velocities
          f(i)%u1=0. ; f(i)%u2=0. ; f(i)%u3=0.
        end do    
        !third loop- yz plane pos side of box 
        do i=old_pcount+inject_size/3+1, old_pcount+inject_size/2
          f(i)%x(1)=box_size/2.01
          f(i)%x(2)=radius*cos(pi*real(2*i-1)/(inject_size/6))
          f(i)%x(3)=-radius*sin(pi*real(2*i-1)/(inject_size/6))
          if (i==(old_pcount+inject_size/3+1)) then
            f(i)%behind=old_pcount+inject_size/2 ; f(i)%infront=i+1
          else if (i==old_pcount+inject_size/2) then 
            f(i)%behind=i-1 ; f(i)%infront=old_pcount+inject_size/3+1
          else
            f(i)%behind=i-1 ; f(i)%infront=i+1
          end if
          !zero the stored velocities
          f(i)%u1=0. ; f(i)%u2=0. ; f(i)%u3=0.
        end do
        !fourth loop- yz plane neg side of box 
        do i=old_pcount+inject_size/2+1, old_pcount+2*inject_size/3
          f(i)%x(1)=-box_size/2.01
          f(i)%x(2)=radius*cos(pi*real(2*i-1)/(inject_size/6))
          f(i)%x(3)=radius*sin(pi*real(2*i-1)/(inject_size/6))
          if (i==(old_pcount+inject_size/2+1)) then
            f(i)%behind=old_pcount+2*inject_size/3 ; f(i)%infront=i+1
          else if (i==old_pcount+2*inject_size/3) then 
            f(i)%behind=i-1 ; f(i)%infront=old_pcount+inject_size/2+1
          else
            f(i)%behind=i-1 ; f(i)%infront=i+1
          end if
          !zero the stored velocities
          f(i)%u1=0. ; f(i)%u2=0. ; f(i)%u3=0.
        end do
        !fifth loop- xz plane pos side box
        do i=old_pcount+2*inject_size/3+1, old_pcount+5*inject_size/6
          f(i)%x(1)=-radius*sin(pi*real(2*i-1)/(inject_size/6))
          f(i)%x(2)=box_size/2.01
          f(i)%x(3)=radius*cos(pi*real(2*i-1)/(inject_size/6))
          if (i==(old_pcount+2*inject_size/3+1)) then
            f(i)%behind=old_pcount+5*inject_size/6 ; f(i)%infront=i+1
          else if (i==old_pcount+5*inject_size/6) then 
            f(i)%behind=i-1 ; f(i)%infront=old_pcount+2*inject_size/3+1
          else
            f(i)%behind=i-1 ; f(i)%infront=i+1
          end if
          !zero the stored velocities
          f(i)%u1=0. ; f(i)%u2=0. ; f(i)%u3=0.
        end do
        !final loop- xz plane neg side of box
        do i=old_pcount+5*inject_size/6+1, pcount
          f(i)%x(1)=radius*sin(pi*real(2*i-1)/(inject_size/6))
          f(i)%x(2)=-box_size/2.01
          f(i)%x(3)=radius*cos(pi*real(2*i-1)/(inject_size/6))
          if (i==(old_pcount+5*inject_size/6+1)) then
            f(i)%behind=pcount ; f(i)%infront=i+1
          else if (i==pcount) then 
            f(i)%behind=i-1 ; f(i)%infront=old_pcount+5*inject_size/6+1
          else
            f(i)%behind=i-1 ; f(i)%infront=i+1
          end if
          !zero the stored velocities
          f(i)%u1=0. ; f(i)%u2=0. ; f(i)%u3=0.
        end do
      case('xy-loop')!loops in xy-plane injected at bottom of box
        radius=(0.75*inject_size*delta)/(2*pi) !75% of potential size
        !loop over particles setting spatial and 'loop' position
        do i=old_pcount+1, pcount
          f(i)%x(1)=radius*sin(pi*real(2*i-1)/inject_size)
          f(i)%x(2)=radius*cos(pi*real(2*i-1)/inject_size)
          f(i)%x(3)=-box_size/2.1
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
      case('yz-loop') !loops in yz plane injected at side of box
        radius=(0.75*inject_size*delta)/(2*pi) !75% of potential size
        !loop over particles setting spatial and 'loop' position
        do i=old_pcount+1, pcount
          f(i)%x(1)=-box_size/2.1
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
      case('rand-yz-loop') !as yz-loop but with a slight pertubation in positions
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
      case('rand-yz-loop-rot') !pertubation in both position and angle
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
      case('rand-xyz-loop') !completetly random loops
        radius=(0.75*inject_size*delta)/(2*pi) !75% of potential size
        !loop over particles setting spatial and 'loop' position
        call random_number(rand1)
        call random_number(rand2)
        call random_number(rand3)
        call random_number(anglex)
        call random_number(angley)
        call random_number(anglez)
        anglex=(2.*anglex-1.)*2*pi*rotation_factor
        angley=(2.*angley-1.)*2*pi*rotation_factor
        anglez=(2.*anglez-1.)*2*pi*rotation_factor
        rand1=box_size*(rand1*2.-1.)/2.
        rand2=box_size*(rand2*2.-1.)/2.
        rand3=box_size*(rand3*2.-1.)/2.
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
    
          f(i)%x(1)=dummy_xp_4(1)+rand1
          f(i)%x(2)=dummy_xp_4(2)+rand2
          f(i)%x(3)=dummy_xp_4(3)+rand3
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
